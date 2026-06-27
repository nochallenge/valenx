//! A small **autonomy simulation harness**: a kinematic vehicle carrying a set
//! of sensors, advanced one `step(dt)` at a time, emitting a [`SensorFrame`]
//! (every sensor reading at time `t`) — the step API an autonomy stack or an RL
//! environment polls.
//!
//! ## The vehicle
//!
//! [`VehicleState`] is a deliberately minimal rigid-body kinematic state in a
//! local ENU world frame (east `+x`, north `+y`, up `+z`): a position, a
//! velocity, an orientation, and a body angular rate. It is **not** the
//! force-based `valenx-vehicle::Car` (a point-mass that computes 0–100 / top
//! speed / braking from power and grip): that crate answers *performance*
//! questions and has no time-stepped pose to attach sensors to, whereas an
//! autonomy harness needs exactly that pose. So the harness composes the
//! sensor models with a small kinematic state of its own. A control input
//! [`Command`] sets body-frame acceleration and angular rate each step; the
//! state integrates them (semi-implicit Euler) and the **world acceleration**
//! needed by the IMU is taken directly from the commanded body acceleration
//! rotated into the world, so the harness ↔ IMU kinematics are self-consistent.
//!
//! ## The frame the sensors see
//!
//! - The **GPS** reads the vehicle's ENU position directly (the world frame *is*
//!   ENU about the GPS datum).
//! - The **LiDAR** and **IMU** use the vehicle's world position and orientation.
//!   The harness assumes each sensor is mounted at the body origin aligned with
//!   the body axes; a fixed mount offset/rotation is a documented extension.
//!
//! Honest scope: the harness is a clean, deterministic *kinematic* testbed for
//! wiring and exercising perception/estimation/control loops. It is **not** a
//! dynamics simulator — there is no mass, tire, or contact model (use
//! `valenx-vehicle` / `valenx-mbd` for forces), and the world is the analytic
//! LiDAR [`Scene`]. Everything is reproducible: each sensor owns a seeded RNG.

use nalgebra::{UnitQuaternion, Vector3};

use crate::camera::Camera;
use crate::error::SensorError;
use crate::gps::{Gps, GpsFix};
use crate::imu::{BodyState, Imu, ImuReading};
use crate::lidar::{Lidar, LidarScan};
use crate::scene::Scene;

/// Minimal rigid-body kinematic state in the local ENU world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VehicleState {
    /// Position in the ENU world frame (m).
    pub position: Vector3<f64>,
    /// Velocity in the ENU world frame (m/s).
    pub velocity: Vector3<f64>,
    /// Orientation (body → world).
    pub orientation: UnitQuaternion<f64>,
    /// Body-frame angular rate (rad/s).
    pub angular_rate: Vector3<f64>,
}

impl Default for VehicleState {
    /// At the origin, at rest, level.
    fn default() -> Self {
        Self {
            position: Vector3::zeros(),
            velocity: Vector3::zeros(),
            orientation: UnitQuaternion::identity(),
            angular_rate: Vector3::zeros(),
        }
    }
}

/// A control input applied over one step: body-frame linear acceleration and
/// body-frame angular rate.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Command {
    /// Commanded linear acceleration in the **body** frame (m/s²).
    pub accel_body: Vector3<f64>,
    /// Commanded angular rate in the **body** frame (rad/s).
    pub angular_rate_body: Vector3<f64>,
}

impl Command {
    /// Coast: no acceleration, no rotation (holds velocity and heading).
    #[must_use]
    pub fn coast() -> Self {
        Self::default()
    }
}

impl VehicleState {
    /// Advance the state by `dt` seconds under `command` (semi-implicit Euler),
    /// returning the **world-frame** linear acceleration applied (the quantity
    /// the IMU measures specific force from).
    ///
    /// Orientation is updated by integrating the body angular rate as a rotation
    /// vector `ω·dt` and composing it on the right (a body-frame increment).
    fn integrate(&mut self, command: &Command, dt: f64) -> Vector3<f64> {
        // World-frame acceleration = R · a_body (R = body→world).
        let accel_world = self.orientation * command.accel_body;
        // Semi-implicit Euler: velocity first, then position with the new v.
        self.velocity += accel_world * dt;
        self.position += self.velocity * dt;
        // Orientation: right-compose the body-frame rotation increment.
        self.angular_rate = command.angular_rate_body;
        let dtheta = command.angular_rate_body * dt;
        let increment = UnitQuaternion::from_scaled_axis(dtheta);
        self.orientation *= increment;
        accel_world
    }
}

/// One synchronised set of sensor readings at simulation time `t`.
///
/// Each field is `Some` only if the harness has the corresponding sensor
/// attached.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorFrame {
    /// Simulation time of this frame (s).
    pub time: f64,
    /// The vehicle state at `time` (ground truth — useful for RL/V&V).
    pub state: VehicleState,
    /// LiDAR scan, if a LiDAR is attached.
    pub lidar: Option<LidarScan>,
    /// IMU reading, if an IMU is attached.
    pub imu: Option<ImuReading>,
    /// GPS fix, if a GPS is attached.
    pub gps: Option<GpsFix>,
}

/// The autonomy harness: a vehicle, an analytic world, attached sensors, and a
/// simulation clock.
#[derive(Debug, Clone)]
pub struct Harness {
    /// The current vehicle state.
    pub state: VehicleState,
    /// The analytic world the LiDAR ranges against.
    pub scene: Scene,
    /// Simulation time (s).
    pub time: f64,
    lidar: Option<Lidar>,
    imu: Option<Imu>,
    gps: Option<Gps>,
    /// World acceleration applied on the most recent step (for the IMU).
    last_accel_world: Vector3<f64>,
    /// A camera is carried for static projection queries (it has no time state);
    /// it is exposed via [`Harness::camera`] rather than sampled into a frame,
    /// since rendering an image needs a scene-to-pixel pipeline beyond the scope
    /// of this kinematic harness.
    camera: Option<Camera>,
}

impl Harness {
    /// A new harness with an initial state and a (possibly empty) scene, no
    /// sensors attached yet.
    #[must_use]
    pub fn new(state: VehicleState, scene: Scene) -> Self {
        Self {
            state,
            scene,
            time: 0.0,
            lidar: None,
            imu: None,
            gps: None,
            last_accel_world: Vector3::zeros(),
            camera: None,
        }
    }

    /// Attach a LiDAR (builder style).
    #[must_use]
    pub fn with_lidar(mut self, lidar: Lidar) -> Self {
        self.lidar = Some(lidar);
        self
    }

    /// Attach an IMU (builder style).
    #[must_use]
    pub fn with_imu(mut self, imu: Imu) -> Self {
        self.imu = Some(imu);
        self
    }

    /// Attach a GPS (builder style).
    #[must_use]
    pub fn with_gps(mut self, gps: Gps) -> Self {
        self.gps = Some(gps);
        self
    }

    /// Attach a camera (builder style). The camera is available via
    /// [`Harness::camera`] for projecting points; it is not sampled into a
    /// [`SensorFrame`].
    #[must_use]
    pub fn with_camera(mut self, camera: Camera) -> Self {
        self.camera = Some(camera);
        self
    }

    /// The attached camera, if any.
    #[must_use]
    pub fn camera(&self) -> Option<&Camera> {
        self.camera.as_ref()
    }

    /// Advance the simulation by `dt` seconds under `command`, then sample all
    /// attached sensors at the new time and return the [`SensorFrame`].
    ///
    /// # Errors
    /// Returns [`SensorError::InvalidConfig`] if `dt` is not finite and
    /// positive, or [`SensorError::NonFinite`] if the command has a non-finite
    /// component (so a bad control input fails loud rather than poisoning the
    /// state with `NaN`).
    pub fn step(&mut self, command: &Command, dt: f64) -> Result<SensorFrame, SensorError> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "dt must be finite and > 0, got {dt}"
            )));
        }
        if !(command.accel_body.iter().all(|c| c.is_finite())
            && command.angular_rate_body.iter().all(|c| c.is_finite()))
        {
            return Err(SensorError::NonFinite("command".into()));
        }

        self.last_accel_world = self.state.integrate(command, dt);
        self.time += dt;
        Ok(self.sample())
    }

    /// Sample every attached sensor at the current state/time without advancing
    /// (uses the world acceleration from the most recent [`Harness::step`], or
    /// zero if none has run yet). Useful for an initial `t = 0` reading.
    pub fn sample(&mut self) -> SensorFrame {
        let lidar = self
            .lidar
            .as_mut()
            .map(|l| l.scan(&self.scene, self.state.position, self.state.orientation));
        let imu = self.imu.as_mut().map(|i| {
            i.sample(&BodyState {
                orientation: self.state.orientation,
                accel_world: self.last_accel_world,
                angular_rate_body: self.state.angular_rate,
            })
        });
        let gps = self.gps.as_mut().map(|g| g.sample(self.state.position));

        SensorFrame {
            time: self.time,
            state: self.state,
            lidar,
            imu,
            gps,
        }
    }

    /// Run `n` steps under a constant `command`, collecting the frames.
    ///
    /// # Errors
    /// Propagates any [`Harness::step`] error (a bad `dt` or command).
    pub fn run(
        &mut self,
        command: &Command,
        dt: f64,
        n: usize,
    ) -> Result<Vec<SensorFrame>, SensorError> {
        let mut frames = Vec::with_capacity(n);
        for _ in 0..n {
            frames.push(self.step(command, dt)?);
        }
        Ok(frames)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::{Distortion, Intrinsics};
    use crate::gps::Geodetic;
    use crate::imu::ImuConfig;
    use crate::lidar::LidarConfig;
    use crate::scene::{Plane, Sphere};

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    fn single_beam_lidar() -> Lidar {
        Lidar::new(
            LidarConfig {
                azimuth_steps: 1,
                elevation_steps: 1,
                h_fov: 0.0,
                v_fov: 0.0,
                min_range: 0.0,
                max_range: 1_000.0,
                range_noise_std: 0.0,
            },
            0,
        )
        .unwrap()
    }

    #[test]
    fn coasting_holds_state_and_advances_clock() {
        let mut h = Harness::new(VehicleState::default(), Scene::new());
        let frame = h.step(&Command::coast(), 0.1).unwrap();
        assert!((frame.time - 0.1).abs() < 1e-12);
        assert!(
            frame.state.position.norm() < 1e-12,
            "at rest should not move"
        );
    }

    #[test]
    fn constant_acceleration_matches_kinematics() {
        // Body +x accel of 2 m/s², level: after the harness integrates with
        // semi-implicit Euler the position is deterministic and forward.
        let mut h = Harness::new(VehicleState::default(), Scene::new());
        let cmd = Command {
            accel_body: v(2.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        let dt = 0.01;
        let frames = h.run(&cmd, dt, 100).unwrap(); // 1 s total
        let last = frames.last().unwrap();
        // Velocity after 1 s of 2 m/s² ≈ 2 m/s (exact for constant accel).
        assert!(
            (last.state.velocity.x - 2.0).abs() < 1e-9,
            "v = {}",
            last.state.velocity.x
        );
        // Semi-implicit Euler position ≈ ½at² + ½a·dt·t; ~1.01 m here.
        assert!(
            (0.99..1.02).contains(&last.state.position.x),
            "x = {}",
            last.state.position.x
        );
        assert!((last.time - 1.0).abs() < 1e-9);
    }

    #[test]
    fn frame_contains_only_attached_sensors() {
        let mut h = Harness::new(VehicleState::default(), Scene::new())
            .with_imu(Imu::new(ImuConfig::default(), 0).unwrap());
        let frame = h.step(&Command::coast(), 0.1).unwrap();
        assert!(frame.imu.is_some());
        assert!(frame.lidar.is_none());
        assert!(frame.gps.is_none());
    }

    #[test]
    fn lidar_in_harness_ranges_the_scene() {
        // Wall at x = 5 (facing −x); a forward single-beam LiDAR on a level,
        // origin-placed vehicle reads 5 m.
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(5.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());
        let mut h = Harness::new(VehicleState::default(), scene).with_lidar(single_beam_lidar());
        let frame = h.sample();
        let scan = frame.lidar.unwrap();
        assert!((scan.beams[0].range.unwrap() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn imu_reads_one_g_at_rest_then_sees_commanded_accel() {
        let mut h = Harness::new(VehicleState::default(), Scene::new())
            .with_imu(Imu::new(ImuConfig::default(), 0).unwrap());
        // t = 0 sample: at rest ⇒ +g up.
        let f0 = h.sample();
        assert!((f0.imu.unwrap().specific_force.z - 9.806_65).abs() < 1e-9);
        // Step with +x body accel ⇒ specific force x-component ≈ +2.
        let cmd = Command {
            accel_body: v(2.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        let f1 = h.step(&cmd, 0.01).unwrap();
        let sf = f1.imu.unwrap().specific_force;
        assert!((sf.x - 2.0).abs() < 1e-9, "ax = {}", sf.x);
        assert!((sf.z - 9.806_65).abs() < 1e-9);
    }

    #[test]
    fn gps_in_harness_tracks_position() {
        let datum = Geodetic::from_degrees(40.0, -75.0, 100.0);
        let gps = Gps::new(datum, 0.0, 0.0, 0).unwrap();
        let mut h = Harness::new(VehicleState::default(), Scene::new()).with_gps(gps);
        // Move the vehicle 100 m east by command, then read the GPS.
        let cmd = Command {
            accel_body: v(1.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        // Place it directly to keep the test simple/deterministic.
        h.state.position = v(0.0, 0.0, 0.0);
        let f0 = h.sample();
        let g0 = f0.gps.unwrap();
        assert!((g0.position.lat - datum.lat).abs() < 1e-12, "at datum");

        h.state.position = v(0.0, 1_000.0, 0.0); // 1 km north
        let _ = cmd; // command unused for the manual-placement check
        let f1 = h.sample();
        assert!(
            f1.gps.unwrap().position.lat > datum.lat,
            "north ⇒ higher lat"
        );
    }

    #[test]
    fn turning_changes_heading_and_lidar_target() {
        // Wall only on +y; drive a yaw rate so that after enough time the
        // forward LiDAR beam swings onto it.
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(0.0, 6.0, 0.0), v(0.0, -1.0, 0.0)).unwrap());
        let mut h = Harness::new(VehicleState::default(), scene).with_lidar(single_beam_lidar());

        // Initially forward (+x) misses the +y wall.
        assert!(h.sample().lidar.unwrap().beams[0].range.is_none());

        // Yaw +90° over 1 s (π/2 rad/s), no translation.
        let cmd = Command {
            accel_body: v(0.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, std::f64::consts::FRAC_PI_2),
        };
        let frames = h.run(&cmd, 0.001, 1000).unwrap();
        // By the end the beam points ~+y and ranges the wall (~6 m, vehicle has
        // not translated).
        let last = frames.last().unwrap();
        let r = last.lidar.as_ref().unwrap().beams[0].range;
        assert!(r.is_some(), "beam should have swung onto the +y wall");
        assert!((r.unwrap() - 6.0).abs() < 0.1, "range = {r:?}");
    }

    #[test]
    fn full_stack_frame_has_all_sensors_and_is_deterministic() {
        let mut scene = Scene::new();
        scene.push_sphere(Sphere::new(v(8.0, 0.0, 0.0), 1.0).unwrap());
        let build = || {
            Harness::new(VehicleState::default(), scene.clone())
                .with_lidar(single_beam_lidar())
                .with_imu(Imu::new(ImuConfig::default(), 1).unwrap())
                .with_gps(Gps::new(Geodetic::from_degrees(0.0, 0.0, 0.0), 1.0, 1.0, 2).unwrap())
                .with_camera(
                    Camera::new(
                        Intrinsics::new(500.0, 500.0, 320.0, 240.0),
                        Distortion::none(),
                        640,
                        480,
                    )
                    .unwrap(),
                )
        };
        let cmd = Command {
            accel_body: v(0.5, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.01),
        };
        let mut a = build();
        let mut b = build();
        let fa = a.run(&cmd, 0.02, 50).unwrap();
        let fb = b.run(&cmd, 0.02, 50).unwrap();
        assert_eq!(fa, fb, "identical seeds ⇒ identical frame stream");
        let last = fa.last().unwrap();
        assert!(last.lidar.is_some() && last.imu.is_some() && last.gps.is_some());
        // Camera is attached but not sampled into the frame.
        assert!(a.camera().is_some());
    }

    #[test]
    fn bad_step_inputs_fail_loud() {
        let mut h = Harness::new(VehicleState::default(), Scene::new());
        assert!(h.step(&Command::coast(), 0.0).is_err());
        assert!(h.step(&Command::coast(), -0.1).is_err());
        assert!(h.step(&Command::coast(), f64::NAN).is_err());
        let bad = Command {
            accel_body: v(f64::NAN, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        assert!(h.step(&bad, 0.01).is_err());
    }
}
