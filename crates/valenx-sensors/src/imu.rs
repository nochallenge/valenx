//! Simulated strap-down **IMU** (inertial measurement unit): a 3-axis
//! accelerometer + 3-axis gyroscope.
//!
//! An IMU does **not** measure acceleration directly — an accelerometer measures
//! **specific force**: the non-gravitational force per unit mass, in the *body*
//! frame. From a rigid-body trajectory with world-frame acceleration `a`, body
//! orientation `R` (body → world), and gravity `g` (a world-frame vector, e.g.
//! `(0, 0, −9.81)`), the ideal accelerometer reading is
//!
//! ```text
//! f_body = Rᵀ · (a − g)
//! ```
//!
//! so a body **at rest** (`a = 0`) on the ground reads `Rᵀ·(−g)`, i.e. `+g`
//! upward along its own up-axis ("the 1 g you feel sitting still"), not zero.
//! The gyroscope measures the body angular rate `ω`, already a body-frame
//! quantity.
//!
//! On top of the ideal kinematics this model adds a constant **bias** and
//! zero-mean **Gaussian noise** to each axis, drawn from the crate's seeded
//! [`crate::SplitMix64`] so a run is reproducible. With zero bias and zero noise
//! the reading is the exact kinematic specific force / angular rate, which is
//! what the tests pin.
//!
//! Honest scope: this is the standard *ideal + bias + white-noise* model. It
//! omits scale-factor and axis-misalignment errors, a random-walk (Brownian)
//! bias drift, temperature dependence, `g`-sensitivity, and quantisation — the
//! second-order effects a full IMU error model (and a Kalman filter's process
//! model) would include.

use nalgebra::{UnitQuaternion, Vector3};

use crate::error::SensorError;
use crate::rng::SplitMix64;

/// Per-axis bias + noise specification for one 3-axis sensor.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AxisNoise {
    /// Constant additive bias applied to every axis (same units as the sensor).
    pub bias: f64,
    /// Standard deviation of zero-mean Gaussian noise per axis (≥ 0).
    pub std: f64,
}

impl AxisNoise {
    /// No bias and no noise — the ideal sensor.
    #[must_use]
    pub fn none() -> Self {
        Self::default()
    }

    /// A bias and noise standard deviation.
    #[must_use]
    pub fn new(bias: f64, std: f64) -> Self {
        Self { bias, std }
    }

    fn validate(&self, what: &str) -> Result<(), SensorError> {
        if !self.bias.is_finite() {
            return Err(SensorError::InvalidNoise(format!(
                "{what} bias must be finite"
            )));
        }
        if !(self.std.is_finite() && self.std >= 0.0) {
            return Err(SensorError::InvalidNoise(format!(
                "{what} std must be finite and ≥ 0, got {}",
                self.std
            )));
        }
        Ok(())
    }
}

/// IMU configuration: gravity vector plus accel/gyro bias-and-noise specs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImuConfig {
    /// World-frame gravity vector (m/s²), e.g. `(0, 0, −9.80665)`.
    pub gravity: Vector3<f64>,
    /// Accelerometer bias + noise (m/s²).
    pub accel: AxisNoise,
    /// Gyroscope bias + noise (rad/s).
    pub gyro: AxisNoise,
}

impl Default for ImuConfig {
    /// Standard gravity pointing down `−z`, ideal (noise-free) sensors.
    fn default() -> Self {
        Self {
            gravity: Vector3::new(0.0, 0.0, -9.806_65),
            accel: AxisNoise::none(),
            gyro: AxisNoise::none(),
        }
    }
}

/// One rigid-body kinematic sample the IMU measures from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyState {
    /// Body orientation (body → world).
    pub orientation: UnitQuaternion<f64>,
    /// World-frame linear acceleration of the body (m/s²).
    pub accel_world: Vector3<f64>,
    /// Body-frame angular rate (rad/s).
    pub angular_rate_body: Vector3<f64>,
}

/// One IMU reading: body-frame specific force and angular rate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImuReading {
    /// Accelerometer: specific force in the body frame (m/s²).
    pub specific_force: Vector3<f64>,
    /// Gyroscope: angular rate in the body frame (rad/s).
    pub angular_rate: Vector3<f64>,
}

/// A simulated strap-down IMU.
#[derive(Debug, Clone)]
pub struct Imu {
    config: ImuConfig,
    rng: SplitMix64,
}

impl Imu {
    /// Build an IMU from a config and a noise seed, validating the noise specs.
    ///
    /// # Errors
    /// - [`SensorError::NonFinite`] if the gravity vector has a non-finite
    ///   component.
    /// - [`SensorError::InvalidNoise`] if an accel/gyro bias or standard
    ///   deviation is non-finite, or a standard deviation is negative.
    pub fn new(config: ImuConfig, seed: u64) -> Result<Self, SensorError> {
        if !config.gravity.iter().all(|c| c.is_finite()) {
            return Err(SensorError::NonFinite("gravity vector".into()));
        }
        config.accel.validate("accel")?;
        config.gyro.validate("gyro")?;
        Ok(Self {
            config,
            rng: SplitMix64::new(seed),
        })
    }

    /// The configuration this IMU was built with.
    #[must_use]
    pub fn config(&self) -> &ImuConfig {
        &self.config
    }

    /// The **ideal** (bias-free, noise-free) specific force for a body state:
    /// `Rᵀ·(a − g)`.
    #[must_use]
    pub fn ideal_specific_force(&self, state: &BodyState) -> Vector3<f64> {
        let world_specific = state.accel_world - self.config.gravity;
        state.orientation.inverse() * world_specific
    }

    /// Sample the IMU at a body state, adding bias and seeded Gaussian noise.
    ///
    /// The accelerometer returns `Rᵀ·(a − g)` + bias + noise; the gyroscope
    /// returns the body angular rate + bias + noise.
    pub fn sample(&mut self, state: &BodyState) -> ImuReading {
        let ideal_f = self.ideal_specific_force(state);
        let specific_force = self.add_noise(ideal_f, self.config.accel);
        let angular_rate = self.add_noise(state.angular_rate_body, self.config.gyro);
        ImuReading {
            specific_force,
            angular_rate,
        }
    }

    /// Add per-axis bias + Gaussian noise to a 3-vector.
    fn add_noise(&mut self, v: Vector3<f64>, spec: AxisNoise) -> Vector3<f64> {
        Vector3::new(
            v.x + spec.bias + self.rng.next_normal(0.0, spec.std),
            v.y + spec.bias + self.rng.next_normal(0.0, spec.std),
            v.z + spec.bias + self.rng.next_normal(0.0, spec.std),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    fn at_rest_level() -> BodyState {
        BodyState {
            orientation: UnitQuaternion::identity(),
            accel_world: v(0.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        }
    }

    #[test]
    fn body_at_rest_reads_plus_one_g_up() {
        // Level, at rest: specific force = Rᵀ·(0 − g) = −g = +9.80665 on +z.
        let mut imu = Imu::new(ImuConfig::default(), 0).unwrap();
        let r = imu.sample(&at_rest_level());
        assert!((r.specific_force.x).abs() < 1e-12);
        assert!((r.specific_force.y).abs() < 1e-12);
        assert!(
            (r.specific_force.z - 9.806_65).abs() < 1e-12,
            "az = {}",
            r.specific_force.z
        );
        assert!(r.angular_rate.norm() < 1e-12);
    }

    #[test]
    fn free_fall_reads_zero_specific_force() {
        // In free fall a = g, so specific force = Rᵀ·(g − g) = 0 (weightless).
        let mut imu = Imu::new(ImuConfig::default(), 0).unwrap();
        let state = BodyState {
            orientation: UnitQuaternion::identity(),
            accel_world: v(0.0, 0.0, -9.806_65),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        let r = imu.sample(&state);
        assert!(
            r.specific_force.norm() < 1e-12,
            "free fall ⇒ ~0, got {}",
            r.specific_force.norm()
        );
    }

    #[test]
    fn horizontal_acceleration_appears_on_the_body_x_axis() {
        // Level body accelerating +2 m/s² along world +x: specific force =
        // Rᵀ·(a − g) = (2, 0, +9.80665).
        let mut imu = Imu::new(ImuConfig::default(), 0).unwrap();
        let state = BodyState {
            orientation: UnitQuaternion::identity(),
            accel_world: v(2.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        let r = imu.sample(&state);
        assert!((r.specific_force.x - 2.0).abs() < 1e-12);
        assert!((r.specific_force.z - 9.806_65).abs() < 1e-12);
    }

    #[test]
    fn orientation_rotates_gravity_into_the_body_frame() {
        // Pitch the body 90° about +y so its body +z points along world +x.
        // Rᵀ·(−g) puts the 1 g on the body −x axis.
        let pitch =
            UnitQuaternion::from_axis_angle(&Vector3::y_axis(), std::f64::consts::FRAC_PI_2);
        let mut imu = Imu::new(ImuConfig::default(), 0).unwrap();
        let state = BodyState {
            orientation: pitch,
            accel_world: v(0.0, 0.0, 0.0),
            angular_rate_body: v(0.0, 0.0, 0.0),
        };
        let r = imu.sample(&state);
        // Magnitude is still 1 g; it now lies on body x.
        assert!((r.specific_force.norm() - 9.806_65).abs() < 1e-9);
        assert!(
            (r.specific_force.x + 9.806_65).abs() < 1e-9,
            "f = {:?}",
            r.specific_force
        );
        assert!(r.specific_force.z.abs() < 1e-9);
    }

    #[test]
    fn gyro_reports_the_body_angular_rate() {
        let mut imu = Imu::new(ImuConfig::default(), 0).unwrap();
        let state = BodyState {
            orientation: UnitQuaternion::identity(),
            accel_world: v(0.0, 0.0, 0.0),
            angular_rate_body: v(0.1, -0.2, 0.3),
        };
        let r = imu.sample(&state);
        assert!((r.angular_rate - v(0.1, -0.2, 0.3)).norm() < 1e-12);
    }

    #[test]
    fn bias_offsets_every_axis() {
        let cfg = ImuConfig {
            accel: AxisNoise::new(0.5, 0.0),
            gyro: AxisNoise::new(0.01, 0.0),
            ..ImuConfig::default()
        };
        let mut imu = Imu::new(cfg, 0).unwrap();
        let r = imu.sample(&at_rest_level());
        // Each accel axis is shifted by +0.5 from the ideal (0,0,+g).
        assert!((r.specific_force.x - 0.5).abs() < 1e-12);
        assert!((r.specific_force.z - (9.806_65 + 0.5)).abs() < 1e-12);
        assert!((r.angular_rate.x - 0.01).abs() < 1e-12);
    }

    #[test]
    fn noise_is_deterministic_and_small() {
        let cfg = ImuConfig {
            accel: AxisNoise::new(0.0, 0.02),
            gyro: AxisNoise::new(0.0, 0.001),
            ..ImuConfig::default()
        };
        let mut a = Imu::new(cfg, 99).unwrap();
        let mut b = Imu::new(cfg, 99).unwrap();
        let ra = a.sample(&at_rest_level());
        let rb = b.sample(&at_rest_level());
        assert_eq!(ra, rb, "same seed ⇒ identical noisy reading");
        // Reading is near the ideal +g on z, perturbed by noise.
        assert!((ra.specific_force.z - 9.806_65).abs() < 0.2);
    }

    #[test]
    fn invalid_noise_is_rejected() {
        // Negative accel noise std.
        let cfg = ImuConfig {
            accel: AxisNoise::new(0.0, -1.0),
            ..ImuConfig::default()
        };
        assert!(Imu::new(cfg, 0).is_err());

        // Non-finite gyro bias.
        let cfg = ImuConfig {
            gyro: AxisNoise::new(f64::NAN, 0.0),
            ..ImuConfig::default()
        };
        assert!(Imu::new(cfg, 0).is_err());

        // Non-finite gravity.
        let cfg = ImuConfig {
            gravity: v(0.0, 0.0, f64::INFINITY),
            ..ImuConfig::default()
        };
        assert!(Imu::new(cfg, 0).is_err());
    }
}
