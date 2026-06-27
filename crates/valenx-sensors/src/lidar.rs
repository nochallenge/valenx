//! Simulated **LiDAR** (light detection and ranging) range sensor.
//!
//! A [`Lidar`] emits a regular angular pattern of beams from a pose and
//! ray-casts each beam against an analytic [`crate::scene::Scene`], returning the
//! range (metres) to the nearest surface along that beam — or **no return** if
//! the beam misses everything or the nearest hit lies outside the
//! `[min_range, max_range]` window. This is the geometric heart of a scanning or
//! flash LiDAR; it omits the radiometry (no beam divergence, no intensity /
//! reflectivity, no multi-echo, no atmospheric attenuation).
//!
//! ## Beam pattern
//!
//! Beams are laid out on a grid of `azimuth_steps × elevation_steps` directions
//! spanning `[−h_fov/2, +h_fov/2]` in azimuth and `[−v_fov/2, +v_fov/2]` in
//! elevation, in the sensor body frame. The forward axis is **+x**, left is
//! **+y**, up is **+z** (a right-handed, REP-103-style body frame). A beam at
//! azimuth `α` and elevation `ε` points along
//! `(cosε·cosα, cosε·sinα, sinε)`. With a single beam each (`1 × 1`) the lone
//! beam points straight forward along +x.
//!
//! ## Noise
//!
//! Optional zero-mean Gaussian range noise (`range_noise_std`, metres) is added
//! to each *returning* beam from the crate's seeded [`crate::SplitMix64`]. With
//! `range_noise_std == 0` the ranges are the exact geometric distances, which is
//! what the tests pin.

use nalgebra::{UnitQuaternion, Vector3};

use crate::error::SensorError;
use crate::rng::SplitMix64;
use crate::scene::{Ray, Scene};

/// Configuration for a [`Lidar`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LidarConfig {
    /// Number of azimuth (horizontal) beam steps (≥ 1).
    pub azimuth_steps: usize,
    /// Number of elevation (vertical) beam steps (≥ 1).
    pub elevation_steps: usize,
    /// Horizontal field of view (radians, ≥ 0). `0` collapses the azimuth fan to
    /// the forward direction.
    pub h_fov: f64,
    /// Vertical field of view (radians, ≥ 0). `0` collapses the elevation fan to
    /// the horizon.
    pub v_fov: f64,
    /// Minimum reportable range (m, ≥ 0). Hits nearer than this are a no-return.
    pub min_range: f64,
    /// Maximum reportable range (m, > `min_range`). Hits beyond this are a
    /// no-return.
    pub max_range: f64,
    /// Standard deviation of additive zero-mean Gaussian range noise (m, ≥ 0).
    pub range_noise_std: f64,
}

impl Default for LidarConfig {
    /// A modest forward-looking scanner: a 32 × 8 beam grid over a 90° × 30°
    /// field of view, 0.2–100 m, noise-free.
    fn default() -> Self {
        Self {
            azimuth_steps: 32,
            elevation_steps: 8,
            h_fov: std::f64::consts::FRAC_PI_2, // 90°
            v_fov: std::f64::consts::FRAC_PI_6, // 30°
            min_range: 0.2,
            max_range: 100.0,
            range_noise_std: 0.0,
        }
    }
}

/// A simulated scanning LiDAR.
#[derive(Debug, Clone)]
pub struct Lidar {
    config: LidarConfig,
    rng: SplitMix64,
}

/// A single beam's measurement: its body-frame unit direction and the measured
/// range, where `range == None` is a **no-return** (miss or out of window).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beam {
    /// Azimuth angle of this beam (rad, in `[−h_fov/2, h_fov/2]`).
    pub azimuth: f64,
    /// Elevation angle of this beam (rad, in `[−v_fov/2, v_fov/2]`).
    pub elevation: f64,
    /// Body-frame unit direction of the beam.
    pub direction: Vector3<f64>,
    /// Measured range (m), or `None` for a no-return.
    pub range: Option<f64>,
}

/// A full LiDAR scan: one [`Beam`] per direction, row-major over
/// `elevation × azimuth`.
#[derive(Debug, Clone, PartialEq)]
pub struct LidarScan {
    /// The beams, in elevation-major then azimuth order.
    pub beams: Vec<Beam>,
}

impl LidarScan {
    /// Number of beams that returned a range (not a no-return).
    #[must_use]
    pub fn num_returns(&self) -> usize {
        self.beams.iter().filter(|b| b.range.is_some()).count()
    }

    /// All returning ranges, dropping no-returns.
    #[must_use]
    pub fn ranges(&self) -> Vec<f64> {
        self.beams.iter().filter_map(|b| b.range).collect()
    }
}

impl Lidar {
    /// Build a LiDAR from a config and a noise seed, validating the config.
    ///
    /// # Errors
    /// - [`SensorError::InvalidConfig`] if either step count is `0`, a field of
    ///   view is negative or non-finite, `min_range < 0`, or
    ///   `max_range <= min_range`.
    /// - [`SensorError::InvalidNoise`] if `range_noise_std` is negative or
    ///   non-finite.
    pub fn new(config: LidarConfig, seed: u64) -> Result<Self, SensorError> {
        if config.azimuth_steps == 0 || config.elevation_steps == 0 {
            return Err(SensorError::InvalidConfig(
                "azimuth_steps and elevation_steps must be ≥ 1".into(),
            ));
        }
        if !(config.h_fov.is_finite() && config.h_fov >= 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "h_fov must be finite and ≥ 0, got {}",
                config.h_fov
            )));
        }
        if !(config.v_fov.is_finite() && config.v_fov >= 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "v_fov must be finite and ≥ 0, got {}",
                config.v_fov
            )));
        }
        if !(config.min_range.is_finite() && config.min_range >= 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "min_range must be finite and ≥ 0, got {}",
                config.min_range
            )));
        }
        if !(config.max_range.is_finite() && config.max_range > config.min_range) {
            return Err(SensorError::InvalidConfig(format!(
                "max_range must be finite and > min_range ({} ≤ {})",
                config.max_range, config.min_range
            )));
        }
        if !(config.range_noise_std.is_finite() && config.range_noise_std >= 0.0) {
            return Err(SensorError::InvalidNoise(format!(
                "range_noise_std must be finite and ≥ 0, got {}",
                config.range_noise_std
            )));
        }
        Ok(Self {
            config,
            rng: SplitMix64::new(seed),
        })
    }

    /// The configuration this LiDAR was built with.
    #[must_use]
    pub fn config(&self) -> &LidarConfig {
        &self.config
    }

    /// The total number of beams (`azimuth_steps × elevation_steps`).
    #[must_use]
    pub fn beam_count(&self) -> usize {
        self.config.azimuth_steps * self.config.elevation_steps
    }

    /// The body-frame unit direction of a beam at azimuth `az` and elevation
    /// `el` (rad). Forward is +x, left +y, up +z.
    #[must_use]
    pub fn beam_direction(az: f64, el: f64) -> Vector3<f64> {
        let (sa, ca) = az.sin_cos();
        let (se, ce) = el.sin_cos();
        Vector3::new(ce * ca, ce * sa, se)
    }

    /// Evenly spaced angles for `steps` beams spanning `[−fov/2, fov/2]`.
    ///
    /// One step is the single centre angle `0`; `n > 1` steps are spaced by
    /// `fov / (n − 1)` so the first and last beams sit exactly on the fan edges.
    fn fan_angles(steps: usize, fov: f64) -> Vec<f64> {
        if steps <= 1 {
            return vec![0.0];
        }
        let start = -0.5 * fov;
        let step = fov / (steps - 1) as f64;
        (0..steps).map(|i| start + step * i as f64).collect()
    }

    /// Scan the `scene` from a sensor placed at `origin` with body orientation
    /// `orientation` (body → world). Returns one [`Beam`] per direction.
    ///
    /// Each beam direction is rotated from the body frame into the world frame,
    /// cast against the scene, and the nearest hit within
    /// `[min_range, max_range]` becomes the beam's range (plus optional noise);
    /// otherwise the beam is a no-return.
    pub fn scan(
        &mut self,
        scene: &Scene,
        origin: Vector3<f64>,
        orientation: UnitQuaternion<f64>,
    ) -> LidarScan {
        let azimuths = Self::fan_angles(self.config.azimuth_steps, self.config.h_fov);
        let elevations = Self::fan_angles(self.config.elevation_steps, self.config.v_fov);
        let mut beams = Vec::with_capacity(self.beam_count());

        for &el in &elevations {
            for &az in &azimuths {
                let body_dir = Self::beam_direction(az, el);
                let world_dir = orientation * body_dir;
                // body_dir is unit and a rotation preserves length, so Ray::new
                // cannot fail here; fall back to a no-return if it ever did.
                let range = match Ray::new(origin, world_dir) {
                    Ok(ray) => self.measure_along(scene, &ray),
                    Err(_) => None,
                };
                beams.push(Beam {
                    azimuth: az,
                    elevation: el,
                    direction: body_dir,
                    range,
                });
            }
        }
        LidarScan { beams }
    }

    /// Cast one prepared world-frame ray and apply the range window + noise.
    fn measure_along(&mut self, scene: &Scene, ray: &Ray) -> Option<f64> {
        let t = scene.cast(ray)?;
        if t < self.config.min_range || t > self.config.max_range {
            return None;
        }
        let measured = if self.config.range_noise_std > 0.0 {
            t + self.rng.next_normal(0.0, self.config.range_noise_std)
        } else {
            t
        };
        Some(measured.max(0.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::{Plane, Sphere};

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    fn single_beam_config() -> LidarConfig {
        LidarConfig {
            azimuth_steps: 1,
            elevation_steps: 1,
            h_fov: 0.0,
            v_fov: 0.0,
            min_range: 0.0,
            max_range: 1_000.0,
            range_noise_std: 0.0,
        }
    }

    #[test]
    fn single_forward_beam_ranges_a_known_plane_exactly() {
        // Wall at x = 7 m, sensor at origin facing +x: the lone beam must read
        // exactly 7 m.
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(7.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());
        let mut lidar = Lidar::new(single_beam_config(), 0).unwrap();
        let scan = lidar.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        assert_eq!(scan.beams.len(), 1);
        let r = scan.beams[0].range.expect("beam should return");
        assert!((r - 7.0).abs() < 1e-12, "range = {r}");
    }

    #[test]
    fn miss_is_a_no_return() {
        // Empty scene ⇒ the beam returns None.
        let scene = Scene::new();
        let mut lidar = Lidar::new(single_beam_config(), 0).unwrap();
        let scan = lidar.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        assert!(scan.beams[0].range.is_none());
        assert_eq!(scan.num_returns(), 0);
    }

    #[test]
    fn out_of_window_hits_are_no_returns() {
        // Plane at x = 7; window [0, 5] excludes it (too far); window starting at
        // 8 also excludes it (too near).
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(7.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());

        let mut cfg = single_beam_config();
        cfg.max_range = 5.0;
        let mut lidar = Lidar::new(cfg, 0).unwrap();
        let scan = lidar.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        assert!(
            scan.beams[0].range.is_none(),
            "7 m should exceed 5 m window"
        );

        cfg.max_range = 1_000.0;
        cfg.min_range = 8.0;
        let mut lidar = Lidar::new(cfg, 0).unwrap();
        let scan = lidar.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        assert!(scan.beams[0].range.is_none(), "7 m should be under 8 m min");
    }

    #[test]
    fn rotating_the_sensor_aims_the_beam() {
        // Plane at y = 4 (facing −y). Forward (+x) misses it; yaw +90° aims +y
        // and the beam reads 4 m.
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(0.0, 4.0, 0.0), v(0.0, -1.0, 0.0)).unwrap());
        let mut lidar = Lidar::new(single_beam_config(), 0).unwrap();

        let fwd = lidar.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        assert!(
            fwd.beams[0].range.is_none(),
            "forward should miss the +y wall"
        );

        let yaw = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), std::f64::consts::FRAC_PI_2);
        let turned = lidar.scan(&scene, v(0.0, 0.0, 0.0), yaw);
        let r = turned.beams[0].range.expect("yawed beam should hit");
        assert!((r - 4.0).abs() < 1e-9, "range = {r}");
    }

    #[test]
    fn full_grid_has_expected_beam_count_and_edges() {
        let cfg = LidarConfig {
            azimuth_steps: 5,
            elevation_steps: 3,
            h_fov: std::f64::consts::FRAC_PI_2,
            v_fov: std::f64::consts::FRAC_PI_6,
            min_range: 0.0,
            max_range: 100.0,
            range_noise_std: 0.0,
        };
        let lidar = Lidar::new(cfg, 0).unwrap();
        assert_eq!(lidar.beam_count(), 15);
        let az = Lidar::fan_angles(5, std::f64::consts::FRAC_PI_2);
        assert!((az[0] + std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        assert!((az[4] - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        assert!(az[2].abs() < 1e-12, "centre beam should be straight ahead");
    }

    #[test]
    fn sphere_in_grid_returns_a_cluster_of_hits() {
        // A sphere dead ahead should be hit by the central beams of a forward
        // grid and missed by the extreme ones.
        let mut scene = Scene::new();
        scene.push_sphere(Sphere::new(v(10.0, 0.0, 0.0), 1.0).unwrap());
        let cfg = LidarConfig {
            azimuth_steps: 21,
            elevation_steps: 1,
            h_fov: std::f64::consts::FRAC_PI_2,
            v_fov: 0.0,
            min_range: 0.0,
            max_range: 100.0,
            range_noise_std: 0.0,
        };
        let mut lidar = Lidar::new(cfg, 0).unwrap();
        let scan = lidar.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        // Some beams hit, but not all (the sphere subtends only a small angle).
        let hits = scan.num_returns();
        assert!(
            hits > 0 && hits < 21,
            "expected a partial cluster, got {hits}"
        );
        // The central beam hits the near face at t = 9.
        let centre = &scan.beams[10];
        assert!(centre.azimuth.abs() < 1e-12);
        assert!((centre.range.unwrap() - 9.0).abs() < 1e-9);
    }

    #[test]
    fn noise_is_deterministic_for_a_seed() {
        let mut scene = Scene::new();
        scene.push_plane(Plane::new(v(7.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());
        let mut cfg = single_beam_config();
        cfg.range_noise_std = 0.05;

        let mut a = Lidar::new(cfg, 42).unwrap();
        let mut b = Lidar::new(cfg, 42).unwrap();
        let ra = a.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        let rb = b.scan(&scene, v(0.0, 0.0, 0.0), UnitQuaternion::identity());
        assert_eq!(ra, rb, "same seed ⇒ identical noisy scan");
        // Noisy range is near, but not exactly, the truth.
        let r = ra.beams[0].range.unwrap();
        assert!((r - 7.0).abs() < 0.5 && (r - 7.0).abs() > 0.0);
    }

    #[test]
    fn invalid_configs_are_rejected() {
        let mut cfg = single_beam_config();
        cfg.azimuth_steps = 0;
        assert!(Lidar::new(cfg, 0).is_err());

        let mut cfg = single_beam_config();
        cfg.max_range = cfg.min_range; // not strictly greater
        assert!(Lidar::new(cfg, 0).is_err());

        let mut cfg = single_beam_config();
        cfg.range_noise_std = -1.0;
        assert!(Lidar::new(cfg, 0).is_err());

        let mut cfg = single_beam_config();
        cfg.h_fov = -0.1;
        assert!(Lidar::new(cfg, 0).is_err());
    }
}
