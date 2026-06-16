//! Closed-form vacuum (drag-free) projectile motion.
//!
//! For a point mass launched from flat level ground at speed `v0` and
//! elevation angle `θ` in a uniform gravity field `g`, classical
//! kinematics give exact closed forms (no integration needed):
//!
//! - Time of flight: `T = 2·v0·sin(θ) / g`
//! - Horizontal range: `R = v0²·sin(2θ) / g`
//! - Apex (maximum) height: `H = v0²·sin²(θ) / (2g)`
//!
//! The range is maximised at `θ = 45°`, where `sin(2θ) = 1`, giving the
//! well-known `R_max = v0² / g`. These identities are the ground-truth
//! checks exercised by the unit tests below. The range relation inverts
//! to the two complementary launch angles that hit a target range
//! ([`vacuum_angles_for_range`]).

use serde::{Deserialize, Serialize};

use crate::error::{ProjectileError, Result};

/// Standard gravity at the Earth's surface, in metres per second squared
/// (the CODATA / ISO 80000 conventional value `g₀ = 9.80665 m/s²`).
pub const STANDARD_GRAVITY: f64 = 9.806_65;

/// A drag-free launch condition on flat, level ground.
///
/// All fields use SI units: speed in m/s, angle in radians, gravity in
/// m/s². Construct via [`VacuumShot::new`], which validates the inputs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VacuumShot {
    /// Launch speed `v0` (m/s), strictly positive.
    pub speed: f64,
    /// Launch elevation angle `θ` (radians) in `[0, π/2]`.
    pub angle_rad: f64,
    /// Gravitational acceleration `g` (m/s²), strictly positive.
    pub gravity: f64,
}

impl VacuumShot {
    /// Build a validated vacuum shot.
    ///
    /// # Errors
    ///
    /// Returns [`ProjectileError::InvalidParameter`] if `speed` or
    /// `gravity` is not finite-and-positive, or
    /// [`ProjectileError::AngleOutOfRange`] if `angle_rad` is outside
    /// `[0, π/2]`.
    pub fn new(speed: f64, angle_rad: f64, gravity: f64) -> Result<Self> {
        let speed = ProjectileError::require_positive("speed", speed)?;
        let gravity = ProjectileError::require_positive("gravity", gravity)?;
        let angle_rad = ProjectileError::require_angle(angle_rad)?;
        Ok(Self {
            speed,
            angle_rad,
            gravity,
        })
    }

    /// Build a vacuum shot using [`STANDARD_GRAVITY`] for `g`.
    ///
    /// # Errors
    ///
    /// As [`VacuumShot::new`], minus the gravity check (the constant is
    /// always valid).
    pub fn at_standard_gravity(speed: f64, angle_rad: f64) -> Result<Self> {
        Self::new(speed, angle_rad, STANDARD_GRAVITY)
    }

    /// Time of flight back to the launch height: `T = 2·v0·sin(θ) / g`.
    #[must_use]
    pub fn time_of_flight(&self) -> f64 {
        2.0 * self.speed * self.angle_rad.sin() / self.gravity
    }

    /// Horizontal range on level ground: `R = v0²·sin(2θ) / g`.
    #[must_use]
    pub fn range(&self) -> f64 {
        self.speed * self.speed * (2.0 * self.angle_rad).sin() / self.gravity
    }

    /// Apex (maximum) height: `H = v0²·sin²(θ) / (2g)`.
    #[must_use]
    pub fn apex_height(&self) -> f64 {
        let s = self.angle_rad.sin();
        self.speed * self.speed * s * s / (2.0 * self.gravity)
    }

    /// Time to reach the apex: `t_apex = v0·sin(θ) / g`.
    ///
    /// By symmetry this is exactly half [`VacuumShot::time_of_flight`].
    #[must_use]
    pub fn time_to_apex(&self) -> f64 {
        self.speed * self.angle_rad.sin() / self.gravity
    }
}

/// The launch angle (radians) that maximises vacuum range on level
/// ground: always exactly `π/4` (45°), independent of speed and gravity,
/// because `R ∝ sin(2θ)` peaks at `2θ = π/2`.
#[must_use]
pub fn optimal_vacuum_angle_rad() -> f64 {
    core::f64::consts::FRAC_PI_4
}

/// Maximum achievable vacuum range on level ground: `R_max = v0² / g`,
/// the value of [`VacuumShot::range`] at the 45° optimum.
///
/// # Errors
///
/// Returns [`ProjectileError::InvalidParameter`] if `speed` or `gravity`
/// is not finite-and-positive.
pub fn max_vacuum_range(speed: f64, gravity: f64) -> Result<f64> {
    let speed = ProjectileError::require_positive("speed", speed)?;
    let gravity = ProjectileError::require_positive("gravity", gravity)?;
    Ok(speed * speed / gravity)
}

/// The two launch angles (radians) that reach a target horizontal range
/// `range_target` in vacuum on level ground, inverting
/// [`VacuumShot::range`]:
///
/// ```text
/// sin(2θ) = R·g / v0²  ->  θ_low  = arcsin(R·g/v0²) / 2,
///                          θ_high = π/2 − θ_low
/// ```
///
/// Returned as `(low, high)` with `low <= π/4 <= high`: the **direct**
/// (flatter) trajectory and the **lofted** (higher) one. They are
/// complementary (`low + high = π/2`, the classic `sin(2θ)` symmetry)
/// and coincide at `π/4` when the target equals the maximum range
/// `v0²/g`.
///
/// # Errors
///
/// Returns [`ProjectileError::InvalidParameter`] if `speed` or `gravity`
/// is not finite-and-positive, if `range_target` is negative or
/// non-finite, or if `range_target` exceeds the maximum vacuum range
/// `v0²/g` (no launch angle can reach it).
pub fn vacuum_angles_for_range(speed: f64, range_target: f64, gravity: f64) -> Result<(f64, f64)> {
    let speed = ProjectileError::require_positive("speed", speed)?;
    let gravity = ProjectileError::require_positive("gravity", gravity)?;
    let range_target = ProjectileError::require_non_negative("range_target", range_target)?;
    let s = range_target * gravity / (speed * speed);
    // Reject genuinely unreachable targets, but tolerate a hair over 1.0
    // from rounding at exactly the maximum range.
    if s > 1.0 + 1e-9 {
        return Err(ProjectileError::InvalidParameter {
            field: "range_target",
            value: range_target,
            reason: "exceeds the maximum vacuum range v0^2/g",
        });
    }
    // Clamp the rounding overshoot so asin stays inside its [-1, 1] domain.
    let low = s.min(1.0).asin() / 2.0;
    let high = core::f64::consts::FRAC_PI_2 - low;
    Ok((low, high))
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    const EPS: f64 = 1e-9;

    /// The 45° range equals the closed-form maximum `v0²/g`, and the
    /// dedicated `optimal_vacuum_angle_rad` helper returns π/4.
    #[test]
    fn optimal_angle_is_45_degrees() {
        let v0 = 30.0;
        let g = STANDARD_GRAVITY;
        let opt = optimal_vacuum_angle_rad();
        assert!((opt - FRAC_PI_4).abs() < EPS);

        let shot = VacuumShot::new(v0, opt, g).unwrap();
        let expected_max = max_vacuum_range(v0, g).unwrap();
        assert!((shot.range() - expected_max).abs() < 1e-6);
        assert!((expected_max - v0 * v0 / g).abs() < EPS);
    }

    /// Vacuum range is genuinely maximised at 45°: sample the whole
    /// quadrant and confirm no angle beats the 45° range.
    #[test]
    fn range_is_maximised_at_45_degrees() {
        let v0 = 42.0;
        let g = STANDARD_GRAVITY;
        let r45 = VacuumShot::new(v0, FRAC_PI_4, g).unwrap().range();

        // 1-degree sweep across [0, 90]; every sample is <= the 45° range.
        for deg in 0..=90 {
            let theta = (deg as f64).to_radians();
            let r = VacuumShot::new(v0, theta, g).unwrap().range();
            assert!(
                r <= r45 + 1e-9,
                "angle {deg}° gave range {r} > 45° range {r45}"
            );
        }
    }

    /// The range matches the analytic `v0²·sin(2θ)/g` at several angles.
    #[test]
    fn range_matches_analytic_formula() {
        let v0 = 25.0;
        let g = 9.81;
        for deg in [10.0_f64, 20.0, 30.0, 45.0, 60.0, 75.0] {
            let theta = deg.to_radians();
            let shot = VacuumShot::new(v0, theta, g).unwrap();
            let analytic = v0 * v0 * (2.0 * theta).sin() / g;
            assert!(
                (shot.range() - analytic).abs() < EPS,
                "range mismatch at {deg}°"
            );
        }
    }

    /// Complementary angles (θ and 90°−θ) give equal range — a classic
    /// property of `sin(2θ)`.
    #[test]
    fn complementary_angles_share_range() {
        let v0 = 18.0;
        let g = STANDARD_GRAVITY;
        let a = 30.0_f64.to_radians();
        let b = FRAC_PI_2 - a; // 60 degrees
        let ra = VacuumShot::new(v0, a, g).unwrap().range();
        let rb = VacuumShot::new(v0, b, g).unwrap().range();
        assert!((ra - rb).abs() < 1e-9);
    }

    /// Apex height equals the analytic `v0²·sin²(θ)/(2g)`.
    #[test]
    fn apex_matches_analytic_formula() {
        let v0 = 50.0;
        let g = 9.81;
        for deg in [15.0_f64, 45.0, 90.0] {
            let theta = deg.to_radians();
            let shot = VacuumShot::new(v0, theta, g).unwrap();
            let s = theta.sin();
            let analytic = v0 * v0 * s * s / (2.0 * g);
            assert!(
                (shot.apex_height() - analytic).abs() < 1e-9,
                "apex mismatch at {deg}°"
            );
        }
    }

    /// A vertical (90°) shot: range is zero, apex is `v0²/(2g)`, and the
    /// total flight time is `2·v0/g`.
    #[test]
    fn vertical_shot_known_values() {
        let v0 = 20.0;
        let g = 10.0;
        let shot = VacuumShot::new(v0, FRAC_PI_2, g).unwrap();
        assert!(shot.range().abs() < 1e-9);
        assert!((shot.apex_height() - v0 * v0 / (2.0 * g)).abs() < EPS);
        assert!((shot.time_of_flight() - 2.0 * v0 / g).abs() < EPS);
    }

    /// Time to apex is exactly half the time of flight (launch/return
    /// symmetry).
    #[test]
    fn time_to_apex_is_half_flight() {
        let shot = VacuumShot::new(33.0, 37.0_f64.to_radians(), 9.81).unwrap();
        assert!((2.0 * shot.time_to_apex() - shot.time_of_flight()).abs() < EPS);
    }

    /// The range-vs-flight-time identity `R = v0·cos(θ)·T` must hold,
    /// since horizontal velocity is constant in vacuum.
    #[test]
    fn range_equals_horizontal_velocity_times_time() {
        let v0 = 27.5;
        let g = 9.80665;
        let theta = 55.0_f64.to_radians();
        let shot = VacuumShot::new(v0, theta, g).unwrap();
        let by_kinematics = v0 * theta.cos() * shot.time_of_flight();
        assert!((shot.range() - by_kinematics).abs() < 1e-9);
    }

    /// `max_vacuum_range` reproduces the full quadrant peak to within the
    /// sampling resolution and equals `v0²/g` exactly.
    #[test]
    fn max_range_helper_matches_sampled_peak() {
        let v0 = 12.3;
        let g = 9.81;
        let analytic = max_vacuum_range(v0, g).unwrap();
        let mut sampled_peak = 0.0_f64;
        // Fine sweep to approach the continuous maximum.
        let n = 9_000;
        for i in 0..=n {
            let theta = PI / 2.0 * (i as f64 / n as f64);
            let r = VacuumShot::new(v0, theta, g).unwrap().range();
            sampled_peak = sampled_peak.max(r);
        }
        assert!((analytic - v0 * v0 / g).abs() < EPS);
        assert!((sampled_peak - analytic).abs() < 1e-3);
    }

    /// Both inverse-solved angles reproduce the target range exactly.
    #[test]
    fn angles_for_range_round_trip_through_range() {
        let v0 = 40.0;
        let g = STANDARD_GRAVITY;
        let target = 0.6 * max_vacuum_range(v0, g).unwrap();
        let (low, high) = vacuum_angles_for_range(v0, target, g).unwrap();
        let r_low = VacuumShot::new(v0, low, g).unwrap().range();
        let r_high = VacuumShot::new(v0, high, g).unwrap().range();
        assert!((r_low - target).abs() < 1e-9, "low angle range {r_low}");
        assert!((r_high - target).abs() < 1e-9, "high angle range {r_high}");
    }

    /// The two solutions are complementary and straddle 45°.
    #[test]
    fn angles_are_complementary_and_straddle_45() {
        let (low, high) = vacuum_angles_for_range(33.0, 50.0, 9.81).unwrap();
        assert!((low + high - FRAC_PI_2).abs() < 1e-12, "low+high");
        assert!(low <= FRAC_PI_4 + 1e-12 && high >= FRAC_PI_4 - 1e-12);
    }

    /// At exactly the maximum range both angles collapse onto 45°.
    #[test]
    fn at_max_range_both_angles_are_45() {
        let v0 = 25.0;
        let g = STANDARD_GRAVITY;
        let r_max = max_vacuum_range(v0, g).unwrap();
        let (low, high) = vacuum_angles_for_range(v0, r_max, g).unwrap();
        assert!((low - FRAC_PI_4).abs() < 1e-6, "low {low}");
        assert!((high - FRAC_PI_4).abs() < 1e-6, "high {high}");
    }

    /// Zero range is reached by a flat (0°) or vertical (90°) shot.
    #[test]
    fn zero_range_gives_horizontal_and_vertical() {
        let (low, high) = vacuum_angles_for_range(20.0, 0.0, 9.81).unwrap();
        assert!(low.abs() < 1e-12, "low {low}");
        assert!((high - FRAC_PI_2).abs() < 1e-12, "high {high}");
    }

    /// Known closed-form solution: for theta_low = 15°, sin(30°) = 0.5, so
    /// R = 0.5 v0^2/g and the complement is 75°.
    #[test]
    fn angles_match_known_low_solution() {
        let v0 = 30.0;
        let g = 10.0;
        let target = 0.5 * v0 * v0 / g; // 45 m
        let (low, high) = vacuum_angles_for_range(v0, target, g).unwrap();
        assert!((low - 15.0_f64.to_radians()).abs() < 1e-9, "low {low}");
        assert!((high - 75.0_f64.to_radians()).abs() < 1e-9, "high {high}");
    }

    /// Unreachable targets and non-physical inputs are rejected; exactly
    /// the maximum range still succeeds (rounding-robust).
    #[test]
    fn angles_for_range_rejects_unreachable_and_bad_inputs() {
        let v0 = 20.0;
        let g = 9.81;
        let r_max = max_vacuum_range(v0, g).unwrap();
        assert!(vacuum_angles_for_range(v0, r_max * 1.01, g).is_err());
        assert!(vacuum_angles_for_range(v0, -1.0, g).is_err());
        assert!(vacuum_angles_for_range(0.0, 10.0, g).is_err());
        assert!(vacuum_angles_for_range(v0, f64::NAN, g).is_err());
        assert!(vacuum_angles_for_range(v0, 10.0, 0.0).is_err());
        assert!(vacuum_angles_for_range(v0, r_max, g).is_ok());
    }

    /// Constructors reject non-physical inputs.
    #[test]
    fn rejects_invalid_inputs() {
        assert!(matches!(
            VacuumShot::new(-1.0, FRAC_PI_4, 9.81),
            Err(ProjectileError::InvalidParameter { field: "speed", .. })
        ));
        assert!(matches!(
            VacuumShot::new(10.0, FRAC_PI_4, 0.0),
            Err(ProjectileError::InvalidParameter {
                field: "gravity",
                ..
            })
        ));
        assert!(matches!(
            VacuumShot::new(10.0, -0.1, 9.81),
            Err(ProjectileError::AngleOutOfRange { .. })
        ));
        assert!(matches!(
            VacuumShot::new(10.0, FRAC_PI_2 + 0.1, 9.81),
            Err(ProjectileError::AngleOutOfRange { .. })
        ));
        assert!(matches!(
            max_vacuum_range(f64::NAN, 9.81),
            Err(ProjectileError::InvalidParameter { field: "speed", .. })
        ));
    }
}
