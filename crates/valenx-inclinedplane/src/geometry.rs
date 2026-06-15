//! Ideal (frictionless) geometry and mechanical advantage of a ramp.
//!
//! This module captures the purely kinematic side of the inclined
//! plane: the slope angle and the velocity-ratio mechanical advantage
//! `MA = 1 / sin(theta)`. Friction and weight live in
//! [`crate::statics`].

use serde::{Deserialize, Serialize};

use crate::error::InclinedPlaneError;

/// Standard gravitational acceleration at the Earth's surface, in
/// metres per second squared (CODATA / ISO 80000 conventional value
/// `g_n`). Provided as a convenience default for callers that work in
/// SI units and want to turn a mass into a weight.
pub const GRAVITY_STANDARD: f64 = 9.806_65;

/// A frictionless inclined plane defined solely by its slope angle.
///
/// The angle is stored in radians and is constrained to the open
/// interval `(0, pi/2)` — a ramp at `0` is flat (infinite mechanical
/// advantage, no lifting) and a ramp at `pi/2` is vertical (mechanical
/// advantage `1`, i.e. a straight lift). Use [`IdealRamp::from_angle`]
/// or [`IdealRamp::from_rise_run`] to build one.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IdealRamp {
    /// Slope angle from the horizontal, in radians, in `(0, pi/2)`.
    angle_rad: f64,
}

impl IdealRamp {
    /// Build a ramp from its slope angle in radians.
    ///
    /// # Errors
    ///
    /// Returns [`InclinedPlaneError::BadParameter`] if `angle_rad` is
    /// not finite or does not lie in the open interval `(0, pi/2)`.
    pub fn from_angle(angle_rad: f64) -> Result<Self, InclinedPlaneError> {
        if !angle_rad.is_finite() {
            return Err(InclinedPlaneError::bad_parameter(
                "angle_rad",
                format!("must be finite, got {angle_rad}"),
            ));
        }
        if angle_rad <= 0.0 || angle_rad >= std::f64::consts::FRAC_PI_2 {
            return Err(InclinedPlaneError::bad_parameter(
                "angle_rad",
                format!("must lie in (0, pi/2) radians, got {angle_rad}"),
            ));
        }
        Ok(Self { angle_rad })
    }

    /// Build a ramp from a vertical rise and horizontal run.
    ///
    /// The slope angle is `atan(rise / run)`. Both lengths must be
    /// finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`InclinedPlaneError::BadParameter`] if either length is
    /// non-finite or non-positive.
    pub fn from_rise_run(rise: f64, run: f64) -> Result<Self, InclinedPlaneError> {
        if !rise.is_finite() || rise <= 0.0 {
            return Err(InclinedPlaneError::bad_parameter(
                "rise",
                format!("must be finite and > 0, got {rise}"),
            ));
        }
        if !run.is_finite() || run <= 0.0 {
            return Err(InclinedPlaneError::bad_parameter(
                "run",
                format!("must be finite and > 0, got {run}"),
            ));
        }
        Self::from_angle(rise.atan2(run))
    }

    /// The slope angle in radians.
    pub fn angle_rad(&self) -> f64 {
        self.angle_rad
    }

    /// The slope angle in degrees.
    pub fn angle_deg(&self) -> f64 {
        self.angle_rad.to_degrees()
    }

    /// Ideal (frictionless) mechanical advantage, the velocity ratio
    /// `MA = 1 / sin(theta)`.
    ///
    /// This is how many times the slope-parallel effort is multiplied
    /// into lifting force when friction is neglected. It always exceeds
    /// `1` for a non-vertical ramp and grows without bound as the ramp
    /// flattens.
    pub fn ideal_mechanical_advantage(&self) -> f64 {
        1.0 / self.angle_rad.sin()
    }

    /// Length of incline (hypotenuse) travelled to gain a given
    /// vertical `height`.
    ///
    /// Moving along the ramp by `height / sin(theta)` raises the load
    /// by `height`; this is the distance the effort acts over, and the
    /// inverse of the mechanical-advantage trade-off (more distance,
    /// less force).
    ///
    /// # Errors
    ///
    /// Returns [`InclinedPlaneError::BadParameter`] if `height` is
    /// non-finite or negative.
    pub fn slope_length_for_height(&self, height: f64) -> Result<f64, InclinedPlaneError> {
        if !height.is_finite() || height < 0.0 {
            return Err(InclinedPlaneError::bad_parameter(
                "height",
                format!("must be finite and >= 0, got {height}"),
            ));
        }
        Ok(height / self.angle_rad.sin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn rejects_flat_and_vertical_and_nonfinite() {
        assert!(IdealRamp::from_angle(0.0).is_err());
        assert!(IdealRamp::from_angle(std::f64::consts::FRAC_PI_2).is_err());
        assert!(IdealRamp::from_angle(-0.1).is_err());
        assert!(IdealRamp::from_angle(f64::NAN).is_err());
        assert!(IdealRamp::from_angle(f64::INFINITY).is_err());
    }

    #[test]
    fn from_rise_run_validates_inputs() {
        assert!(IdealRamp::from_rise_run(0.0, 1.0).is_err());
        assert!(IdealRamp::from_rise_run(1.0, 0.0).is_err());
        assert!(IdealRamp::from_rise_run(-1.0, 1.0).is_err());
        assert!(IdealRamp::from_rise_run(1.0, f64::NAN).is_err());
    }

    #[test]
    fn rise_run_recovers_angle() {
        // rise == run -> 45 degrees.
        let ramp = IdealRamp::from_rise_run(2.0, 2.0).unwrap();
        assert!((ramp.angle_rad() - std::f64::consts::FRAC_PI_4).abs() < EPS);
        assert!((ramp.angle_deg() - 45.0).abs() < 1e-9);
    }

    #[test]
    fn mechanical_advantage_is_one_over_sin() {
        // 30 degrees: sin = 0.5, so MA = 2 exactly.
        let ramp = IdealRamp::from_angle(std::f64::consts::FRAC_PI_6).unwrap();
        assert!((ramp.ideal_mechanical_advantage() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn steeper_ramp_has_less_mechanical_advantage() {
        let shallow = IdealRamp::from_angle(0.2).unwrap();
        let steep = IdealRamp::from_angle(1.0).unwrap();
        assert!(steep.angle_rad() > shallow.angle_rad());
        assert!(steep.ideal_mechanical_advantage() < shallow.ideal_mechanical_advantage());
    }

    #[test]
    fn vertical_limit_approaches_unity_ma() {
        // As theta -> pi/2, MA -> 1. Check just below vertical.
        let almost_vertical = IdealRamp::from_angle(std::f64::consts::FRAC_PI_2 - 1e-6).unwrap();
        assert!((almost_vertical.ideal_mechanical_advantage() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn slope_length_matches_height_over_sin() {
        // 30 degrees, raise by 1.0: slope length = 1 / 0.5 = 2.0.
        let ramp = IdealRamp::from_angle(std::f64::consts::FRAC_PI_6).unwrap();
        let len = ramp.slope_length_for_height(1.0).unwrap();
        assert!((len - 2.0).abs() < 1e-12);
        assert!(ramp.slope_length_for_height(-1.0).is_err());
    }
}
