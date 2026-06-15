//! Drive-power models for a belt conveyor.
//!
//! # Model
//!
//! Two complementary closed-form views of conveyor power, both at
//! steady state:
//!
//! Lift power — the rate of work done raising material against gravity
//! over a vertical rise `h` (m):
//!
//! ```text
//! P_lift = mdot * g * h          [W]
//! ```
//!
//! where `mdot` is the mass flow rate (kg/s) and `g` is gravitational
//! acceleration. A horizontal conveyor (`h = 0`) does zero lift work; an
//! inclined conveyor adds exactly this term.
//!
//! Drive power from belt tension — the mechanical power the drive pulley
//! must deliver to overcome the effective belt tension `F` (N) at belt
//! speed `v` (m/s):
//!
//! ```text
//! P = F * v          [W]
//! ```
//!
//! For a pure lift with no friction the effective tension is the weight
//! component along the incline, and the two formulas agree (see the
//! tests). When friction is present, `F` is the sum of the friction
//! resistance and the lift component, and `P = F * v` is the more
//! general statement.
//!
//! # Honest scope
//!
//! These are first-principles relations. They do not model the detailed
//! resistance build-up of a real installation — idler rolling
//! resistance, belt flexure, skirtboard drag, pulley wrap, take-up and
//! drive efficiency — which standards such as CEMA and ISO 5048 fold
//! into the effective tension `F`. Supply your own `F` (or friction
//! coefficient) to use those; this crate keeps the algebra transparent.

use crate::belt::G;
use crate::error::ConveyorError;
use serde::{Deserialize, Serialize};

/// Lift power `P_lift = mdot * g * h` in watts, using standard gravity.
///
/// `mass_flow` is in kg/s and `lift_height` in m. A zero lift height is
/// allowed and yields zero power (a horizontal conveyor does no lift
/// work); a negative height is rejected.
///
/// # Errors
///
/// Returns [`ConveyorError`] if `mass_flow` is non-finite or not
/// strictly positive, or if `lift_height` is non-finite or negative.
pub fn lift_power(mass_flow: f64, lift_height: f64) -> Result<f64, ConveyorError> {
    lift_power_g(mass_flow, lift_height, G)
}

/// Lift power with a caller-supplied gravitational acceleration `g`
/// (m/s^2), for non-Earth or sensitivity studies.
///
/// # Errors
///
/// Returns [`ConveyorError`] if `mass_flow` or `g` is non-finite or not
/// strictly positive, or if `lift_height` is non-finite or negative.
pub fn lift_power_g(mass_flow: f64, lift_height: f64, g: f64) -> Result<f64, ConveyorError> {
    let mass_flow = ConveyorError::require_positive("mass_flow", mass_flow)?;
    let g = ConveyorError::require_positive("g", g)?;
    let lift_height = ConveyorError::require_range("lift_height", lift_height, 0.0, f64::INFINITY)?;
    Ok(mass_flow * g * lift_height)
}

/// Drive power `P = F * v` in watts from an effective belt tension `F`
/// (N) and belt speed `v` (m/s).
///
/// # Errors
///
/// Returns [`ConveyorError`] if `tension` is non-finite or negative, or
/// if `speed` is non-finite or not strictly positive.
pub fn drive_power(tension: f64, speed: f64) -> Result<f64, ConveyorError> {
    let tension = ConveyorError::require_range("tension", tension, 0.0, f64::INFINITY)?;
    let speed = ConveyorError::require_positive("speed", speed)?;
    Ok(tension * speed)
}

/// Vertical rise of an incline of slope `angle_rad` (radians) over a
/// belt path of length `length` (m): `h = length * sin(angle)`.
///
/// Useful for turning an inclined-conveyor geometry into the `lift_height`
/// argument of [`lift_power`].
///
/// # Errors
///
/// Returns [`ConveyorError`] if `length` is non-finite or not strictly
/// positive, or if `angle_rad` is non-finite or outside
/// `[0, PI/2]` radians (0 = horizontal, PI/2 = vertical).
pub fn incline_rise(length: f64, angle_rad: f64) -> Result<f64, ConveyorError> {
    let length = ConveyorError::require_positive("length", length)?;
    let angle_rad =
        ConveyorError::require_range("angle_rad", angle_rad, 0.0, std::f64::consts::FRAC_PI_2)?;
    Ok(length * angle_rad.sin())
}

/// A resolved conveyor power breakdown returned by [`PowerBreakdown::solve`].
///
/// Holds the effective tension actually used and the resulting drive and
/// lift powers, so a caller can show how much of the demand is lift work
/// versus friction work.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PowerBreakdown {
    /// Effective belt tension used for the drive-power figure, N.
    pub tension: f64,
    /// Lift power `mdot * g * h`, W.
    pub lift: f64,
    /// Total drive power `F * v`, W.
    pub total: f64,
}

impl PowerBreakdown {
    /// Resolve the power demand of a horizontal-plus-incline conveyor.
    ///
    /// The effective tension is built from a horizontal friction
    /// resistance `friction_force` (N) plus the along-incline weight
    /// component of the conveyed load, `mdot * g * sin(angle) * (L / v)`
    /// expressed as a force. To keep the model transparent and
    /// self-consistent with [`lift_power`], the lift contribution to the
    /// tension is taken as `m_on_belt * g * sin(angle)` where the mass
    /// resident on the belt is `m_on_belt = mdot * L / v` (mass flow
    /// times residence time). The drive power is then `F * v`, which
    /// equals `friction_force * v + mdot * g * h` with `h = L sin(angle)`.
    ///
    /// # Errors
    ///
    /// Returns [`ConveyorError`] if any of `mass_flow`, `speed` or
    /// `length` is non-finite or not strictly positive, if
    /// `friction_force` is non-finite or negative, or if `angle_rad` is
    /// non-finite or outside `[0, PI/2]`.
    pub fn solve(
        mass_flow: f64,
        speed: f64,
        length: f64,
        angle_rad: f64,
        friction_force: f64,
    ) -> Result<Self, ConveyorError> {
        let mass_flow = ConveyorError::require_positive("mass_flow", mass_flow)?;
        let speed = ConveyorError::require_positive("speed", speed)?;
        let length = ConveyorError::require_positive("length", length)?;
        let friction_force =
            ConveyorError::require_range("friction_force", friction_force, 0.0, f64::INFINITY)?;
        let h = incline_rise(length, angle_rad)?;

        let lift = mass_flow * G * h;
        // P = F * v, and the lift part of F*v is exactly `lift`, so the
        // effective tension is the friction force plus lift / v.
        let tension = friction_force + lift / speed;
        let total = tension * speed;
        Ok(Self {
            tension,
            lift,
            total,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn lift_power_matches_mdot_g_h() {
        // mdot = 600 kg/s, h = 10 m -> 600 * 9.80665 * 10 = 58839.9 W.
        let p = lift_power(600.0, 10.0).expect("valid");
        assert!((p - 600.0 * G * 10.0).abs() < 1e-6);
        assert!((p - 58_839.9).abs() < 1e-1);
    }

    #[test]
    fn horizontal_conveyor_has_zero_lift_power() {
        let p = lift_power(600.0, 0.0).expect("valid");
        assert!(p.abs() < EPS);
    }

    #[test]
    fn lift_power_scales_linearly_with_height() {
        let p1 = lift_power(600.0, 5.0).expect("valid");
        let p2 = lift_power(600.0, 10.0).expect("valid");
        assert!((p2 - 2.0 * p1).abs() < 1e-6);
    }

    #[test]
    fn drive_power_matches_f_times_v() {
        // F = 5000 N, v = 2.5 m/s -> 12500 W.
        let p = drive_power(5000.0, 2.5).expect("valid");
        assert!((p - 12_500.0).abs() < 1e-9);
    }

    #[test]
    fn drive_power_scales_linearly_with_speed() {
        let p1 = drive_power(5000.0, 2.0).expect("valid");
        let p2 = drive_power(5000.0, 4.0).expect("valid");
        assert!((p2 - 2.0 * p1).abs() < 1e-9);
    }

    #[test]
    fn incline_rise_matches_l_sin_theta() {
        // 30 deg over 100 m -> 100 * 0.5 = 50 m.
        let h = incline_rise(100.0, std::f64::consts::FRAC_PI_6).expect("valid");
        assert!((h - 50.0).abs() < 1e-9);
    }

    #[test]
    fn horizontal_incline_has_zero_rise() {
        let h = incline_rise(100.0, 0.0).expect("valid");
        assert!(h.abs() < EPS);
    }

    #[test]
    fn vertical_incline_rise_equals_length() {
        let h = incline_rise(100.0, std::f64::consts::FRAC_PI_2).expect("valid");
        assert!((h - 100.0).abs() < 1e-9);
    }

    #[test]
    fn lift_via_drive_power_agrees_for_pure_lift() {
        // A pure (frictionless) lift: the along-incline weight component
        // of the resident mass equals mdot*g*sin(theta)*(L/v), so the
        // tension F = mdot*g*h / v and F*v == mdot*g*h.
        let mdot = 600.0;
        let v = 2.0;
        let length = 100.0;
        let angle = std::f64::consts::FRAC_PI_6; // 30 deg, h = 50 m.
        let h = incline_rise(length, angle).expect("valid");
        let p_lift = lift_power(mdot, h).expect("valid");

        let bd = PowerBreakdown::solve(mdot, v, length, angle, 0.0).expect("valid");
        // With zero friction, total drive power equals the lift power.
        assert!((bd.total - p_lift).abs() < 1e-6);
        assert!((bd.lift - p_lift).abs() < 1e-6);
        // And F * v reproduces the total.
        assert!((bd.tension * v - bd.total).abs() < 1e-6);
    }

    #[test]
    fn incline_adds_lift_over_horizontal() {
        // Same friction, but an incline must demand strictly more power
        // than the horizontal case (extra term is the lift power).
        let flat = PowerBreakdown::solve(600.0, 2.0, 100.0, 0.0, 5000.0).expect("valid");
        let up = PowerBreakdown::solve(600.0, 2.0, 100.0, std::f64::consts::FRAC_PI_6, 5000.0)
            .expect("valid");
        assert!(up.total > flat.total);
        // The flat case has no lift; its drive power is pure friction*v.
        assert!(flat.lift.abs() < EPS);
        assert!((flat.total - 5000.0 * 2.0).abs() < 1e-6);
        // The incline's extra power equals its lift power exactly.
        let extra = up.total - flat.total;
        assert!((extra - up.lift).abs() < 1e-6);
    }

    #[test]
    fn alt_gravity_changes_lift_power() {
        // Lunar g ~ 1.62 m/s^2 gives a smaller lift power than Earth.
        let earth = lift_power(600.0, 10.0).expect("valid");
        let moon = lift_power_g(600.0, 10.0, 1.62).expect("valid");
        assert!(moon < earth);
        assert!((moon - 600.0 * 1.62 * 10.0).abs() < 1e-6);
    }

    #[test]
    fn power_functions_reject_bad_inputs() {
        assert!(lift_power(0.0, 10.0).is_err());
        assert!(lift_power(600.0, -1.0).is_err());
        assert!(drive_power(-1.0, 2.0).is_err());
        assert!(drive_power(5000.0, 0.0).is_err());
        assert!(incline_rise(100.0, -0.1).is_err());
        assert!(incline_rise(100.0, std::f64::consts::PI).is_err());
        assert!(PowerBreakdown::solve(600.0, 2.0, 100.0, 0.0, -1.0).is_err());
    }
}
