//! **Occupant / vehicle survivability** — a simple acceleration-injury screen.
//!
//! Blast and impact protection is ultimately about keeping *people* alive, so
//! the final survivability readout compares the acceleration the occupant
//! experiences against a human tolerance threshold. This is the same physics as
//! automotive crash safety: a protective structure that survives the load is
//! only useful if the *occupant* also survives the resulting accelerations.
//!
//! ## Model
//!
//! The screen is deliberately simple and transparent:
//!
//! - **Peak acceleration in g** — convert the peak structural/seat acceleration
//!   `a` (m/s²) to gravities `g_peak = a / g₀` (`g₀ = 9.80665 m/s²`).
//! - **Tolerance margin** — compare `g_peak` to a supplied **tolerance limit**
//!   `g_tol` (the survivable peak for the relevant duration and direction). The
//!   margin is `g_tol / g_peak`: `≥ 1` is survivable, `< 1` exceeds tolerance.
//!
//! Human acceleration tolerance is duration- and direction-dependent (Eiband,
//! NASA Memo 5-19-59E, 1959; the basis of the Eiband whole-body tolerance
//! curves). Rather than hard-code a single number, the caller passes the
//! `g_tol` appropriate to the pulse duration and loading axis. For very short
//! pulses, tolerance is far higher than for sustained loads — that is the
//! caller's `g_tol` choice, kept out of this screen on purpose.
//!
//! This is a research/educational survivability *screen*, validation-pending —
//! not an injury-prediction or medical-certification model. A "survivable"
//! result means only that the peak acceleration is below the supplied
//! whole-body tolerance, not that no injury occurs.

use crate::error::SurvivabilityError;
use serde::{Deserialize, Serialize};

/// Standard gravity `g₀` (m/s²).
pub const G0: f64 = 9.806_65;

/// The occupant acceleration-injury screen result.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OccupantAssessment {
    /// Peak acceleration in gravities, `g_peak = a / g₀`.
    pub peak_g: f64,
    /// The whole-body tolerance limit `g_tol` (g) supplied for this duration /
    /// direction.
    pub tolerance_g: f64,
    /// Survivability margin `g_tol / g_peak`. `≥ 1` ⇒ within tolerance.
    pub margin: f64,
    /// `true` when `g_peak ≤ g_tol` (the occupant is within the supplied
    /// whole-body tolerance).
    pub survivable: bool,
}

/// Assess occupant survivability from a peak acceleration `peak_accel_m_s2`
/// (m/s²) against a whole-body tolerance limit `tolerance_g` (g).
///
/// Typically the peak acceleration comes straight from a structural-response
/// solve — [`crate::response::SdofResponse::peak_acceleration_m_s2`].
///
/// # Errors
///
/// [`SurvivabilityError::InvalidParameter`] if the peak acceleration is not
/// finite-and-non-negative, or the tolerance is not finite-and-positive.
pub fn assess_occupant(
    peak_accel_m_s2: f64,
    tolerance_g: f64,
) -> Result<OccupantAssessment, SurvivabilityError> {
    if !(peak_accel_m_s2.is_finite() && peak_accel_m_s2 >= 0.0) {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "peak acceleration must be finite and >= 0, got {peak_accel_m_s2}"
        )));
    }
    if !(tolerance_g.is_finite() && tolerance_g > 0.0) {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "tolerance must be finite and > 0 g, got {tolerance_g}"
        )));
    }

    let peak_g = peak_accel_m_s2 / G0;
    // peak_g may be 0 (no load) ⇒ margin is +∞, which is correctly "survivable".
    let margin = if peak_g > 0.0 {
        tolerance_g / peak_g
    } else {
        f64::INFINITY
    };
    let survivable = peak_g <= tolerance_g;

    Ok(OccupantAssessment {
        peak_g,
        tolerance_g,
        margin,
        survivable,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peak_g_conversion() {
        let a = assess_occupant(98.0665, 50.0).unwrap();
        assert!((a.peak_g - 10.0).abs() < 1e-9); // 98.0665 / 9.80665 = 10 g
    }

    #[test]
    fn survivable_when_below_tolerance() {
        let a = assess_occupant(20.0 * G0, 40.0).unwrap();
        assert!(a.survivable);
        assert!(a.margin > 1.0);
        assert!((a.margin - 2.0).abs() < 1e-9);
    }

    #[test]
    fn not_survivable_when_above_tolerance() {
        let a = assess_occupant(60.0 * G0, 40.0).unwrap();
        assert!(!a.survivable);
        assert!(a.margin < 1.0);
    }

    #[test]
    fn margin_decreases_as_acceleration_grows() {
        // PIN: higher peak acceleration ⇒ smaller (worse) margin.
        let mut prev = f64::INFINITY;
        for g in [5.0, 10.0, 20.0, 40.0, 80.0] {
            let a = assess_occupant(g * G0, 50.0).unwrap();
            assert!(a.margin < prev, "margin not decreasing at {g} g");
            prev = a.margin;
        }
    }

    #[test]
    fn zero_acceleration_is_survivable() {
        let a = assess_occupant(0.0, 30.0).unwrap();
        assert!(a.survivable);
        assert!(a.peak_g == 0.0);
        assert!(a.margin.is_infinite());
    }

    #[test]
    fn degenerate_inputs_error_not_panic() {
        assert!(assess_occupant(-1.0, 30.0).is_err()); // negative accel
        assert!(assess_occupant(f64::NAN, 30.0).is_err());
        assert!(assess_occupant(10.0, 0.0).is_err()); // zero tolerance
        assert!(assess_occupant(10.0, -5.0).is_err());
        assert!(assess_occupant(f64::INFINITY, 30.0).is_err());
    }

    #[test]
    fn serde_round_trip() {
        let a = assess_occupant(30.0 * G0, 50.0).unwrap();
        let json = serde_json::to_string(&a).unwrap();
        let back: OccupantAssessment = serde_json::from_str(&json).unwrap();
        // Compare with a tolerance: JSON text round-trips f64 to within ~1 ULP,
        // so an exact `assert_eq!` on the derived `margin` is too strict.
        assert!((a.peak_g - back.peak_g).abs() < 1e-9);
        assert!((a.tolerance_g - back.tolerance_g).abs() < 1e-9);
        assert!((a.margin - back.margin).abs() < 1e-9);
        assert_eq!(a.survivable, back.survivable);
    }
}
