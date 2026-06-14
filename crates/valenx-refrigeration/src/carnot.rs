//! Reversible (Carnot) coefficient-of-performance limits.
//!
//! No refrigerator or heat pump operating between a cold reservoir at
//! absolute temperature `Tc` and a hot reservoir at absolute temperature
//! `Th` (with `Th > Tc`) can exceed the reversible Carnot device that
//! works between the same two reservoirs. The Carnot limits are
//!
//! ```text
//! COP_cool,Carnot = Tc / (Th - Tc)
//! COP_heat,Carnot = Th / (Th - Tc)
//! ```
//!
//! and they differ by exactly one, mirroring the energy-balance identity
//! for real cycles. Both diverge as the temperature lift `Th - Tc`
//! shrinks to zero (free heat pumping at zero lift) and both fall as the
//! lift grows, so **a larger temperature lift always lowers the
//! achievable COP**.
//!
//! Temperatures here are *absolute* (kelvin); a non-positive temperature
//! or a non-positive lift is rejected.

use crate::error::{RefrigError, Result};

/// Carnot (reversible) cooling COP between two reservoirs,
/// `Tc / (Th - Tc)`.
///
/// `t_cold` and `t_hot` are absolute temperatures (kelvin) of the cold
/// (refrigerated) and hot (heat-rejection) reservoirs respectively.
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if either temperature is not a
/// finite, strictly positive value, and [`RefrigError::InvalidLift`] if
/// `t_hot` does not strictly exceed `t_cold`.
pub fn carnot_cop_cool(t_cold: f64, t_hot: f64) -> Result<f64> {
    let (tc, th) = validated_lift(t_cold, t_hot)?;
    Ok(tc / (th - tc))
}

/// Carnot (reversible) heating COP between two reservoirs,
/// `Th / (Th - Tc)`.
///
/// `t_cold` and `t_hot` are absolute temperatures (kelvin).
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if either temperature is not a
/// finite, strictly positive value, and [`RefrigError::InvalidLift`] if
/// `t_hot` does not strictly exceed `t_cold`.
pub fn carnot_cop_heat(t_cold: f64, t_hot: f64) -> Result<f64> {
    let (tc, th) = validated_lift(t_cold, t_hot)?;
    Ok(th / (th - tc))
}

/// Second-law (exergetic) efficiency of a real cooling cycle: the ratio
/// of its actual cooling COP to the Carnot limit for the same reservoir
/// temperatures, `eta_II = COP_cool / COP_cool,Carnot`.
///
/// The result lies in `(0, 1]`: it equals one only for a reversible
/// cycle and is strictly less than one for any real cycle. A value above
/// one signals an input inconsistency (a real cycle cannot beat Carnot).
///
/// `cop_cool` is the actual cooling coefficient of performance;
/// `t_cold` and `t_hot` are the absolute reservoir temperatures.
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if `cop_cool` is negative or
/// non-finite, propagates the temperature-validation errors of
/// [`carnot_cop_cool`], and returns
/// [`RefrigError::InconsistentCycle`] if the actual COP exceeds the
/// Carnot limit (which would violate the second law).
pub fn second_law_efficiency_cool(cop_cool: f64, t_cold: f64, t_hot: f64) -> Result<f64> {
    let cop = RefrigError::require_non_negative("cop_cool", cop_cool)?;
    let carnot = carnot_cop_cool(t_cold, t_hot)?;
    let eta = cop / carnot;
    // Allow a tiny epsilon so an exactly-Carnot input does not trip on
    // round-off, but reject genuine super-Carnot performance.
    if eta > 1.0 + 1e-9 {
        return Err(RefrigError::InconsistentCycle(
            "actual cooling COP exceeds the Carnot limit (violates the second law)",
        ));
    }
    Ok(eta)
}

/// Second-law (exergetic) efficiency of a real heat-pump cycle: the
/// ratio of its actual heating COP to the Carnot heating limit,
/// `eta_II = COP_heat / COP_heat,Carnot`.
///
/// Like [`second_law_efficiency_cool`], the result lies in `(0, 1]` and
/// a value above one signals an input that violates the second law.
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if `cop_heat` is negative or
/// non-finite, propagates the temperature-validation errors of
/// [`carnot_cop_heat`], and returns
/// [`RefrigError::InconsistentCycle`] if the actual COP exceeds the
/// Carnot limit.
pub fn second_law_efficiency_heat(cop_heat: f64, t_cold: f64, t_hot: f64) -> Result<f64> {
    let cop = RefrigError::require_non_negative("cop_heat", cop_heat)?;
    let carnot = carnot_cop_heat(t_cold, t_hot)?;
    let eta = cop / carnot;
    if eta > 1.0 + 1e-9 {
        return Err(RefrigError::InconsistentCycle(
            "actual heating COP exceeds the Carnot limit (violates the second law)",
        ));
    }
    Ok(eta)
}

/// Validate a temperature pair and return `(t_cold, t_hot)` on success.
///
/// Both must be finite and strictly positive, and `t_hot` must strictly
/// exceed `t_cold` so the lift `Th - Tc` is positive.
fn validated_lift(t_cold: f64, t_hot: f64) -> Result<(f64, f64)> {
    let tc = RefrigError::require_positive("t_cold", t_cold)?;
    let th = RefrigError::require_positive("t_hot", t_hot)?;
    if th <= tc {
        return Err(RefrigError::InvalidLift {
            t_hot: th,
            t_cold: tc,
        });
    }
    Ok((tc, th))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in this module.
    const EPS: f64 = 1e-9;

    #[test]
    fn carnot_cool_known_value() {
        // Tc = 250 K, Th = 300 K → 250 / 50 = 5.
        let cop = carnot_cop_cool(250.0, 300.0).unwrap();
        assert!((cop - 5.0).abs() < EPS, "got {cop}");
    }

    #[test]
    fn carnot_heat_known_value() {
        // Tc = 250 K, Th = 300 K → 300 / 50 = 6.
        let cop = carnot_cop_heat(250.0, 300.0).unwrap();
        assert!((cop - 6.0).abs() < EPS, "got {cop}");
    }

    #[test]
    fn carnot_heating_is_cooling_plus_one() {
        // The reversible limits obey the same +1 identity as real cycles.
        let tc = 268.15;
        let th = 313.15;
        let cc = carnot_cop_cool(tc, th).unwrap();
        let ch = carnot_cop_heat(tc, th).unwrap();
        assert!((ch - (cc + 1.0)).abs() < EPS, "cc={cc} ch={ch}");
    }

    #[test]
    fn larger_lift_lowers_cop() {
        // Fix the cold side, raise the hot side: the lift grows and the
        // cooling COP must fall monotonically.
        let tc = 270.0;
        let cop_small_lift = carnot_cop_cool(tc, 290.0).unwrap();
        let cop_mid_lift = carnot_cop_cool(tc, 310.0).unwrap();
        let cop_big_lift = carnot_cop_cool(tc, 330.0).unwrap();
        assert!(
            cop_small_lift > cop_mid_lift && cop_mid_lift > cop_big_lift,
            "{cop_small_lift} > {cop_mid_lift} > {cop_big_lift}"
        );
    }

    #[test]
    fn carnot_bounds_a_real_cycle() {
        // A real cooling COP below the Carnot limit yields a second-law
        // efficiency strictly inside (0, 1).
        let tc = 250.0;
        let th = 300.0;
        let carnot = carnot_cop_cool(tc, th).unwrap();
        let real_cop = 0.4 * carnot; // a plausible 40 % of reversible
        assert!(real_cop < carnot);

        let eta = second_law_efficiency_cool(real_cop, tc, th).unwrap();
        assert!((eta - 0.4).abs() < EPS, "eta={eta}");
        assert!(eta > 0.0 && eta < 1.0);
    }

    #[test]
    fn reversible_cycle_has_unit_efficiency() {
        let tc = 255.0;
        let th = 305.0;
        let carnot = carnot_cop_cool(tc, th).unwrap();
        let eta = second_law_efficiency_cool(carnot, tc, th).unwrap();
        assert!((eta - 1.0).abs() < EPS, "eta={eta}");
    }

    #[test]
    fn super_carnot_input_is_rejected() {
        let tc = 250.0;
        let th = 300.0;
        let carnot = carnot_cop_cool(tc, th).unwrap();
        // 10 % above the reversible limit is impossible.
        let bogus = carnot * 1.1;
        assert!(matches!(
            second_law_efficiency_cool(bogus, tc, th),
            Err(RefrigError::InconsistentCycle(_))
        ));
    }

    #[test]
    fn invalid_lift_is_rejected() {
        // Th must exceed Tc.
        assert!(matches!(
            carnot_cop_cool(300.0, 300.0),
            Err(RefrigError::InvalidLift { .. })
        ));
        assert!(matches!(
            carnot_cop_heat(310.0, 300.0),
            Err(RefrigError::InvalidLift { .. })
        ));
        // Non-positive absolute temperature.
        assert!(matches!(
            carnot_cop_cool(0.0, 300.0),
            Err(RefrigError::BadParameter { .. })
        ));
        assert!(matches!(
            carnot_cop_cool(-5.0, 300.0),
            Err(RefrigError::BadParameter { .. })
        ));
    }
}
