//! Shunt-capacitor power-factor correction.
//!
//! A lagging (inductive) load draws reactive power `Q1 = P * tan(phi1)`.
//! Adding a shunt capacitor in parallel supplies leading reactive power
//! that cancels part of that draw, lifting the power factor toward unity
//! without changing the real power `P`. To reach a better target angle
//! `phi2 < phi1` the capacitor must provide
//!
//! ```text
//!   Qc = P * (tan(phi1) - tan(phi2))     [var]
//! ```
//!
//! Because the real power is unchanged, the new reactive power is
//! `Q2 = P * tan(phi2) = Q1 - Qc`, and the apparent power drops from
//! `S1 = P / cos(phi1)` to `S2 = P / cos(phi2)`. The lower apparent
//! power (and hence lower current) is the practical motivation for
//! correction.
//!
//! ## Scope reminder
//!
//! This is the textbook single-frequency model: it sizes an *ideal*
//! capacitor's reactive rating and assumes a purely sinusoidal lagging
//! load. It does not account for harmonics, over-correction resonance,
//! switching transients, or capacitor tolerance/derating, so it is not
//! a substitute for engineered capacitor-bank design.

use crate::error::PowerError;
use crate::triangle::{Phase, PowerTriangle};
use serde::{Deserialize, Serialize};

/// The result of a power-factor correction calculation.
///
/// Captures the before/after reactive powers and the required capacitor
/// reactive rating, so the caller can both size the capacitor and verify
/// the outcome.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Correction {
    /// Real power `P`, unchanged by correction, in watts (W).
    pub real_w: f64,
    /// Reactive power before correction, `Q1 = P * tan(phi1)`, in var.
    pub reactive_before_var: f64,
    /// Reactive power after correction, `Q2 = P * tan(phi2)`, in var.
    pub reactive_after_var: f64,
    /// Capacitor reactive rating required,
    /// `Qc = Q1 - Q2 = P * (tan(phi1) - tan(phi2))`, in var. Always
    /// `> 0` for a genuine improvement.
    pub capacitor_var: f64,
    /// Power factor before correction (the starting `cos(phi1)`).
    pub power_factor_before: f64,
    /// Power factor after correction (the target `cos(phi2)`).
    pub power_factor_after: f64,
    /// Apparent power before correction, `S1 = P / cos(phi1)`, in VA.
    pub apparent_before_va: f64,
    /// Apparent power after correction, `S2 = P / cos(phi2)`, in VA.
    pub apparent_after_va: f64,
}

impl Correction {
    /// Size a shunt capacitor to raise a lagging load's power factor.
    ///
    /// Given the real power `real_w` and the present and target power
    /// factors (both lagging), returns the capacitor reactive rating and
    /// the full before/after picture.
    ///
    /// `real_w` must be strictly positive (there is nothing to correct
    /// for a load with no real power). Both power factors must lie in
    /// `(0, 1]`, and `power_factor_target` must be strictly greater than
    /// `power_factor_present` — correction only ever raises the power
    /// factor.
    ///
    /// # Errors
    ///
    /// Returns [`PowerError::NonPositive`] if `real_w` is not positive,
    /// [`PowerError::PowerFactorOutOfRange`] if either power factor is
    /// outside `[0, 1]` (or is exactly `0`, where `tan(phi)` diverges),
    /// and [`PowerError::NoCorrectionNeeded`] if the target does not
    /// improve on the present power factor.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_powerfactor::Correction;
    /// // Raise a 10 kW, 0.7 PF load to 0.95.
    /// let c = Correction::for_target_pf(10_000.0, 0.7, 0.95).unwrap();
    /// assert!(c.capacitor_var > 0.0);
    /// assert!(c.reactive_after_var < c.reactive_before_var);
    /// assert!((c.power_factor_after - 0.95).abs() < 1e-12);
    /// ```
    pub fn for_target_pf(
        real_w: f64,
        power_factor_present: f64,
        power_factor_target: f64,
    ) -> Result<Correction, PowerError> {
        let real_w = PowerError::positive("real_w", real_w)?;
        let pf1 = PowerError::power_factor("power_factor_present", power_factor_present)?;
        let pf2 = PowerError::power_factor("power_factor_target", power_factor_target)?;
        // tan(phi) = sin/cos = sqrt(1 - PF^2) / PF diverges at PF = 0.
        if pf1 == 0.0 {
            return Err(PowerError::PowerFactorOutOfRange { value: pf1 });
        }
        if pf2 == 0.0 {
            return Err(PowerError::PowerFactorOutOfRange { value: pf2 });
        }
        if pf2 <= pf1 {
            return Err(PowerError::NoCorrectionNeeded {
                present: pf1,
                target: pf2,
            });
        }

        let tan_phi1 = tan_from_pf(pf1);
        let tan_phi2 = tan_from_pf(pf2);
        let reactive_before_var = real_w * tan_phi1;
        let reactive_after_var = real_w * tan_phi2;
        let capacitor_var = reactive_before_var - reactive_after_var;
        Ok(Correction {
            real_w,
            reactive_before_var,
            reactive_after_var,
            capacitor_var,
            power_factor_before: pf1,
            power_factor_after: pf2,
            apparent_before_va: real_w / pf1,
            apparent_after_va: real_w / pf2,
        })
    }

    /// Size a shunt capacitor from an existing power triangle and a
    /// target power factor.
    ///
    /// A convenience wrapper over [`Correction::for_target_pf`] that
    /// pulls the real power and present power factor out of an already
    /// resolved [`PowerTriangle`]. The triangle must describe a lagging
    /// load: correcting a leading load with a (leading) capacitor would
    /// only make it worse, so [`Phase::Leading`] is rejected.
    ///
    /// # Errors
    ///
    /// In addition to the errors from [`Correction::for_target_pf`],
    /// returns [`PowerError::NoCorrectionNeeded`] when `triangle` is not
    /// a lagging load (a leading or unity load cannot be improved by a
    /// shunt capacitor).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_powerfactor::{Correction, Phase, PowerTriangle};
    /// let load = PowerTriangle::from_vi_pf(400.0, 50.0, 0.75, Phase::Lagging).unwrap();
    /// let c = Correction::for_triangle(&load, 0.98).unwrap();
    /// assert!(c.capacitor_var > 0.0);
    /// ```
    pub fn for_triangle(
        triangle: &PowerTriangle,
        power_factor_target: f64,
    ) -> Result<Correction, PowerError> {
        if triangle.phase != Phase::Lagging {
            return Err(PowerError::NoCorrectionNeeded {
                present: triangle.power_factor,
                target: power_factor_target,
            });
        }
        Correction::for_target_pf(triangle.real_w, triangle.power_factor, power_factor_target)
    }
}

/// Tangent of the phase angle for a given lagging power factor.
///
/// For `PF = cos(phi)` with `phi` in `(0, pi/2]`, the tangent is
/// `tan(phi) = sin(phi) / cos(phi) = sqrt(1 - PF^2) / PF`. The caller
/// must ensure `power_factor` is in `(0, 1]`; the radicand is clamped to
/// zero to absorb round-off when `power_factor == 1`.
fn tan_from_pf(power_factor: f64) -> f64 {
    let sin_phi = (1.0 - power_factor * power_factor).max(0.0).sqrt();
    sin_phi / power_factor
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in tests.
    const EPS: f64 = 1e-9;

    #[test]
    fn tan_from_pf_known_angles() {
        // PF = cos 60 = 0.5 => tan 60 = sqrt(3).
        assert!((tan_from_pf(0.5) - (3.0_f64).sqrt()).abs() < EPS);
        // PF = cos 45 => tan 45 = 1.
        let pf45 = (45_f64).to_radians().cos();
        assert!((tan_from_pf(pf45) - 1.0).abs() < EPS);
        // PF = 1 => tan 0 = 0.
        assert!(tan_from_pf(1.0).abs() < EPS);
    }

    #[test]
    fn capacitor_matches_closed_form() {
        // Qc = P * (tan phi1 - tan phi2). Cross-check against an
        // independent hand computation for P = 10 kW, 0.8 -> 0.95.
        let p = 10_000.0;
        let pf1 = 0.8;
        let pf2 = 0.95;
        let c = Correction::for_target_pf(p, pf1, pf2).unwrap();
        let expected_qc = p * (tan_from_pf(pf1) - tan_from_pf(pf2));
        assert!((c.capacitor_var - expected_qc).abs() < 1e-6);
        // Ground-truth numbers (textbook): tan(acos 0.8)=0.75,
        // tan(acos 0.95)=0.328684. Qc ~= 10000*(0.75-0.328684)=4213.16.
        assert!((c.capacitor_var - 4213.16).abs() < 1.0);
    }

    #[test]
    fn correction_lowers_reactive_power() {
        // VALIDATE: correction lowers Q. Q2 < Q1 strictly for pf2 > pf1.
        let c = Correction::for_target_pf(5000.0, 0.6, 0.9).unwrap();
        assert!(c.reactive_after_var < c.reactive_before_var);
        assert!(c.capacitor_var > 0.0);
        // Qc is exactly the reduction in Q.
        assert!((c.capacitor_var - (c.reactive_before_var - c.reactive_after_var)).abs() < EPS);
    }

    #[test]
    fn real_power_is_preserved() {
        // Correction must not change P.
        let p = 7500.0;
        let c = Correction::for_target_pf(p, 0.65, 0.99).unwrap();
        assert!((c.real_w - p).abs() < EPS);
    }

    #[test]
    fn apparent_power_drops_after_correction() {
        // S2 = P/pf2 < S1 = P/pf1 because pf2 > pf1.
        let c = Correction::for_target_pf(12_000.0, 0.7, 0.95).unwrap();
        assert!(c.apparent_after_va < c.apparent_before_va);
        assert!((c.apparent_before_va - 12_000.0 / 0.7).abs() < 1e-6);
        assert!((c.apparent_after_va - 12_000.0 / 0.95).abs() < 1e-6);
    }

    #[test]
    fn corrected_triangle_has_target_pf() {
        // Rebuild the post-correction triangle from P and Q2 and confirm
        // its power factor equals the requested target.
        let c = Correction::for_target_pf(9000.0, 0.72, 0.96).unwrap();
        let after = PowerTriangle::from_p_q(c.real_w, c.reactive_after_var).unwrap();
        assert!((after.power_factor - 0.96).abs() < 1e-9);
    }

    #[test]
    fn unity_target_drives_reactive_to_zero() {
        // Correcting all the way to PF = 1 cancels all reactive power:
        // Qc == Q1 and Q2 == 0.
        let c = Correction::for_target_pf(4000.0, 0.8, 1.0).unwrap();
        assert!(c.reactive_after_var.abs() < EPS);
        assert!((c.capacitor_var - c.reactive_before_var).abs() < EPS);
    }

    #[test]
    fn for_triangle_agrees_with_scalar_form() {
        let load = PowerTriangle::from_vi_pf(400.0, 50.0, 0.75, Phase::Lagging).unwrap();
        let a = Correction::for_triangle(&load, 0.98).unwrap();
        let b = Correction::for_target_pf(load.real_w, 0.75, 0.98).unwrap();
        assert!((a.capacitor_var - b.capacitor_var).abs() < 1e-6);
        assert!((a.real_w - b.real_w).abs() < EPS);
    }

    #[test]
    fn rejects_no_improvement() {
        // Equal target -> nothing to do.
        assert!(matches!(
            Correction::for_target_pf(1000.0, 0.9, 0.9),
            Err(PowerError::NoCorrectionNeeded { .. })
        ));
        // Worse target -> rejected.
        assert!(matches!(
            Correction::for_target_pf(1000.0, 0.9, 0.8),
            Err(PowerError::NoCorrectionNeeded { .. })
        ));
    }

    #[test]
    fn rejects_zero_power_factor() {
        // tan diverges at PF = 0 on either side, so a present or target
        // power factor of exactly zero is rejected as out of range
        // before any improvement comparison.
        assert!(matches!(
            Correction::for_target_pf(1000.0, 0.0, 0.9),
            Err(PowerError::PowerFactorOutOfRange { .. })
        ));
        assert!(matches!(
            Correction::for_target_pf(1000.0, 0.5, 0.0),
            Err(PowerError::PowerFactorOutOfRange { .. })
        ));
    }

    #[test]
    fn rejects_non_positive_real_power() {
        assert!(matches!(
            Correction::for_target_pf(0.0, 0.7, 0.9),
            Err(PowerError::NonPositive { name: "real_w", .. })
        ));
    }

    #[test]
    fn for_triangle_rejects_leading_load() {
        let leading = PowerTriangle::from_vi_pf(230.0, 5.0, 0.8, Phase::Leading).unwrap();
        assert!(matches!(
            Correction::for_triangle(&leading, 0.95),
            Err(PowerError::NoCorrectionNeeded { .. })
        ));
    }

    #[test]
    fn for_triangle_rejects_unity_load() {
        let unity = PowerTriangle::from_p_q(1000.0, 0.0).unwrap();
        assert!(matches!(
            Correction::for_triangle(&unity, 0.99),
            Err(PowerError::NoCorrectionNeeded { .. })
        ));
    }
}
