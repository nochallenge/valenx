//! Constrained (restrained) thermal stress.
//!
//! When a uniform bar is heated by `dT` it *wants* to grow by the free
//! thermal strain `eps_thermal = alpha * dT`. If both ends are rigidly
//! held so the bar cannot change length, an equal and opposite mechanical
//! strain `-eps_thermal` is forced on it, and Hooke's law turns that into
//! a stress
//!
//! ```text
//! sigma = E * alpha * dT
//! ```
//!
//! where `E` is Young's modulus. By the usual sign convention a positive
//! `dT` (heating a fully restrained bar) produces a *compressive* stress;
//! this function returns the signed value `E * alpha * dT`, so heating
//! gives a positive number whose physical interpretation is compression
//! and cooling gives a negative number interpreted as tension. The
//! magnitude is what design checks compare against an allowable.
//!
//! `E` is in pascals (Pa), `alpha` in inverse kelvin (1/K) and `dT` in
//! kelvin, so `sigma` comes out in pascals.

use crate::error::{require_finite, require_positive, ThermalError};
use crate::expansion::LinearCoefficient;
use serde::{Deserialize, Serialize};

/// Young's modulus (modulus of elasticity), in pascals (Pa).
///
/// A newtype so a modulus is not confused with a stress or a pressure.
/// Construct with [`YoungsModulus::new`], which rejects non-positive and
/// non-finite values.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct YoungsModulus(f64);

impl YoungsModulus {
    /// Build a validated Young's modulus.
    ///
    /// # Errors
    ///
    /// Returns [`ThermalError::NonPositive`] if `modulus <= 0` and
    /// [`ThermalError::NonFinite`] if `modulus` is `NaN` or infinite.
    pub fn new(modulus: f64) -> Result<Self, ThermalError> {
        Ok(Self(require_positive("youngs_modulus", modulus)?))
    }

    /// Build a Young's modulus from a value in gigapascals (GPa).
    ///
    /// Convenience for the usual way moduli are quoted; `200.0` GPa for
    /// structural steel becomes `200e9` Pa.
    ///
    /// # Errors
    ///
    /// Same conditions as [`YoungsModulus::new`].
    pub fn from_gpa(gpa: f64) -> Result<Self, ThermalError> {
        Self::new(gpa * 1.0e9)
    }

    /// The modulus value in pascals.
    pub fn pascals(self) -> f64 {
        self.0
    }
}

/// The fully-constrained thermal stress `sigma = E * alpha * dT`, in
/// pascals.
///
/// The returned value is signed: with the convention used here a positive
/// result corresponds to compression (heating a restrained bar) and a
/// negative result to tension (cooling). A zero `delta_t` gives exactly
/// `0.0`.
///
/// This is the *fully* constrained case (restraint factor 1). Real joints
/// are rarely perfectly rigid; scale the result by a restraint factor in
/// `[0, 1]` with [`constrained_thermal_stress_restrained`] if you have one.
///
/// # Errors
///
/// Returns [`ThermalError::NonFinite`] if `delta_t` is not finite. `alpha`
/// and `youngs_modulus` are already validated by their newtypes.
pub fn constrained_thermal_stress(
    youngs_modulus: YoungsModulus,
    alpha: LinearCoefficient,
    delta_t: f64,
) -> Result<f64, ThermalError> {
    let delta_t = require_finite("delta_t", delta_t)?;
    Ok(youngs_modulus.pascals() * alpha.per_kelvin() * delta_t)
}

/// The constrained thermal stress scaled by a partial restraint factor:
/// `sigma = restraint * E * alpha * dT`.
///
/// A `restraint` of `1.0` reproduces [`constrained_thermal_stress`] (both
/// ends rigidly fixed); `0.0` means the part is free to expand and the
/// stress is zero; intermediate values model a finite support stiffness.
///
/// # Errors
///
/// Returns [`ThermalError::NonFinite`] if `delta_t` or `restraint` is not
/// finite, and [`ThermalError::NonPositive`] if `restraint < 0` (a
/// negative restraint factor has no physical meaning). `restraint == 0.0`
/// is permitted and yields zero stress.
pub fn constrained_thermal_stress_restrained(
    youngs_modulus: YoungsModulus,
    alpha: LinearCoefficient,
    delta_t: f64,
    restraint: f64,
) -> Result<f64, ThermalError> {
    let restraint = require_finite("restraint", restraint)?;
    if restraint < 0.0 {
        return Err(ThermalError::NonPositive {
            name: "restraint",
            value: restraint,
        });
    }
    Ok(restraint * constrained_thermal_stress(youngs_modulus, alpha, delta_t)?)
}

/// The free thermal strain `eps = alpha * dT` (dimensionless).
///
/// This is the strain a part would experience if it were unconstrained;
/// the constrained stress is just `E` times this strain in magnitude.
///
/// # Errors
///
/// Returns [`ThermalError::NonFinite`] if `delta_t` is not finite.
pub fn free_thermal_strain(alpha: LinearCoefficient, delta_t: f64) -> Result<f64, ThermalError> {
    let delta_t = require_finite("delta_t", delta_t)?;
    Ok(alpha.per_kelvin() * delta_t)
}

/// The temperature rise that first closes a free-expansion `gap` (metres)
/// for a bar of length `length` (metres) and linear coefficient `alpha`,
/// inverting the free expansion `dL = alpha * L0 * dT`:
///
/// ```text
/// dT_gap = gap / (alpha * L0)
/// ```
///
/// Below this rise the bar grows freely into the gap and carries no
/// stress; at and above it the gap is closed and further heating is
/// resisted (see [`gap_constrained_thermal_stress`]). A zero gap closes
/// at `dT_gap = 0`.
///
/// # Errors
///
/// Returns [`ThermalError::NonPositive`] if `length` is not strictly
/// positive or if `gap` is negative, and [`ThermalError::NonFinite`] if
/// either is not finite. `alpha` is validated by its newtype.
pub fn gap_closure_temperature(
    alpha: LinearCoefficient,
    length: f64,
    gap: f64,
) -> Result<f64, ThermalError> {
    let length = require_positive("length", length)?;
    let gap = require_finite("gap", gap)?;
    if gap < 0.0 {
        return Err(ThermalError::NonPositive {
            name: "gap",
            value: gap,
        });
    }
    Ok(gap / (alpha.per_kelvin() * length))
}

/// The thermal stress in a bar that expands freely until it closes a
/// `gap` (metres) to a rigid wall, then is constrained — the classic
/// "bar with a gap" problem:
///
/// ```text
/// sigma = E * max(alpha * dT - gap / L0, 0)
/// ```
///
/// The bar of length `length` first takes up the gap (its free expansion
/// `alpha * L0 * dT` growing to `gap`); only the *excess* strain
/// `alpha * dT - gap / L0` beyond that is suppressed by the wall, and
/// Hooke's law turns it into a compressive stress. The result is
/// therefore non-negative: it is exactly zero while the gap is still open
/// (`alpha * dT <= gap / L0`, including the
/// [closure temperature](gap_closure_temperature)) and, because a gap can
/// only resist *growth* into the wall, also zero under any cooling
/// (`dT < 0`). This is the one difference from
/// [`constrained_thermal_stress`], which (modelling a bar held at both
/// ends) goes into tension on cooling.
///
/// Two exact identities tie it back to the fully-constrained case: a zero
/// gap reproduces `constrained_thermal_stress` for heating, and once the
/// gap is closed the stress equals the fully-constrained stress of the
/// *excess* temperature, `sigma(dT) = E * alpha * (dT - dT_gap)`.
///
/// # Errors
///
/// Returns [`ThermalError::NonPositive`] if `length` is not strictly
/// positive or if `gap` is negative, and [`ThermalError::NonFinite`] if
/// `delta_t`, `length`, or `gap` is not finite. `alpha` and
/// `youngs_modulus` are validated by their newtypes.
pub fn gap_constrained_thermal_stress(
    youngs_modulus: YoungsModulus,
    alpha: LinearCoefficient,
    length: f64,
    delta_t: f64,
    gap: f64,
) -> Result<f64, ThermalError> {
    let length = require_positive("length", length)?;
    let delta_t = require_finite("delta_t", delta_t)?;
    let gap = require_finite("gap", gap)?;
    if gap < 0.0 {
        return Err(ThermalError::NonPositive {
            name: "gap",
            value: gap,
        });
    }
    // Free thermal strain minus the gap expressed as a strain; the wall
    // only resists growth into it, so a still-open gap or any cooling
    // gives zero.
    let excess_strain = (alpha.per_kelvin() * delta_t - gap / length).max(0.0);
    Ok(youngs_modulus.pascals() * excess_strain)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expansion::linear_expansion;

    const EPS_STRAIN: f64 = 1e-12;
    /// Stresses are O(1e6)-O(1e8) Pa here, so an absolute tolerance scaled
    /// to that magnitude (1e-3 Pa, i.e. micro-Pa relative) is appropriate.
    const EPS_STRESS: f64 = 1e-3;

    fn steel_e() -> YoungsModulus {
        // Structural steel, 200 GPa.
        YoungsModulus::from_gpa(200.0).unwrap()
    }

    fn steel_alpha() -> LinearCoefficient {
        // Structural steel CTE, 12e-6 / K.
        LinearCoefficient::new(12.0e-6).unwrap()
    }

    #[test]
    fn stress_matches_e_alpha_dt() {
        let e = steel_e();
        let alpha = steel_alpha();
        let dt = 50.0;
        let got = constrained_thermal_stress(e, alpha, dt).unwrap();
        let expected = 200.0e9 * 12.0e-6 * dt;
        assert!(
            (got - expected).abs() < EPS_STRESS,
            "sigma: got {got}, expected {expected}"
        );
        // Hand-checked: 200e9 * 12e-6 * 50 = 1.2e8 Pa = 120 MPa.
        assert!(
            (got - 120.0e6).abs() < EPS_STRESS,
            "sigma magnitude wrong: {got}"
        );
    }

    #[test]
    fn stress_is_e_times_free_strain() {
        let e = steel_e();
        let alpha = steel_alpha();
        let dt = 73.0;
        let strain = free_thermal_strain(alpha, dt).unwrap();
        let sigma = constrained_thermal_stress(e, alpha, dt).unwrap();
        assert!(
            (sigma - e.pascals() * strain).abs() < EPS_STRESS,
            "sigma != E * strain"
        );
        // Free strain itself is alpha * dT.
        assert!(
            (strain - 12.0e-6 * dt).abs() < EPS_STRAIN,
            "strain != alpha dT"
        );
    }

    #[test]
    fn heating_and_cooling_are_opposite_signs() {
        let e = steel_e();
        let alpha = steel_alpha();
        let hot = constrained_thermal_stress(e, alpha, 60.0).unwrap();
        let cold = constrained_thermal_stress(e, alpha, -60.0).unwrap();
        assert!(
            hot > 0.0,
            "heating restrained bar should be positive: {hot}"
        );
        assert!(
            cold < 0.0,
            "cooling restrained bar should be negative: {cold}"
        );
        assert!((hot + cold).abs() < EPS_STRESS, "not antisymmetric");
    }

    #[test]
    fn zero_delta_t_gives_zero_stress() {
        let e = steel_e();
        let alpha = steel_alpha();
        assert!(constrained_thermal_stress(e, alpha, 0.0).unwrap().abs() < EPS_STRESS);
        assert!(free_thermal_strain(alpha, 0.0).unwrap().abs() < EPS_STRAIN);
    }

    #[test]
    fn restraint_factor_scales_linearly() {
        let e = steel_e();
        let alpha = steel_alpha();
        let dt = 50.0;
        let full = constrained_thermal_stress(e, alpha, dt).unwrap();
        // restraint = 1 reproduces the fully constrained case.
        let r1 = constrained_thermal_stress_restrained(e, alpha, dt, 1.0).unwrap();
        assert!((r1 - full).abs() < EPS_STRESS, "restraint 1 != full");
        // restraint = 0 means free expansion, zero stress.
        let r0 = constrained_thermal_stress_restrained(e, alpha, dt, 0.0).unwrap();
        assert!(r0.abs() < EPS_STRESS, "restraint 0 should be zero: {r0}");
        // restraint = 0.5 is exactly half.
        let rhalf = constrained_thermal_stress_restrained(e, alpha, dt, 0.5).unwrap();
        assert!(
            (rhalf - 0.5 * full).abs() < EPS_STRESS,
            "restraint 0.5 != half"
        );
    }

    #[test]
    fn restraint_rejects_negative_and_nonfinite() {
        let e = steel_e();
        let alpha = steel_alpha();
        assert!(matches!(
            constrained_thermal_stress_restrained(e, alpha, 10.0, -0.1),
            Err(ThermalError::NonPositive {
                name: "restraint",
                ..
            })
        ));
        assert!(matches!(
            constrained_thermal_stress_restrained(e, alpha, 10.0, f64::INFINITY),
            Err(ThermalError::NonFinite {
                name: "restraint",
                ..
            })
        ));
    }

    #[test]
    fn stress_rejects_nonfinite_delta_t() {
        let e = steel_e();
        let alpha = steel_alpha();
        assert!(matches!(
            constrained_thermal_stress(e, alpha, f64::NAN),
            Err(ThermalError::NonFinite {
                name: "delta_t",
                ..
            })
        ));
    }

    #[test]
    fn modulus_constructors_validate() {
        assert!(matches!(
            YoungsModulus::new(0.0),
            Err(ThermalError::NonPositive {
                name: "youngs_modulus",
                ..
            })
        ));
        assert!(matches!(
            YoungsModulus::from_gpa(-5.0),
            Err(ThermalError::NonPositive { .. })
        ));
        // from_gpa scales by 1e9.
        assert!((YoungsModulus::from_gpa(70.0).unwrap().pascals() - 70.0e9).abs() < 1.0);
    }

    // --- bar with a gap ---------------------------------------------------

    #[test]
    fn gap_closure_temperature_matches_inverse_of_free_expansion() {
        // alpha = 12e-6 /K, L0 = 2 m, gap = 1 mm -> dT_gap = 41.6667 K.
        let alpha = steel_alpha();
        let dt_gap = gap_closure_temperature(alpha, 2.0, 1.0e-3).unwrap();
        let expected = 1.0e-3 / (12.0e-6 * 2.0);
        assert!(
            (dt_gap - expected).abs() < 1e-9,
            "dT_gap {dt_gap} vs {expected}"
        );
        assert!(
            (dt_gap - 41.666_666_666_666_664).abs() < 1e-9,
            "dT_gap {dt_gap}"
        );
        // At exactly dT_gap the free expansion equals the gap.
        let dl = linear_expansion(alpha, 2.0, dt_gap).unwrap();
        assert!((dl - 1.0e-3).abs() < 1e-15, "free expansion closes the gap");
        // A zero gap closes immediately.
        assert!(gap_closure_temperature(alpha, 2.0, 0.0).unwrap().abs() < 1e-12);
    }

    #[test]
    fn zero_gap_reproduces_full_constraint_when_heating() {
        // GOLD identity: gap = 0 -> gap stress equals the fully-constrained
        // stress for any heating dT.
        let e = steel_e();
        let alpha = steel_alpha();
        for &dt in &[5.0, 40.0, 120.0] {
            let gap = gap_constrained_thermal_stress(e, alpha, 1.5, dt, 0.0).unwrap();
            let full = constrained_thermal_stress(e, alpha, dt).unwrap();
            assert!((gap - full).abs() < EPS_STRESS, "g=0 != full at dT={dt}");
        }
    }

    #[test]
    fn closed_gap_stress_equals_constraint_of_excess_temperature() {
        // Once the gap is closed the stress equals the fully-constrained
        // stress of the temperature beyond closure: sigma = E alpha (dT - dT_gap).
        let e = steel_e();
        let alpha = steel_alpha();
        let (l0, gap, dt) = (2.0, 1.0e-3, 80.0);
        let dt_gap = gap_closure_temperature(alpha, l0, gap).unwrap();
        let sigma = gap_constrained_thermal_stress(e, alpha, l0, dt, gap).unwrap();
        let via_excess = constrained_thermal_stress(e, alpha, dt - dt_gap).unwrap();
        assert!((sigma - via_excess).abs() < EPS_STRESS, "excess identity");
        // Independent closed form: E (alpha dT - gap/L0)
        //   = 200e9 (12e-6*80 - 1e-3/2) = 200e9 * 4.6e-4 = 92 MPa.
        assert!(
            (sigma - 92.0e6).abs() < EPS_STRESS,
            "sigma = {sigma}, want 92 MPa"
        );
    }

    #[test]
    fn no_stress_until_the_gap_closes() {
        let e = steel_e();
        let alpha = steel_alpha();
        let (l0, gap) = (2.0, 1.0e-3);
        let dt_gap = gap_closure_temperature(alpha, l0, gap).unwrap();
        // Below closure: free expansion still inside the gap -> zero stress.
        let below = gap_constrained_thermal_stress(e, alpha, l0, 0.9 * dt_gap, gap).unwrap();
        assert!(below.abs() < EPS_STRESS, "stress below closure: {below}");
        // Exactly at closure: just touching, still zero.
        let at = gap_constrained_thermal_stress(e, alpha, l0, dt_gap, gap).unwrap();
        assert!(at.abs() < EPS_STRESS, "stress at closure: {at}");
        // Just above closure: positive (compressive).
        let above = gap_constrained_thermal_stress(e, alpha, l0, 1.01 * dt_gap, gap).unwrap();
        assert!(above > 0.0, "stress just above closure: {above}");
    }

    #[test]
    fn cooling_a_gapped_bar_gives_zero_stress() {
        // A one-sided gap resists only growth, so cooling never stresses
        // the bar -- even with a zero gap (unlike the both-ends-fixed
        // constrained case, which would go into tension).
        let e = steel_e();
        let alpha = steel_alpha();
        let cold = gap_constrained_thermal_stress(e, alpha, 2.0, -60.0, 1.0e-3).unwrap();
        assert!(cold.abs() < EPS_STRESS, "cooling with gap: {cold}");
        let cold0 = gap_constrained_thermal_stress(e, alpha, 2.0, -60.0, 0.0).unwrap();
        assert!(cold0.abs() < EPS_STRESS, "cooling with zero gap: {cold0}");
        // The both-ends-fixed model, by contrast, is in tension here.
        assert!(constrained_thermal_stress(e, alpha, -60.0).unwrap() < 0.0);
    }

    #[test]
    fn stress_rises_with_heating_and_falls_with_a_wider_gap() {
        let e = steel_e();
        let alpha = steel_alpha();
        let (l0, gap) = (2.0, 1.0e-3);
        // Monotonic increase past closure.
        let s60 = gap_constrained_thermal_stress(e, alpha, l0, 60.0, gap).unwrap();
        let s90 = gap_constrained_thermal_stress(e, alpha, l0, 90.0, gap).unwrap();
        assert!(s90 > s60 && s60 > 0.0, "stress grows with heating");
        // A wider gap closes later, so less stress at the same dT.
        let wide = gap_constrained_thermal_stress(e, alpha, l0, 90.0, 2.0e-3).unwrap();
        assert!(wide < s90, "wider gap -> lower stress: {wide} < {s90}");
    }

    #[test]
    fn gap_stress_rejects_bad_inputs() {
        let e = steel_e();
        let alpha = steel_alpha();
        assert!(matches!(
            gap_constrained_thermal_stress(e, alpha, 0.0, 50.0, 1.0e-3),
            Err(ThermalError::NonPositive { name: "length", .. })
        ));
        assert!(matches!(
            gap_constrained_thermal_stress(e, alpha, 2.0, 50.0, -1.0e-3),
            Err(ThermalError::NonPositive { name: "gap", .. })
        ));
        assert!(matches!(
            gap_constrained_thermal_stress(e, alpha, 2.0, f64::NAN, 1.0e-3),
            Err(ThermalError::NonFinite {
                name: "delta_t",
                ..
            })
        ));
        assert!(matches!(
            gap_closure_temperature(alpha, -2.0, 1.0e-3),
            Err(ThermalError::NonPositive { name: "length", .. })
        ));
        assert!(matches!(
            gap_closure_temperature(alpha, 2.0, f64::INFINITY),
            Err(ThermalError::NonFinite { name: "gap", .. })
        ));
    }
}
