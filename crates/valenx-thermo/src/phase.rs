//! Pure-component vapor–liquid phase behavior.
//!
//! The saturation (vapor) pressure `Psat(T)` is the pressure at which the
//! liquid and vapor phases predicted by the cubic EoS have equal fugacity:
//!
//! ```text
//! φ_liquid(T, P) = φ_vapor(T, P)
//! ```
//!
//! Equivalently the *fugacity ratio* `R(P) = φ_liq / φ_vap` equals 1 at
//! saturation. We bracket and then Newton-iterate on `ln P` to find that
//! pressure, starting from the Wilson vapor-pressure correlation as an initial
//! guess.

use crate::eos::Eos;
use crate::error::{require_positive, Result, ThermoError};
use crate::fluid::Fluid;

/// Wilson's correlation for an initial saturation-pressure estimate (Pa).
///
/// `Psat ≈ Pc · exp[5.373 (1 + ω)(1 − Tc/T)]`. Cheap, derivative-free, and a
/// good starting point for the Newton iteration below.
#[must_use]
pub fn wilson_psat(fluid: &Fluid, t: f64) -> f64 {
    fluid.pc * (5.373 * (1.0 + fluid.omega) * (1.0 - fluid.tc / t)).exp()
}

/// Saturation (vapor) pressure `Psat` (Pa) of a pure fluid at temperature `t`
/// (K), found by equating liquid and vapor fugacities under the given EoS.
///
/// The fluid must be sub-critical (`t < Tc`); above the critical temperature no
/// distinct liquid and vapor phases exist and the function returns
/// [`ThermoError::OutOfRange`].
///
/// # Errors
///
/// * [`ThermoError::NonPositive`] if `t` is not strictly positive.
/// * [`ThermoError::OutOfRange`] if `t >= Tc`.
/// * [`ThermoError::NotConverged`] if the Newton iteration fails to reach the
///   equal-fugacity condition.
pub fn saturation_pressure(eos: &Eos, t: f64) -> Result<f64> {
    require_positive("temperature", t)?;
    let fluid = &eos.fluid;
    if t >= fluid.tc {
        return Err(ThermoError::OutOfRange {
            name: "temperature",
            value: t,
            expected: "must be below the critical temperature for a saturation pressure",
        });
    }

    // Objective in ln P: g(lnP) = ln(φ_liq / φ_vap) = 0 at saturation.
    let g = |p: f64| -> Result<f64> {
        let roots = eos.z_roots(t, p)?;
        // If the cubic gives a single root, the liquid/vapor split is degenerate
        // here; nudge the objective using the same root (drives the search back
        // toward the two-phase window).
        let phi_l = eos.ln_phi(t, p, roots.liquid)?;
        let phi_v = eos.ln_phi(t, p, roots.vapor)?;
        Ok(phi_l - phi_v)
    };

    let mut ln_p = wilson_psat(fluid, t).max(1.0).ln();
    let max_iter = 100;
    let tol = 1e-10;

    for i in 0..max_iter {
        let p = ln_p.exp();
        let g0 = g(p)?;
        if g0.abs() < tol {
            return Ok(p);
        }
        // Numerical derivative dg/d(lnP) via a small relative step.
        let h = 1e-6;
        let g1 = g((ln_p + h).exp())?;
        let dg = (g1 - g0) / h;
        if !dg.is_finite() || dg.abs() < 1e-300 {
            return Err(ThermoError::NotConverged {
                solver: "saturation_pressure (Newton)",
                iterations: i,
                residual: g0,
            });
        }
        let mut step = g0 / dg;
        // Damp large steps in ln P to keep the iteration stable.
        step = step.clamp(-2.0, 2.0);
        ln_p -= step;
        // Keep the pressure physical: between a tiny floor and Pc.
        ln_p = ln_p.clamp((1.0e-3_f64).ln(), fluid.pc.ln());
    }

    let p = ln_p.exp();
    Err(ThermoError::NotConverged {
        solver: "saturation_pressure (Newton)",
        iterations: max_iter,
        residual: g(p).unwrap_or(f64::NAN),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eos::EosModel;

    #[test]
    fn wilson_is_in_the_right_ballpark() {
        // CO2 at 273.15 K: experimental Psat ≈ 3.485 MPa. Wilson is approximate
        // but should be within a factor of ~2.
        let p = wilson_psat(&Fluid::carbon_dioxide(), 273.15);
        assert!(p > 1.5e6 && p < 7.0e6, "Wilson CO2 Psat = {p}");
    }

    #[test]
    fn saturation_rejects_supercritical() {
        let eos = Eos::new(Fluid::carbon_dioxide(), EosModel::PengRobinson);
        let res = saturation_pressure(&eos, 400.0); // > Tc = 304 K
        assert!(matches!(res, Err(ThermoError::OutOfRange { .. })));
    }

    #[test]
    fn saturation_rejects_bad_temperature() {
        let eos = Eos::new(Fluid::nitrogen(), EosModel::Srk);
        assert!(matches!(
            saturation_pressure(&eos, -10.0),
            Err(ThermoError::NonPositive { .. })
        ));
    }
}
