//! Δv ↔ propellant budgeting via the rocket equation.
//!
//! Ties the [`crate::maneuver`] Δv numbers to actual propellant mass for
//! a given stage `Isp`, and evaluates a whole **mission sequence** of
//! burns to report propellant used, mass remaining, and whether the
//! stage can afford the plan. This is the bookkeeping that turns "the
//! mission needs 4.2 km/s" into "you need 11 t of propellant" — or
//! "you're 600 kg short."

use serde::{Deserialize, Serialize};

use crate::constants::G0;
use crate::error::{AstroError, Result};

/// Validate a specific impulse used as the rocket-equation divisor. A
/// non-finite or non-positive `isp` drives `Isp·g₀` to zero/NaN and would
/// otherwise leak a silent `NaN`/`Inf` through the exponential or the
/// logarithm, so it is rejected up front.
fn check_isp(isp: f64) -> Result<()> {
    if !isp.is_finite() || isp <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "specific impulse (isp) must be finite and positive",
        ));
    }
    Ok(())
}

/// Propellant mass (kg) needed to produce `delta_v` (m/s) starting from
/// `initial_mass` (kg) at the given `isp` (s), via Tsiolkovsky:
/// `m_f = m₀·e^(−Δv/(Isp·g₀))`, propellant = `m₀ − m_f`.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `isp` is non-finite or
/// non-positive — the `Isp·g₀` divisor would otherwise be zero/NaN and
/// the result a silent `NaN`/`Inf`.
pub fn propellant_for_delta_v(delta_v: f64, isp: f64, initial_mass: f64) -> Result<f64> {
    check_isp(isp)?;
    let final_mass = initial_mass * (-delta_v / (isp * G0)).exp();
    Ok(initial_mass - final_mass)
}

/// Δv (m/s) a burn delivers when it consumes `propellant` (kg) from an
/// `initial_mass` (kg) at the given `isp` (s).
///
/// A burn that empties (or over-draws) the tanks — `final_mass ≤ 0` —
/// returns `f64::INFINITY` as a **documented sentinel**: the requested
/// propellant implies an unreachable, unbounded Δv. That sentinel is
/// deliberate and preserved.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `isp` is non-finite or
/// non-positive — that is an *unintended* `NaN` source (`Isp·g₀ = 0`/NaN
/// scaling the logarithm), distinct from the `final_mass ≤ 0` Inf
/// sentinel above.
pub fn delta_v_for_propellant(propellant: f64, isp: f64, initial_mass: f64) -> Result<f64> {
    check_isp(isp)?;
    let final_mass = initial_mass - propellant;
    if final_mass <= 0.0 {
        // Documented sentinel: emptying/over-draining the tanks is an
        // infinite, unreachable Δv. Preserved exactly (not an error).
        return Ok(f64::INFINITY);
    }
    Ok(isp * G0 * (initial_mass / final_mass).ln())
}

/// The outcome of evaluating a mission Δv sequence on a stage.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BudgetResult {
    /// Total Δv requested across all burns (m/s).
    pub total_delta_v: f64,
    /// Total propellant consumed to deliver it (kg).
    pub propellant_used: f64,
    /// Stage mass after the planned burns (kg). Equals dry mass when the
    /// plan exactly empties the tanks.
    pub final_mass: f64,
    /// Whether the stage's usable propellant covers the plan.
    pub feasible: bool,
    /// Δv the stage can still deliver after the plan (m/s); zero when
    /// infeasible or empty.
    pub delta_v_margin: f64,
}

/// Evaluate a sequence of Δv burns (m/s) on a stage of `dry_mass` +
/// `propellant_mass` (kg) at the given `isp` (s).
///
/// Burns are applied in order; each consumes propellant from the running
/// mass. The plan is feasible if the burns never exhaust the usable
/// propellant.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `isp` is non-finite or
/// non-positive — the rocket equation would otherwise return a silent
/// `NaN`/`Inf` budget.
pub fn evaluate_sequence(
    dry_mass: f64,
    propellant_mass: f64,
    isp: f64,
    burns: &[f64],
) -> Result<BudgetResult> {
    check_isp(isp)?;
    let initial = dry_mass + propellant_mass;
    let mut mass = initial;
    let mut feasible = true;
    for &dv in burns {
        let used = propellant_for_delta_v(dv, isp, mass)?;
        if mass - used < dry_mass - 1e-6 {
            feasible = false;
            mass = dry_mass; // tanks dry
            break;
        }
        mass -= used;
    }
    let propellant_used = initial - mass;
    // Remaining Δv from whatever propellant is left.
    let remaining_prop = (mass - dry_mass).max(0.0);
    let delta_v_margin = if feasible {
        delta_v_for_propellant(remaining_prop, isp, mass)?
    } else {
        0.0
    };
    Ok(BudgetResult {
        total_delta_v: burns.iter().sum(),
        propellant_used,
        final_mass: mass,
        feasible,
        delta_v_margin,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn propellant_delta_v_round_trip() {
        let (isp, m0) = (340.0, 50_000.0);
        let prop = propellant_for_delta_v(3_000.0, isp, m0).expect("valid isp");
        let dv = delta_v_for_propellant(prop, isp, m0).expect("valid isp");
        assert!((dv - 3_000.0).abs() < 1e-6, "dv {dv}");
    }

    #[test]
    fn tsiolkovsky_matches_closed_form() {
        // m0=10t, isp=300 -> dv to burn 9t (mf=1t) is Isp·g0·ln(10).
        let prop = 9_000.0;
        let dv = delta_v_for_propellant(prop, 300.0, 10_000.0).expect("valid isp");
        assert!((dv - 300.0 * G0 * 10.0_f64.ln()).abs() < 1e-6);
    }

    #[test]
    fn feasible_sequence_reports_positive_margin() {
        // A stage with generous propellant covers a modest plan.
        let r = evaluate_sequence(4_000.0, 100_000.0, 348.0, &[2_400.0, 1_500.0]).expect("valid");
        assert!(r.feasible);
        assert!((r.total_delta_v - 3_900.0).abs() < 1e-9);
        assert!(r.propellant_used > 0.0 && r.final_mass > 4_000.0);
        assert!(r.delta_v_margin > 0.0);
    }

    #[test]
    fn infeasible_sequence_is_flagged() {
        // Tiny tank, huge plan -> not feasible.
        let r = evaluate_sequence(4_000.0, 5_000.0, 348.0, &[6_000.0]).expect("valid");
        assert!(!r.feasible);
        assert_eq!(r.delta_v_margin, 0.0);
    }

    // A hand-built stage with isp=0 (or non-finite) drives the rocket
    // equation's `isp·g₀` divisor to zero/NaN. Pre-fix this leaked a
    // silent NaN/Inf budget; the guard now rejects it cleanly.
    #[test]
    fn zero_or_nonfinite_isp_is_rejected_not_silent_nan() {
        // propellant_for_delta_v: dv=0, isp=0 -> 0/0 -> was NaN, now Err.
        assert!(
            matches!(
                propellant_for_delta_v(0.0, 0.0, 50_000.0),
                Err(AstroError::NonPhysicalState(_))
            ),
            "isp=0 must be rejected"
        );
        // A negative-Δv burn with isp=0 was -Inf; non-finite isp was NaN.
        assert!(propellant_for_delta_v(-100.0, 0.0, 50_000.0).is_err());
        assert!(propellant_for_delta_v(3_000.0, f64::NAN, 50_000.0).is_err());

        // evaluate_sequence propagates the rejection instead of a
        // non-finite BudgetResult.
        assert!(
            matches!(
                evaluate_sequence(4_000.0, 100_000.0, 0.0, &[0.0]),
                Err(AstroError::NonPhysicalState(_))
            ),
            "isp=0 sequence must be rejected"
        );

        // delta_v_for_propellant: isp=NaN with a valid burn was an
        // UNINTENDED NaN -> now Err.
        assert!(delta_v_for_propellant(9_000.0, f64::NAN, 10_000.0).is_err());
        assert!(delta_v_for_propellant(9_000.0, 0.0, 10_000.0).is_err());
    }

    // The `final_mass <= 0` Inf return of delta_v_for_propellant is a
    // DOCUMENTED sentinel (unreachable, unbounded Δv) and must survive the
    // isp guard untouched.
    #[test]
    fn delta_v_for_propellant_infinity_sentinel_preserved() {
        // Over-draining the tanks (propellant >= initial_mass) with a
        // valid isp still returns the documented +Inf, not an error.
        let dv = delta_v_for_propellant(12_000.0, 300.0, 10_000.0).expect("valid isp");
        assert!(dv.is_infinite() && dv > 0.0, "expected +Inf sentinel, got {dv}");
        // Exactly emptying the tanks (final_mass == 0) also hits it.
        let dv0 = delta_v_for_propellant(10_000.0, 300.0, 10_000.0).expect("valid isp");
        assert!(dv0.is_infinite() && dv0 > 0.0, "expected +Inf sentinel, got {dv0}");
    }

    #[test]
    fn exact_full_burn_lands_on_dry_mass() {
        // Burn exactly the stage's total ideal Δv -> ends ~ at dry mass.
        let (dry, prop, isp) = (2_000.0, 18_000.0, 320.0);
        let dv = delta_v_for_propellant(prop, isp, dry + prop).expect("valid isp");
        let r = evaluate_sequence(dry, prop, isp, &[dv]).expect("valid");
        assert!(r.feasible);
        assert!((r.final_mass - dry).abs() < 1.0, "final {}", r.final_mass);
        assert!(r.delta_v_margin < 1.0);
    }
}
