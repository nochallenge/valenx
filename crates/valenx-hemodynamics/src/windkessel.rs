//! Two-element (RC) Windkessel diastolic pressure decay.
//!
//! The classic 2-element Windkessel lumps the arterial tree into a
//! single resistance `R` (the peripheral / systemic vascular
//! resistance) in parallel with a single compliance `C` (the elastic
//! storage of the large arteries). During diastole the aortic valve is
//! shut, so no new flow enters; the compliance discharges its stored
//! volume through the resistance and the pressure relaxes exponentially.
//!
//! With inflow set to zero the governing balance `C dP/dt + P/R = 0`
//! integrates to:
//!
//! ```text
//! P(t) = P0 * exp(-t / (R * C))
//! ```
//!
//! where `P0` is the pressure at the start of diastole and `tau = R*C`
//! is the time constant of the decay. After one time constant the
//! pressure has fallen to `P0 / e` (about 37 % of `P0`).
//!
//! # Units
//!
//! The exponent must be dimensionless, so `R*C` and the elapsed time
//! `t` must share a time unit. In SI, `R` is Pa·s/m^3 and `C` is
//! m^3/Pa, so `R*C` is in seconds; `P0` may be in any pressure unit and
//! the result comes back in that same unit.

use crate::error::{require_non_negative, require_positive};
use crate::HemodynamicsError;

/// Diastolic arterial pressure at time `t` under the 2-element
/// Windkessel decay `P(t) = P0 * exp(-t / (R * C))`.
///
/// Models the exponential relaxation of arterial pressure during
/// diastole, when the elastic arteries (compliance `c`) discharge
/// through the peripheral resistance (`r`) with no new inflow.
///
/// # Units
///
/// `r` and `c` must combine to a time in the same unit as `t` (SI:
/// `r` in Pa·s/m^3, `c` in m^3/Pa, `t` in s). The result is in the same
/// pressure unit as `p0`.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] / [`HemodynamicsError::NotFinite`]
/// if `r` or `c` is not strictly positive, or
/// [`HemodynamicsError::Negative`] / [`HemodynamicsError::NotFinite`]
/// if `t` is negative; `p0` is only required to be finite.
pub fn windkessel_pressure(p0: f64, t: f64, r: f64, c: f64) -> Result<f64, HemodynamicsError> {
    if !p0.is_finite() {
        return Err(HemodynamicsError::NotFinite {
            name: "p0",
            value: p0,
        });
    }
    let t = require_non_negative("t", t)?;
    let r = require_positive("r", r)?;
    let c = require_positive("c", c)?;
    let tau = r * c;
    Ok(p0 * (-t / tau).exp())
}

/// The Windkessel time constant `tau = R * C`.
///
/// After this much elapsed time the diastolic pressure has decayed to
/// `P0 / e`. Exposed as its own function because the time constant is
/// the single parameter that characterises the whole decay.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] / [`HemodynamicsError::NotFinite`]
/// if `r` or `c` is not a strictly-positive finite number.
pub fn time_constant(r: f64, c: f64) -> Result<f64, HemodynamicsError> {
    let r = require_positive("r", r)?;
    let c = require_positive("c", c)?;
    Ok(r * c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;

    const P0: f64 = 80.0; // e.g. 80 mmHg at end-systole
    const R: f64 = 1.2e8; // Pa·s/m^3
    const C: f64 = 1.0e-8; // m^3/Pa  -> tau = R*C = 1.2 s

    #[test]
    fn pressure_at_zero_time_is_p0() {
        let p = windkessel_pressure(P0, 0.0, R, C).expect("valid");
        assert!((p - P0).abs() < 1e-12);
    }

    #[test]
    fn decays_to_p0_over_e_at_one_time_constant() {
        // VALIDATE: at t = R*C the pressure is P0 / e.
        let tau = time_constant(R, C).expect("valid");
        let p = windkessel_pressure(P0, tau, R, C).expect("valid");
        let expected = P0 / E;
        assert!(
            (p - expected).abs() < 1e-9,
            "at t=tau expected {expected}, got {p}"
        );
        // And the decay factor itself is 1/e ~ 0.3679.
        assert!((p / P0 - 1.0 / E).abs() < 1e-12);
    }

    #[test]
    fn decays_to_p0_over_e_squared_at_two_time_constants() {
        let tau = time_constant(R, C).expect("valid");
        let p = windkessel_pressure(P0, 2.0 * tau, R, C).expect("valid");
        let expected = P0 / (E * E);
        assert!((p - expected).abs() < 1e-9, "expected {expected}, got {p}");
    }

    #[test]
    fn pressure_is_monotonically_decreasing() {
        let mut prev = windkessel_pressure(P0, 0.0, R, C).expect("valid");
        for k in 1..=20 {
            let t = 0.1 * f64::from(k);
            let p = windkessel_pressure(P0, t, R, C).expect("valid");
            assert!(p < prev, "not decreasing at t={t}: {p} >= {prev}");
            assert!(p > 0.0);
            prev = p;
        }
    }

    #[test]
    fn time_constant_is_r_times_c() {
        let tau = time_constant(R, C).expect("valid");
        assert!((tau - R * C).abs() < 1e-9 * (R * C));
        assert!((tau - 1.2).abs() < 1e-6);
    }

    #[test]
    fn matches_closed_form_at_arbitrary_time() {
        let t = 0.37;
        let p = windkessel_pressure(P0, t, R, C).expect("valid");
        let expected = P0 * (-t / (R * C)).exp();
        assert!((p - expected).abs() < 1e-12 * expected);
    }

    #[test]
    fn larger_time_constant_decays_slower() {
        // Doubling compliance doubles tau; at the same elapsed time the
        // pressure is higher (slower decay).
        let t = 1.2;
        let p_fast = windkessel_pressure(P0, t, R, C).expect("valid");
        let p_slow = windkessel_pressure(P0, t, R, 2.0 * C).expect("valid");
        assert!(p_slow > p_fast);
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(windkessel_pressure(P0, -1.0, R, C).is_err());
        assert!(windkessel_pressure(P0, 1.0, 0.0, C).is_err());
        assert!(windkessel_pressure(P0, 1.0, R, -1.0).is_err());
        assert!(windkessel_pressure(f64::NAN, 1.0, R, C).is_err());
        assert!(time_constant(0.0, C).is_err());
    }
}
