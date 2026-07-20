//! Power delivered to a resistive load by a Thevenin source, and the
//! maximum power transfer theorem.
//!
//! ## Model
//!
//! A Thevenin source (`v_th`, `r_th`) drives a resistive load `r_load`.
//! The series current is `i = v_th / (r_th + r_load)`, so the power
//! dissipated in the load is
//!
//! ```text
//! P(r_load) = i^2 * r_load = v_th^2 * r_load / (r_th + r_load)^2
//! ```
//!
//! Differentiating with respect to `r_load` and setting `dP/dr_load = 0`
//! gives the stationary point at `r_load = r_th` (the maximum power
//! transfer theorem). Substituting back:
//!
//! ```text
//! P_max = v_th^2 / (4 * r_th)        at  r_load = r_th
//! ```
//!
//! At that matched load the efficiency is exactly 50%: half the delivered
//! power is dissipated internally in `r_th`. These closed-form facts are
//! the analytic ground truths checked by the tests, including a numeric
//! sweep confirming no load resistance beats `P_max`.

use crate::equivalent::Thevenin;
use crate::error::{non_negative, TheveninError};

/// Power `P` (watts) dissipated in a resistive load `r_load` (>= 0 ohms)
/// driven by Thevenin source `src`.
///
/// Uses `P = v_th^2 * r_load / (r_th + r_load)^2`.
///
/// # Errors
///
/// Returns [`TheveninError::Negative`] / [`TheveninError::NonFinite`] if
/// `r_load` is negative or non-finite.
pub fn load_power(src: Thevenin, r_load: f64) -> Result<f64, TheveninError> {
    let r_load = non_negative("r_load", r_load)?;
    let denom = src.r_th + r_load;
    // `r_th > 0` and `r_load >= 0`, so `denom > 0`; no division-by-zero.
    Ok(src.v_th * src.v_th * r_load / (denom * denom))
}

/// The maximum power `P_max = v_th^2 / (4 r_th)` (watts) deliverable to a
/// matched resistive load.
///
/// This is the value of [`load_power`] at `r_load = r_th`.
pub fn max_power(src: Thevenin) -> f64 {
    src.v_th * src.v_th / (4.0 * src.r_th)
}

/// The load resistance that maximises delivered power, namely
/// `r_load = r_th` (the maximum power transfer theorem).
pub fn matched_load(src: Thevenin) -> f64 {
    src.r_th
}

/// Power-transfer efficiency `eta = P_load / P_total` for a given load,
/// where `P_total` is the power supplied by the ideal source. With a
/// purely resistive divider this reduces to `eta = r_load / (r_th + r_load)`.
///
/// Returns a dimensionless ratio in `[0, 1)`; it equals exactly `0.5` at
/// the matched load `r_load = r_th`.
///
/// # Errors
///
/// Returns an error if `r_load` is negative or non-finite.
pub fn efficiency(src: Thevenin, r_load: f64) -> Result<f64, TheveninError> {
    let r_load = non_negative("r_load", r_load)?;
    Ok(r_load / (src.r_th + r_load))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    #[test]
    fn max_power_textbook_value() {
        // V_th = 12 V, R_th = 6 ohm  =>  P_max = 144 / 24 = 6 W.
        let src = Thevenin::new(12.0, 6.0).unwrap();
        assert!(close(max_power(src), 6.0));
        assert!(close(matched_load(src), 6.0));
    }

    #[test]
    fn load_power_at_matched_load_equals_max_power() {
        // Ground truth: P(r_th) == P_max for several sources.
        for &(v, r) in &[(12.0, 6.0), (5.0, 2.0), (1.0, 1.0), (100.0, 50.0)] {
            let src = Thevenin::new(v, r).unwrap();
            let p_at_match = load_power(src, r).unwrap();
            assert!(
                close(p_at_match, max_power(src)),
                "P(r_th)={p_at_match} != P_max={} for v={v}, r={r}",
                max_power(src)
            );
        }
    }

    #[test]
    fn load_power_formula_spot_check() {
        // V_th = 10, R_th = 5, R_load = 15:
        // i = 10/20 = 0.5 A ; P = 0.5^2 * 15 = 3.75 W.
        let src = Thevenin::new(10.0, 5.0).unwrap();
        assert!(close(load_power(src, 15.0).unwrap(), 3.75));
    }

    #[test]
    fn load_power_is_zero_at_open_and_short() {
        let src = Thevenin::new(8.0, 4.0).unwrap();
        // R_load = 0 (short): all volts dropped internally, load power 0.
        assert!(close(load_power(src, 0.0).unwrap(), 0.0));
        // Very large load (approaching open circuit): power -> 0.
        assert!(load_power(src, 1.0e15).unwrap() < 1.0e-9);
    }

    #[test]
    fn matched_load_maximises_power_sweep() {
        // Independent confirmation of the theorem: sweep R_load over a wide
        // range and assert NOTHING beats P(r_th). Step deliberately avoids
        // landing exactly on r_th so the comparison is strict elsewhere.
        let src = Thevenin::new(12.0, 6.0).unwrap();
        let p_match = max_power(src);
        let mut r = 0.01_f64;
        while r < 600.0 {
            let p = load_power(src, r).unwrap();
            // Allow a tiny tolerance so the near-match samples don't trip.
            assert!(
                p <= p_match + 1.0e-9,
                "found r_load={r} with power {p} exceeding P_max={p_match}"
            );
            r += 0.37;
        }
    }

    #[test]
    fn power_is_symmetric_about_match_in_log_resistance() {
        // Classic property: P(k*r_th) == P(r_th/k). Check k = 3.
        let src = Thevenin::new(10.0, 5.0).unwrap();
        let r = src.r_th;
        let hi = load_power(src, 3.0 * r).unwrap();
        let lo = load_power(src, r / 3.0).unwrap();
        assert!(close(hi, lo));
    }

    #[test]
    fn efficiency_is_half_at_match() {
        // Ground truth: efficiency is exactly 50% when R_load = R_th.
        let src = Thevenin::new(20.0, 10.0).unwrap();
        assert!(close(efficiency(src, 10.0).unwrap(), 0.5));
        // Light load (R_load >> R_th) -> efficiency -> 1.
        assert!((efficiency(src, 1.0e9).unwrap() - 1.0).abs() < 1.0e-6);
        // Short (R_load = 0) -> efficiency 0.
        assert!(close(efficiency(src, 0.0).unwrap(), 0.0));
    }

    #[test]
    fn power_and_efficiency_reject_negative_load() {
        let src = Thevenin::new(5.0, 2.0).unwrap();
        assert!(matches!(
            load_power(src, -1.0),
            Err(TheveninError::Negative { .. })
        ));
        assert!(matches!(
            efficiency(src, -0.5),
            Err(TheveninError::Negative { .. })
        ));
        assert!(matches!(
            load_power(src, f64::NAN),
            Err(TheveninError::NonFinite { .. })
        ));
    }
}
