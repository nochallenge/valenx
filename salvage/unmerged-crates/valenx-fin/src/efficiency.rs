//! Fin efficiency and effectiveness for the insulated (adiabatic) tip.
//!
//! For a straight fin of constant cross-section with an insulated tip, the
//! temperature excess profile is `theta(x) = theta_b * cosh(m (L - x)) / cosh(mL)`
//! and the **fin efficiency** — actual heat transfer divided by the ideal
//! heat transfer if the whole fin sat at the base temperature — is
//!
//! ```text
//! eta = tanh(mL) / (mL).
//! ```
//!
//! ## Limiting cases (the analytic ground truth)
//!
//! As `mL -> 0` (a short, highly conductive fin) the series
//! `tanh(x)/x = 1 - x^2/3 + 2 x^4/15 - ...` gives `eta -> 1`: the fin is
//! isothermal and perfectly effective per unit area. As `mL` grows large,
//! `tanh(mL) -> 1` so `eta -> 1 / (mL)`: a long fin's tip is dead weight and
//! efficiency falls off like `1/(mL)`.
//!
//! The **fin effectiveness** `epsilon_f` compares the finned heat rate to the
//! bare-base heat rate; for the insulated tip it reduces to
//! `epsilon_f = sqrt(k P / (h A)) * tanh(mL)`, and adding a fin is only
//! worthwhile when `epsilon_f > 1` (a rule of thumb favours `epsilon_f >= 2`).

use crate::error::{FinError, Result};

/// Series cross-over below which `tanh(x)/x` is evaluated from its Taylor
/// expansion to avoid the `0/0` form and catastrophic cancellation near zero.
const SMALL_ML: f64 = 1.0e-4;

/// Fin efficiency `eta = tanh(mL) / (mL)` for an insulated tip.
///
/// `mlength` is the dimensionless fin length `m * L` (see
/// [`crate::dimensionless_length`]); it must be finite and non-negative.
/// At `mlength == 0` the function returns the analytic limit `1.0`, and for
/// very small `mlength` it uses the Taylor series `1 - (mL)^2/3 + ...` so the
/// result stays accurate where the naive quotient would lose precision.
///
/// The returned value always lies in `(0, 1]`.
///
/// # Errors
///
/// Returns [`FinError::Negative`] / [`FinError::NotFinite`] if `mlength` is
/// negative or non-finite.
///
/// # Example
///
/// ```
/// use valenx_fin::efficiency;
/// // mL -> 0 limit.
/// assert!((efficiency(0.0).unwrap() - 1.0).abs() < 1e-12);
/// // mL = 1: tanh(1)/1 = 0.7615941559557649.
/// assert!((efficiency(1.0).unwrap() - 1.0_f64.tanh()).abs() < 1e-12);
/// ```
pub fn efficiency(mlength: f64) -> Result<f64> {
    let x = FinError::non_negative("mL", mlength)?;
    if x < SMALL_ML {
        // tanh(x)/x = 1 - x^2/3 + 2 x^4/15 - ... ; two terms are far more
        // than enough below 1e-4 and avoid the 0/0 form at x == 0.
        Ok(1.0 - x * x / 3.0)
    } else {
        Ok(x.tanh() / x)
    }
}

/// Fin effectiveness `epsilon_f = sqrt(k P / (h A)) * tanh(mL)` (insulated tip).
///
/// Effectiveness is the ratio of the fin's heat rate to the heat that the
/// (now covered) base area would have rejected with no fin attached. Inputs
/// `h`, `perimeter` (P), `k`, and `area` (A) must be finite and strictly
/// positive; `mlength` (`m * L`) must be finite and non-negative.
///
/// A value `<= 1` means the fin does not help (it can even hurt); designers
/// typically want `epsilon_f >= 2`.
///
/// # Errors
///
/// Propagates the positivity checks on `h`, `perimeter`, `k`, `area` and the
/// non-negativity check on `mlength`.
///
/// # Example
///
/// ```
/// use valenx_fin::effectiveness;
/// // kP/hA = (200*0.1)/(100*1e-4) = 20/0.01 = 2000; sqrt = 44.72136.
/// // At mL = 1, tanh(1) = 0.76159..., so eps = 44.72136 * 0.76159 ~= 34.06.
/// let eps = effectiveness(100.0, 0.1, 200.0, 1.0e-4, 1.0).unwrap();
/// let expected = 2000.0_f64.sqrt() * 1.0_f64.tanh();
/// assert!((eps - expected).abs() < 1e-6);
/// ```
pub fn effectiveness(h: f64, perimeter: f64, k: f64, area: f64, mlength: f64) -> Result<f64> {
    let h = FinError::positive("h", h)?;
    let perimeter = FinError::positive("perimeter", perimeter)?;
    let k = FinError::positive("k", k)?;
    let area = FinError::positive("area", area)?;
    let x = FinError::non_negative("mL", mlength)?;
    Ok((k * perimeter / (h * area)).sqrt() * x.tanh())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn efficiency_zero_limit_is_one() {
        // Ground truth: mL -> 0 gives eta -> 1.
        assert!((efficiency(0.0).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn efficiency_small_ml_matches_series() {
        // Below the cross-over the Taylor form 1 - x^2/3 is used; check it
        // agrees with a hand-evaluated number at x = 5e-5.
        let x = 5.0e-5;
        let eta = efficiency(x).unwrap();
        assert!((eta - (1.0 - x * x / 3.0)).abs() < 1e-15);
        // ... and that it is essentially 1 there.
        assert!((eta - 1.0).abs() < 1e-8);
    }

    #[test]
    fn efficiency_at_one_matches_tanh() {
        // tanh(1) = 0.7615941559557649 (independent reference value).
        let eta = efficiency(1.0).unwrap();
        assert!((eta - 0.761_594_155_955_764_9).abs() < 1e-12);
    }

    #[test]
    fn efficiency_at_two_matches_reference() {
        // tanh(2)/2 = 0.9640275800758169 / 2 = 0.48201379003790845.
        let eta = efficiency(2.0).unwrap();
        assert!((eta - 0.482_013_790_037_908_4).abs() < 1e-12);
    }

    #[test]
    fn efficiency_large_ml_approaches_one_over_ml() {
        // Ground truth: large mL gives eta -> 1/(mL) because tanh -> 1.
        let x = 50.0;
        let eta = efficiency(x).unwrap();
        assert!((eta - 1.0 / x).abs() < 1e-9);
    }

    #[test]
    fn efficiency_is_continuous_across_crossover() {
        // The series branch and the tanh branch must agree at the seam.
        let below = efficiency(SMALL_ML * 0.999_999).unwrap();
        let above = efficiency(SMALL_ML * 1.000_001).unwrap();
        assert!((below - above).abs() < 1e-9);
    }

    #[test]
    fn efficiency_is_monotonically_decreasing() {
        // eta(mL) strictly decreases from 1; sample a few points.
        let a = efficiency(0.5).unwrap();
        let b = efficiency(1.0).unwrap();
        let c = efficiency(2.0).unwrap();
        assert!(a > b && b > c);
        assert!(a <= 1.0 && c > 0.0);
    }

    #[test]
    fn efficiency_rejects_negative_and_nan() {
        assert_eq!(efficiency(-0.1).unwrap_err().code(), "fin.negative");
        assert_eq!(efficiency(f64::NAN).unwrap_err().code(), "fin.not-finite");
    }

    #[test]
    fn effectiveness_matches_hand_worked_value() {
        // kP/hA = (200*0.1)/(100*1e-4) = 2000; sqrt(2000)*tanh(1).
        let eps = effectiveness(100.0, 0.1, 200.0, 1.0e-4, 1.0).unwrap();
        let expected = 2000.0_f64.sqrt() * 1.0_f64.tanh();
        assert!((eps - expected).abs() < 1e-6);
    }

    #[test]
    fn effectiveness_is_zero_at_zero_length() {
        // tanh(0) = 0  ->  epsilon_f = 0 (a fin of zero length adds nothing).
        let eps = effectiveness(100.0, 0.1, 200.0, 1.0e-4, 0.0).unwrap();
        assert!(eps.abs() < EPS);
    }

    #[test]
    fn effectiveness_rejects_bad_inputs() {
        assert!(effectiveness(0.0, 0.1, 200.0, 1e-4, 1.0).is_err());
        assert!(effectiveness(100.0, 0.1, 200.0, 1e-4, -1.0).is_err());
    }
}
