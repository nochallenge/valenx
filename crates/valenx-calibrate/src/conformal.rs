//! Split-conformal prediction: distribution-free intervals with a
//! finite-sample coverage guarantee.

use crate::error::CalibrateError;

/// The split-conformal quantile of a set of held-out nonconformity scores
/// (e.g. absolute residuals) at miscoverage level `alpha`.
///
/// Returns the `k`-th smallest score, where `k = ceil((n + 1)(1 - alpha))`.
/// Adding this radius around a point prediction yields an interval that covers
/// a new exchangeable point with probability at least `1 - alpha`. When
/// `k > n` (too few calibration points for the requested coverage) the largest
/// score is returned and the `1 - alpha` guarantee cannot be met — widen the
/// calibration set or raise `alpha`.
pub fn split_conformal_quantile(residuals: &[f64], alpha: f64) -> Result<f64, CalibrateError> {
    if residuals.is_empty() {
        return Err(CalibrateError::Empty { what: "residuals" });
    }
    if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 {
        return Err(CalibrateError::AlphaOutOfRange { alpha });
    }
    for &r in residuals {
        if !r.is_finite() {
            return Err(CalibrateError::NonFinite { what: "residual" });
        }
    }
    let n = residuals.len();
    let mut sorted: Vec<f64> = residuals.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let k = ((n as f64 + 1.0) * (1.0 - alpha)).ceil() as usize;
    let idx = k.clamp(1, n) - 1;
    Ok(sorted[idx])
}

/// A symmetric split-conformal interval `(prediction - q, prediction + q)`,
/// where `q` is the [`split_conformal_quantile`] of the calibration
/// `residuals`.
pub fn conformal_interval(
    prediction: f64,
    residuals: &[f64],
    alpha: f64,
) -> Result<(f64, f64), CalibrateError> {
    if !prediction.is_finite() {
        return Err(CalibrateError::NonFinite { what: "prediction" });
    }
    let q = split_conformal_quantile(residuals, alpha)?;
    Ok((prediction - q, prediction + q))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_matches_rank_formula() {
        let r: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        // k = ceil(11 * 0.9) = 10 -> 10th smallest = 10
        assert!((split_conformal_quantile(&r, 0.1).unwrap() - 10.0).abs() < 1e-12);
        // k = ceil(11 * 0.5) = 6 -> 6th smallest = 6
        assert!((split_conformal_quantile(&r, 0.5).unwrap() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn interval_is_symmetric_around_prediction() {
        let r = vec![1.0; 10];
        let (lo, hi) = conformal_interval(5.0, &r, 0.1).unwrap();
        assert!((lo - 4.0).abs() < 1e-12);
        assert!((hi - 6.0).abs() < 1e-12);
    }

    #[test]
    fn higher_coverage_gives_wider_radius() {
        let r: Vec<f64> = (1..=100).map(|i| i as f64).collect();
        let tight = split_conformal_quantile(&r, 0.5).unwrap();
        let wide = split_conformal_quantile(&r, 0.05).unwrap();
        assert!(wide >= tight);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            split_conformal_quantile(&[], 0.1).unwrap_err().code(),
            "empty"
        );
        assert_eq!(
            split_conformal_quantile(&[1.0], 0.0).unwrap_err().code(),
            "alpha_out_of_range"
        );
        assert_eq!(
            split_conformal_quantile(&[1.0], 1.0).unwrap_err().code(),
            "alpha_out_of_range"
        );
    }
}
