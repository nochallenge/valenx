//! Catenary versus parabola: how far the two idealisations diverge.
//!
//! Both models describe the same physical cable; they differ only in how
//! the weight is distributed — along the arc (catenary) or along the
//! horizontal (parabola). For the *same* `w`, `L`, and `H` the catenary
//! sags slightly more, because per-arc-length weight piles up where the
//! cable is steepest near the supports.
//!
//! Expanding the catenary sag `a(cosh(L/2a) - 1)` in powers of
//! `t = L / 2a` gives `a (t^2/2 + t^4/24 + ...) = w L^2/(8H) (1 + t^2/12
//! + ...)`, whose leading term is exactly the parabolic sag
//! `w L^2 / (8 H)`. So the two agree as `t -> 0` (shallow cable) and the
//! relative gap grows like `t^2 / 12`.
//!
//! This module exposes both sags side by side and their ratio so a
//! caller can judge when the simpler parabola is good enough.

use crate::error::CableError;
use crate::{catenary, parabola};

/// Side-by-side sag comparison for a common `(w, L, H)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SagComparison {
    /// Parabolic mid-span sag, `w L^2 / (8 H)`.
    pub parabola_sag: f64,
    /// Catenary mid-span sag, `a (cosh(L/2a) - 1)`.
    pub catenary_sag: f64,
    /// Ratio `catenary_sag / parabola_sag` (>= 1, -> 1 as the cable
    /// becomes shallow).
    pub ratio: f64,
}

/// Compute both sags and their ratio for the same physical inputs.
///
/// Because the catenary always sags at least as much as the parabola,
/// the returned [`SagComparison::ratio`] is `>= 1` and approaches `1`
/// in the shallow-cable limit.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive
/// (delegated to the underlying model functions).
///
/// # Example
///
/// ```
/// use valenx_cablesag::compare::sag_comparison;
/// // Deep cable: a = 100, L = 100  ->  catenary sags ~2% more.
/// let c = sag_comparison(10.0, 100.0, 1000.0).unwrap();
/// assert!(c.catenary_sag > c.parabola_sag);
/// assert!(c.ratio > 1.0 && c.ratio < 1.05);
/// ```
pub fn sag_comparison(w: f64, span_l: f64, h_tension: f64) -> Result<SagComparison, CableError> {
    let parabola_sag = parabola::sag(w, span_l, h_tension)?;
    let catenary_sag = catenary::sag(w, span_l, h_tension)?;
    Ok(SagComparison {
        parabola_sag,
        catenary_sag,
        ratio: catenary_sag / parabola_sag,
    })
}

/// Leading-order relative sag gap, `t^2 / 12` with `t = L / 2a = w L / 2H`.
///
/// This is the first non-vanishing term of the relative excess
/// `(catenary_sag / parabola_sag) - 1`. It is a cheap estimate of how
/// much the parabola underestimates the catenary sag without evaluating
/// `cosh`, accurate for shallow cables.
///
/// # Errors
///
/// Returns [`CableError`] if any input is non-finite or non-positive.
///
/// # Example
///
/// ```
/// use valenx_cablesag::compare::leading_relative_sag_gap;
/// // t = w L / 2H = 10*100/2000 = 0.5  ->  gap ~ 0.25/12 ~ 0.0208.
/// let g = leading_relative_sag_gap(10.0, 100.0, 1000.0).unwrap();
/// assert!((g - 0.020_833_333).abs() < 1e-6);
/// ```
pub fn leading_relative_sag_gap(w: f64, span_l: f64, h_tension: f64) -> Result<f64, CableError> {
    let w = crate::error::positive("w", w)?;
    let l = crate::error::positive("span_l", span_l)?;
    let h = crate::error::positive("h_tension", h_tension)?;
    let t = w * l / (2.0 * h); // = L / 2a
    Ok(t * t / 12.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::CableError;

    #[test]
    fn catenary_sags_more_than_parabola() {
        let c = sag_comparison(10.0, 100.0, 1000.0).unwrap();
        assert!(c.catenary_sag > c.parabola_sag);
        assert!(c.ratio > 1.0);
    }

    #[test]
    fn ratio_matches_known_values() {
        // parabola: 10*100^2/(8*1000) = 12.5 m
        // catenary: 12.762596520638076 m
        // ratio: 1.0210077216510461
        let c = sag_comparison(10.0, 100.0, 1000.0).unwrap();
        assert!(
            (c.parabola_sag - 12.5).abs() < 1e-9,
            "got {}",
            c.parabola_sag
        );
        assert!(
            (c.catenary_sag - 12.762_596_520_638_076).abs() < 1e-9,
            "got {}",
            c.catenary_sag
        );
        assert!(
            (c.ratio - 1.021_007_721_651_046).abs() < 1e-9,
            "got {}",
            c.ratio
        );
    }

    #[test]
    fn ratio_tends_to_one_for_shallow_cable() {
        // High tension -> nearly straight -> ratio ~ 1. With
        // t = w L / 2H = 100/20000 = 5e-3 the gap is t^2/12 ~ 2.1e-6,
        // so the ratio sits just above 1. (A far larger H is avoided
        // here only because cosh(t) - 1 loses f64 precision for tiny t,
        // which is an arithmetic artefact, not a model effect.)
        let c = sag_comparison(1.0, 100.0, 1.0e4).unwrap();
        assert!(c.ratio > 1.0);
        assert!((c.ratio - 1.0).abs() < 1e-5, "got {}", c.ratio);
    }

    #[test]
    fn leading_gap_approximates_true_gap_for_shallow_cable() {
        // For a shallow cable the series gap t^2/12 should match the
        // actual (ratio - 1) closely.
        let w = 1.0;
        let l = 100.0;
        let h = 5000.0; // t = 100/10000 = 0.01, very shallow
        let g = leading_relative_sag_gap(w, l, h).unwrap();
        let c = sag_comparison(w, l, h).unwrap();
        let true_gap = c.ratio - 1.0;
        assert!((g - true_gap).abs() < 1e-6, "series {g} vs true {true_gap}");
    }

    #[test]
    fn leading_gap_known_value() {
        // t = 0.5  ->  t^2/12 = 0.25/12 = 0.0208333...
        let g = leading_relative_sag_gap(10.0, 100.0, 1000.0).unwrap();
        assert!((g - 0.020_833_333_333_333).abs() < 1e-9, "got {g}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(matches!(
            sag_comparison(-1.0, 100.0, 1000.0),
            Err(CableError::NonPositive { name: "w", .. })
        ));
        assert!(matches!(
            leading_relative_sag_gap(10.0, 100.0, 0.0),
            Err(CableError::NonPositive {
                name: "h_tension",
                ..
            })
        ));
    }
}
