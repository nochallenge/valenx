//! Error taxonomy for `valenx-regression`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, RegressionError>`]. The variants are intentionally coarse —
//! a regression caller usually only cares about three things:
//!
//! 1. Did the two input vectors disagree on length
//!    ([`RegressionError::LengthMismatch`])?
//! 2. Did the caller pass too few points for the requested model — fewer
//!    than two points for a line, or fewer than `degree + 1` points for a
//!    degree-`degree` polynomial ([`RegressionError::TooFewPoints`])?
//! 3. Is the fit degenerate — every `x` identical (zero variance in the
//!    predictor) so the slope is undefined, or a rank-deficient normal
//!    matrix that the linear solve could not factor
//!    ([`RegressionError::Degenerate`])?
//!
//! Use [`RegressionError::code`] for stable log / telemetry tagging and
//! [`RegressionError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-springs`' `SpringsError`.

use thiserror::Error;

/// Errors produced by `valenx-regression`.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum RegressionError {
    /// The `x` and `y` input slices had different lengths. A regression
    /// needs one response value per predictor value.
    #[error("length mismatch: x has {x_len} points but y has {y_len}")]
    LengthMismatch {
        /// Number of predictor (`x`) values supplied.
        x_len: usize,
        /// Number of response (`y`) values supplied.
        y_len: usize,
    },

    /// Too few data points for the requested model. A simple line needs
    /// at least two points; a degree-`required - 1` polynomial needs
    /// `required` points (one per coefficient) to be determined.
    #[error("too few points: need at least {required}, got {got}")]
    TooFewPoints {
        /// Minimum number of points the model requires.
        required: usize,
        /// Number of points actually supplied.
        got: usize,
    },

    /// The fit is degenerate and has no unique solution: the predictor
    /// has zero variance (every `x` is identical, so a line has infinite
    /// slope) or the polynomial normal matrix is rank-deficient and could
    /// not be factored. `reason` is a human-readable explanation.
    #[error("degenerate fit: {reason}")]
    Degenerate {
        /// Human-readable reason surfaced in the UI / logs.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// Caller-supplied input is wrong: mismatched lengths or too few
    /// points.
    Input,
    /// The numerical problem is ill-posed: zero-variance predictor or a
    /// rank-deficient normal matrix.
    Degenerate,
}

impl RegressionError {
    /// Construct a [`RegressionError::LengthMismatch`].
    pub fn length_mismatch(x_len: usize, y_len: usize) -> Self {
        RegressionError::LengthMismatch { x_len, y_len }
    }

    /// Construct a [`RegressionError::TooFewPoints`].
    pub fn too_few_points(required: usize, got: usize) -> Self {
        RegressionError::TooFewPoints { required, got }
    }

    /// Construct a [`RegressionError::Degenerate`].
    pub fn degenerate(reason: impl Into<String>) -> Self {
        RegressionError::Degenerate {
            reason: reason.into(),
        }
    }

    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"regression.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            RegressionError::LengthMismatch { .. } => "regression.length_mismatch",
            RegressionError::TooFewPoints { .. } => "regression.too_few_points",
            RegressionError::Degenerate { .. } => "regression.degenerate",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            RegressionError::LengthMismatch { .. } | RegressionError::TooFewPoints { .. } => {
                ErrorCategory::Input
            }
            RegressionError::Degenerate { .. } => ErrorCategory::Degenerate,
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, RegressionError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = RegressionError::length_mismatch(3, 4);
        assert_eq!(err.code(), "regression.length_mismatch");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = RegressionError::too_few_points(2, 1);
        assert_eq!(err.code(), "regression.too_few_points");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = RegressionError::degenerate("zero variance in x");
        assert_eq!(err.code(), "regression.degenerate");
        assert_eq!(err.category(), ErrorCategory::Degenerate);
    }

    #[test]
    fn display_is_informative() {
        let msg = RegressionError::length_mismatch(3, 4).to_string();
        assert!(msg.contains('3') && msg.contains('4'), "got: {msg}");

        let msg = RegressionError::too_few_points(5, 2).to_string();
        assert!(msg.contains('5') && msg.contains('2'), "got: {msg}");

        let msg = RegressionError::degenerate("all x identical").to_string();
        assert!(msg.contains("all x identical"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(RegressionError::degenerate("x"));
        assert!(err.to_string().contains('x'));
    }
}
