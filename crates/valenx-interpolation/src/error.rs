//! Interpolation error taxonomy.
//!
//! Every fallible constructor in this crate funnels its rejection
//! through [`InterpolationError`] so callers get a stable
//! kebab-cased [`InterpolationError::code`] and a coarse
//! [`ErrorCategory`] for routing.

use thiserror::Error;

/// Errors raised while building or evaluating an interpolant.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum InterpolationError {
    /// No data points were supplied.
    ///
    /// At least one point is required to define any interpolant, and
    /// the cubic spline needs at least two.
    #[error("empty data: at least {needed} point(s) required, got 0")]
    EmptyData {
        /// Minimum number of points the requested method needs.
        needed: usize,
    },

    /// The `x` and `y` slices had different lengths.
    #[error("length mismatch: {x_len} x-values but {y_len} y-values")]
    LengthMismatch {
        /// Number of supplied abscissae.
        x_len: usize,
        /// Number of supplied ordinates.
        y_len: usize,
    },

    /// Not enough points for the requested method.
    ///
    /// For example a natural cubic spline needs at least two knots.
    #[error("too few points: method needs at least {needed}, got {got}")]
    TooFewPoints {
        /// Minimum the method requires.
        needed: usize,
        /// Number actually supplied.
        got: usize,
    },

    /// The abscissae were not in strictly ascending order.
    ///
    /// The offending index `i` is the first position whose `x` is not
    /// strictly greater than its predecessor.
    #[error("unsorted x at index {index}: x[{index}] = {x} is not > x[{prev_index}] = {prev}")]
    Unsorted {
        /// Index of the first non-ascending abscissa.
        index: usize,
        /// Index of its predecessor.
        prev_index: usize,
        /// The non-ascending value.
        x: f64,
        /// The predecessor value.
        prev: f64,
    },

    /// Two abscissae were equal (a vertical segment is undefined).
    ///
    /// Reported separately from [`InterpolationError::Unsorted`] so a
    /// caller can distinguish "your data needs sorting" from "your
    /// data has a duplicate `x`".
    #[error("duplicate x at index {index}: x[{index}] equals x[{prev_index}] = {x}")]
    DuplicateAbscissa {
        /// Index of the duplicate abscissa.
        index: usize,
        /// Index of the earlier equal abscissa.
        prev_index: usize,
        /// The repeated value.
        x: f64,
    },

    /// A non-finite (NaN or infinite) coordinate was supplied.
    #[error("non-finite coordinate in {axis} at index {index}")]
    NonFinite {
        /// Which axis the bad value sits on: `"x"` or `"y"`.
        axis: &'static str,
        /// Index of the offending coordinate.
        index: usize,
    },

    /// An evaluation abscissa fell outside the supported domain.
    ///
    /// All interpolants in this crate are defined only on the closed
    /// span `[x_first, x_last]`; extrapolation is rejected rather than
    /// silently returning a meaningless value.
    #[error("evaluation point {at} is outside the data domain [{lo}, {hi}]")]
    OutOfDomain {
        /// The requested evaluation abscissa.
        at: f64,
        /// Lower bound of the data domain (first abscissa).
        lo: f64,
        /// Upper bound of the data domain (last abscissa).
        hi: f64,
    },
}

/// Coarse error category, mirroring the taxonomy used by sibling
/// workbench crates.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied malformed input data.
    Input,
    /// The requested query lies outside the supported domain.
    Domain,
}

impl InterpolationError {
    /// Stable kebab-cased identifier, suitable for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            InterpolationError::EmptyData { .. } => "interpolation.empty_data",
            InterpolationError::LengthMismatch { .. } => "interpolation.length_mismatch",
            InterpolationError::TooFewPoints { .. } => "interpolation.too_few_points",
            InterpolationError::Unsorted { .. } => "interpolation.unsorted",
            InterpolationError::DuplicateAbscissa { .. } => "interpolation.duplicate_abscissa",
            InterpolationError::NonFinite { .. } => "interpolation.non_finite",
            InterpolationError::OutOfDomain { .. } => "interpolation.out_of_domain",
        }
    }

    /// Coarse category for routing / metrics.
    pub fn category(&self) -> ErrorCategory {
        match self {
            InterpolationError::OutOfDomain { .. } => ErrorCategory::Domain,
            _ => ErrorCategory::Input,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable_and_distinct() {
        let errs = [
            InterpolationError::EmptyData { needed: 1 },
            InterpolationError::LengthMismatch { x_len: 2, y_len: 3 },
            InterpolationError::TooFewPoints { needed: 2, got: 1 },
            InterpolationError::Unsorted {
                index: 1,
                prev_index: 0,
                x: 0.0,
                prev: 1.0,
            },
            InterpolationError::DuplicateAbscissa {
                index: 1,
                prev_index: 0,
                x: 1.0,
            },
            InterpolationError::NonFinite {
                axis: "x",
                index: 0,
            },
            InterpolationError::OutOfDomain {
                at: 9.0,
                lo: 0.0,
                hi: 1.0,
            },
        ];
        // All codes start with the crate prefix.
        for e in &errs {
            assert!(e.code().starts_with("interpolation."));
        }
        // All codes are distinct.
        for i in 0..errs.len() {
            for j in (i + 1)..errs.len() {
                assert_ne!(errs[i].code(), errs[j].code());
            }
        }
    }

    #[test]
    fn out_of_domain_is_the_only_domain_category() {
        let dom = InterpolationError::OutOfDomain {
            at: 9.0,
            lo: 0.0,
            hi: 1.0,
        };
        assert_eq!(dom.category(), ErrorCategory::Domain);

        let inp = InterpolationError::EmptyData { needed: 1 };
        assert_eq!(inp.category(), ErrorCategory::Input);
    }

    #[test]
    fn display_includes_indices() {
        let e = InterpolationError::Unsorted {
            index: 3,
            prev_index: 2,
            x: 1.0,
            prev: 2.0,
        };
        let s = format!("{e}");
        assert!(s.contains("index 3"));
        assert!(s.contains('2'));
    }
}
