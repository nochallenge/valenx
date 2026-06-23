//! The crate-wide error type.

use thiserror::Error;

/// Errors produced by the UQ toolkit.
///
/// Constructors and routines that can be misused (an invalid distribution, a
/// percentile outside `[0, 100]`, mismatched sample/value lengths, a
/// rank-deficient surrogate fit, …) return [`Result<_, UqError>`] rather than
/// panicking, so a bad input is a recoverable error, never a crash or a `NaN`.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum UqError {
    /// A distribution parameter was out of range — e.g. a non-positive
    /// standard deviation, `lo >= hi`, or a triangular `mode` outside
    /// `[lo, hi]`.
    #[error("invalid distribution parameter: {0}")]
    InvalidDistribution(String),

    /// A percentile / probability level was outside its valid range.
    #[error("value out of range: {0}")]
    OutOfRange(String),

    /// An input was empty where at least one element is required (for example,
    /// computing statistics over an empty sample).
    #[error("empty input: {0}")]
    EmptyInput(String),

    /// Two slices that must have matching lengths did not — e.g. the number of
    /// sample rows differs from the number of observed values, or a sample
    /// row's dimension does not match the model / distribution count.
    #[error("dimension mismatch: {0}")]
    DimensionMismatch(String),

    /// A linear-algebra step failed: the least-squares system for the
    /// surrogate could not be solved (typically a rank-deficient design
    /// matrix — too few or collinear samples for the requested polynomial
    /// degree).
    #[error("linear-algebra failure: {0}")]
    LinearAlgebra(String),
}
