//! Pipe-flow error taxonomy.
//!
//! Every fallible entry point in this crate validates its inputs up
//! front and returns a [`PipeFlowError`] describing exactly which
//! physical quantity was out of range, rather than silently producing a
//! `NaN`/`Inf` or a negative friction factor downstream.

use thiserror::Error;

/// Errors raised by the pipe-flow correlations.
#[derive(Debug, Error)]
pub enum PipeFlowError {
    /// A physical input was non-positive where a strictly positive value
    /// is required (density, velocity, diameter, viscosity, length, …).
    #[error("non-positive `{name}`: expected > 0, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A dimensionless input fell outside its admissible range.
    #[error("`{name}` out of range: expected {expected}, got {value}")]
    OutOfRange {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable description of the admissible range.
        expected: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A non-finite (`NaN` or infinite) input was supplied.
    #[error("`{name}` is not finite: got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },
}

/// Coarse error category, mirroring the taxonomy used across sibling
/// Valenx crates.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied a bad input value.
    Input,
}

impl PipeFlowError {
    /// Stable, kebab-cased identifier for the error variant. Useful for
    /// logging and for matching in tests without depending on the
    /// human-readable [`std::fmt::Display`] string.
    pub fn code(&self) -> &'static str {
        match self {
            PipeFlowError::NonPositive { .. } => "pipeflow.non_positive",
            PipeFlowError::OutOfRange { .. } => "pipeflow.out_of_range",
            PipeFlowError::NotFinite { .. } => "pipeflow.not_finite",
        }
    }

    /// Coarse category. Every variant in this crate is an input error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PipeFlowError::NonPositive { .. }
            | PipeFlowError::OutOfRange { .. }
            | PipeFlowError::NotFinite { .. } => ErrorCategory::Input,
        }
    }
}

/// Validate that `value` is finite and strictly positive, returning it
/// unchanged on success. Used by the public constructors so every
/// physical quantity is screened at the boundary.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, PipeFlowError> {
    if !value.is_finite() {
        return Err(PipeFlowError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(PipeFlowError::NonPositive { name, value });
    }
    Ok(value)
}

/// Validate that `value` is finite and lies in the closed interval
/// `[lo, hi]`, returning it unchanged on success.
pub(crate) fn require_in_closed(
    name: &'static str,
    value: f64,
    lo: f64,
    hi: f64,
    expected: &'static str,
) -> Result<f64, PipeFlowError> {
    if !value.is_finite() {
        return Err(PipeFlowError::NotFinite { name, value });
    }
    if value < lo || value > hi {
        return Err(PipeFlowError::OutOfRange {
            name,
            expected,
            value,
        });
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        let v = require_positive("d", 0.25).expect("0.25 is positive");
        assert!((v - 0.25).abs() < 1e-12);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        let zero = require_positive("d", 0.0).unwrap_err();
        assert_eq!(zero.code(), "pipeflow.non_positive");
        let neg = require_positive("d", -1.0).unwrap_err();
        assert_eq!(neg.code(), "pipeflow.non_positive");
        assert_eq!(neg.category(), ErrorCategory::Input);
    }

    #[test]
    fn require_positive_rejects_nan_and_inf() {
        let nan = require_positive("d", f64::NAN).unwrap_err();
        assert_eq!(nan.code(), "pipeflow.not_finite");
        let inf = require_positive("d", f64::INFINITY).unwrap_err();
        assert_eq!(inf.code(), "pipeflow.not_finite");
    }

    #[test]
    fn require_in_closed_enforces_bounds() {
        assert!(require_in_closed("x", 0.0, 0.0, 1.0, "[0,1]").is_ok());
        assert!(require_in_closed("x", 1.0, 0.0, 1.0, "[0,1]").is_ok());
        let too_big = require_in_closed("x", 1.5, 0.0, 1.0, "[0,1]").unwrap_err();
        assert_eq!(too_big.code(), "pipeflow.out_of_range");
        let too_small = require_in_closed("x", -0.1, 0.0, 1.0, "[0,1]").unwrap_err();
        assert_eq!(too_small.code(), "pipeflow.out_of_range");
    }
}
