//! Error taxonomy for strain-rosette reduction.
//!
//! Every fallible entry point in the crate returns
//! [`Result<T, RosetteError>`]. Errors carry a stable kebab-cased
//! [`RosetteError::code`] and a coarse [`ErrorCategory`] so callers can
//! branch on the failure class without string-matching the `Display`
//! text.

use thiserror::Error;

/// Errors raised while validating material properties or reducing a
/// rosette.
#[derive(Debug, Error)]
pub enum RosetteError {
    /// A material constant fell outside its physically admissible range
    /// (for example a non-positive Young's modulus).
    #[error("invalid material parameter `{name}`: {reason} (got {value})")]
    InvalidMaterial {
        /// The offending parameter name (`"youngs_modulus"`,
        /// `"poisson_ratio"`).
        name: &'static str,
        /// Why the value is rejected.
        reason: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A supplied number was not finite (`NaN` or an infinity), which
    /// would silently poison every downstream computation.
    #[error("non-finite value for `{name}`: {value}")]
    NonFinite {
        /// The offending field name.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },
}

/// Coarse classification of a [`RosetteError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-range or non-finite input.
    Input,
    /// A tunable material constant was rejected.
    Material,
}

impl RosetteError {
    /// Stable kebab-cased identifier for programmatic matching.
    ///
    /// The string is part of the crate's public contract and will not
    /// change for a given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            RosetteError::InvalidMaterial { .. } => "rosette.invalid-material",
            RosetteError::NonFinite { .. } => "rosette.non-finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RosetteError::InvalidMaterial { .. } => ErrorCategory::Material,
            RosetteError::NonFinite { .. } => ErrorCategory::Input,
        }
    }
}

/// Internal helper: reject a non-finite input with a [`RosetteError`].
///
/// Used by the validated constructors so the same guard wording is
/// shared everywhere.
pub(crate) fn ensure_finite(name: &'static str, value: f64) -> Result<(), RosetteError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RosetteError::NonFinite { name, value })
    }
}
