//! Inclined-plane error taxonomy.
//!
//! Every fallible constructor in this crate returns
//! [`InclinedPlaneError`]. The variants carry stable
//! [`code`](InclinedPlaneError::code) and
//! [`category`](InclinedPlaneError::category) accessors so callers can
//! branch or log without matching on the human-readable message.

use thiserror::Error;

/// Errors raised when validating ramp parameters.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum InclinedPlaneError {
    /// A scalar parameter fell outside its valid range.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// The offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// Why the supplied value was rejected.
        reason: String,
    },
}

impl InclinedPlaneError {
    /// Construct a [`InclinedPlaneError::BadParameter`].
    ///
    /// `name` is the stable parameter identifier; `reason` is a
    /// human-readable explanation of why the value was rejected.
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        InclinedPlaneError::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Stable kebab-cased identifier for telemetry / matching.
    pub fn code(&self) -> &'static str {
        match self {
            InclinedPlaneError::BadParameter { .. } => "inclinedplane.bad_parameter",
        }
    }

    /// Coarse error category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            InclinedPlaneError::BadParameter { .. } => ErrorCategory::Input,
        }
    }
}

/// Coarse classification of an [`InclinedPlaneError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied input was invalid.
    Input,
}
