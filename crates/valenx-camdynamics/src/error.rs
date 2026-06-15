//! Cam-dynamics error taxonomy.

use thiserror::Error;

/// Errors raised when constructing or evaluating a cam motion law.
#[derive(Debug, Error)]
pub enum CamError {
    /// A scalar parameter was outside its valid domain.
    ///
    /// Raised, for example, when the lift is negative or the rise angle
    /// `beta` is not strictly positive.
    #[error("bad parameter `{name}`: {reason} (got {value})")]
    BadParameter {
        /// Parameter name (a stable, kebab-friendly identifier).
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: &'static str,
        /// The offending value, formatted for display.
        value: f64,
    },

    /// A parameter was required to be finite but was `NaN` or infinite.
    #[error("parameter `{name}` must be finite (got {value})")]
    NotFinite {
        /// Parameter name.
        name: &'static str,
        /// The offending value, formatted for display.
        value: f64,
    },
}

impl CamError {
    /// A stable kebab-cased identifier for this error variant.
    ///
    /// Intended for logging and programmatic matching; it never changes
    /// for a given variant even if the human-readable message does.
    pub fn code(&self) -> &'static str {
        match self {
            CamError::BadParameter { .. } => "camdynamics.bad_parameter",
            CamError::NotFinite { .. } => "camdynamics.not_finite",
        }
    }

    /// The coarse category this error belongs to.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CamError::BadParameter { .. } => ErrorCategory::Input,
            CamError::NotFinite { .. } => ErrorCategory::Input,
        }
    }

    /// Construct a [`CamError::BadParameter`] for a non-positive value.
    ///
    /// This is a small helper used by the validated constructors so the
    /// reason string stays consistent across call sites.
    pub(crate) fn non_positive(name: &'static str, value: f64) -> Self {
        CamError::BadParameter {
            name,
            reason: "must be strictly positive",
            value,
        }
    }

    /// Construct a [`CamError::BadParameter`] for a negative value.
    pub(crate) fn negative(name: &'static str, value: f64) -> Self {
        CamError::BadParameter {
            name,
            reason: "must be non-negative",
            value,
        }
    }

    /// Construct a [`CamError::NotFinite`] for a `NaN` / infinite value.
    pub(crate) fn not_finite(name: &'static str, value: f64) -> Self {
        CamError::NotFinite { name, value }
    }
}

/// Coarse classification of a [`CamError`].
///
/// Useful for deciding how to surface an error: an [`ErrorCategory::Input`]
/// problem is the caller's to fix, whereas other categories would signal
/// an internal contract violation.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-domain or non-finite input.
    Input,
    /// A tunable configuration knob was invalid.
    Config,
    /// An algorithm-internal domain violation occurred.
    Algorithm,
}
