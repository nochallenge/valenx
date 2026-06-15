//! Error taxonomy for BJT biasing.
//!
//! Every fallible constructor and analysis routine in this crate
//! returns a [`BjtError`]. The variants carry the offending parameter
//! name and a human-readable reason so a caller (or a GUI form) can
//! point at the bad field. [`BjtError::code`] gives a stable,
//! kebab-cased identifier for programmatic matching and
//! [`BjtError::category`] a coarse [`ErrorCategory`].

use thiserror::Error;

/// Errors raised while building BJT models or solving bias networks.
#[derive(Debug, Error)]
pub enum BjtError {
    /// A scalar parameter was outside its allowed domain (for example a
    /// non-positive resistance, a negative gain, or a `VBE` that is not
    /// finite).
    #[error("bad parameter `{name}`: {reason} (got {value})")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"beta"`, `"rc"`).
        name: &'static str,
        /// Why the value is rejected.
        reason: &'static str,
        /// The rejected value, echoed back for diagnostics.
        value: f64,
    },

    /// The network as specified has no forward-biased operating point —
    /// for example the Thevenin base voltage does not exceed `VBE`, so
    /// the base current would be zero or negative and the transistor is
    /// cut off rather than conducting.
    #[error("device does not conduct: {reason}")]
    CutOff {
        /// Why no conducting Q-point exists.
        reason: &'static str,
    },
}

impl BjtError {
    /// Construct a [`BjtError::BadParameter`].
    ///
    /// Internal helper used by the validated constructors; exposed so
    /// downstream code building its own derived parameters can raise a
    /// consistent error.
    pub fn bad_parameter(name: &'static str, reason: &'static str, value: f64) -> Self {
        BjtError::BadParameter {
            name,
            reason,
            value,
        }
    }

    /// Construct a [`BjtError::CutOff`].
    pub fn cut_off(reason: &'static str) -> Self {
        BjtError::CutOff { reason }
    }

    /// Stable kebab-cased identifier, suitable for logs or test matching.
    pub fn code(&self) -> &'static str {
        match self {
            BjtError::BadParameter { .. } => "bjt.bad_parameter",
            BjtError::CutOff { .. } => "bjt.cut_off",
        }
    }

    /// Coarse category of the failure.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BjtError::BadParameter { .. } => ErrorCategory::Input,
            BjtError::CutOff { .. } => ErrorCategory::Operating,
        }
    }
}

/// Coarse classification of a [`BjtError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// A user-supplied input value was invalid.
    Input,
    /// The network is valid but has no conducting operating point.
    Operating,
}
