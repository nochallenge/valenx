//! Belt-drive error taxonomy.
//!
//! A single [`BeltError`] enum covers every failure mode in the crate.
//! Inputs are validated up front by the constructors in the topic
//! modules ([`crate::geometry`], [`crate::friction`], [`crate::power`]),
//! so a returned [`BeltError`] always names the offending parameter and
//! the reason it was rejected.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors raised by belt-drive calculations.
#[derive(Debug, Error)]
pub enum BeltError {
    /// A scalar parameter fell outside its physically valid domain
    /// (for example a non-positive diameter, a negative coefficient of
    /// friction, or a wrap angle that is not strictly positive).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the rejected parameter.
        name: &'static str,
        /// Human-readable explanation of why it was rejected.
        reason: String,
    },

    /// The supplied geometry cannot form a real open-belt drive — for
    /// instance the centre distance is too small to clear both pulleys,
    /// which would require an imaginary tangent length.
    #[error("degenerate geometry: {0}")]
    DegenerateGeometry(String),
}

/// Coarse classification of a [`BeltError`], useful for routing a
/// failure to a user-facing message versus an internal log without
/// matching on every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ErrorCategory {
    /// The caller passed an out-of-domain value.
    Input,
    /// The combination of otherwise-valid inputs has no real solution.
    Geometry,
}

impl BeltError {
    /// Stable, kebab-cased identifier for this error, suitable for
    /// logging or matching in tests without depending on the
    /// human-readable [`Display`](std::fmt::Display) text.
    pub fn code(&self) -> &'static str {
        match self {
            BeltError::BadParameter { .. } => "beltdrive.bad_parameter",
            BeltError::DegenerateGeometry(_) => "beltdrive.degenerate_geometry",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BeltError::BadParameter { .. } => ErrorCategory::Input,
            BeltError::DegenerateGeometry(_) => ErrorCategory::Geometry,
        }
    }

    /// Construct a [`BeltError::BadParameter`] from a static parameter
    /// name and an owned reason string. Used by the validated
    /// constructors throughout the crate.
    pub(crate) fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        BeltError::BadParameter {
            name,
            reason: reason.into(),
        }
    }
}
