//! Gears workbench error taxonomy.

use thiserror::Error;

/// Errors raised by gear generation.
#[derive(Debug, Error)]
pub enum GearsError {
    /// Bad parameter (non-positive module, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Tooth count is below the practical minimum (5).
    #[error("teeth count too low: got {0}, need >= 5")]
    TeethTooFew(u32),

    /// Gear kind is not implemented yet.
    #[error("kind `{0}` not yet supported")]
    KindNotImplemented(&'static str),

    /// Geometric inversion (helical wedge / worm self-intersection).
    #[error("degenerate geometry: {0}")]
    Degenerate(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Algorithm domain.
    Algorithm,
}

impl GearsError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            GearsError::BadParameter { .. } => "gears.bad_parameter",
            GearsError::TeethTooFew(_) => "gears.teeth_too_few",
            GearsError::KindNotImplemented(_) => "gears.kind_not_implemented",
            GearsError::Degenerate(_) => "gears.degenerate",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            GearsError::BadParameter { .. } => ErrorCategory::Config,
            GearsError::TeethTooFew(_) => ErrorCategory::Input,
            GearsError::KindNotImplemented(_) | GearsError::Degenerate(_) => {
                ErrorCategory::Algorithm
            }
        }
    }
}
