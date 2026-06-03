//! Parts-library error taxonomy.

use thiserror::Error;

/// Errors raised by the parts library.
#[derive(Debug, Error)]
pub enum PartLibError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// I/O error wrapping a path.
    #[error("io ({path}): {reason}")]
    Io {
        /// Affected path (display string).
        path: String,
        /// Reason.
        reason: String,
    },

    /// Network fetch requested but v1 only supports local install.
    /// The plan is to route remote fetching through the Phase 22
    /// Add-on Manager pipeline (gh-release pattern) in a later phase.
    #[error("network fetch is not supported in v1 (url = {url})")]
    FetchRequiresNetwork {
        /// URL that was requested.
        url: String,
    },

    /// Part already installed under the same name.
    #[error("part `{name}` already installed")]
    DuplicatePart {
        /// Name that collided.
        name: String,
    },

    /// RON ser / de error wrapping the index file.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// I/O environment.
    Io,
    /// Network capability.
    Capability,
    /// Algorithm domain.
    Algorithm,
}

impl PartLibError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            PartLibError::BadParameter { .. } => "partlib.bad_parameter",
            PartLibError::Io { .. } => "partlib.io",
            PartLibError::FetchRequiresNetwork { .. } => "partlib.fetch_requires_network",
            PartLibError::DuplicatePart { .. } => "partlib.duplicate",
            PartLibError::Ron(_) => "partlib.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PartLibError::BadParameter { .. } => ErrorCategory::Input,
            PartLibError::Io { .. } => ErrorCategory::Io,
            PartLibError::FetchRequiresNetwork { .. } => ErrorCategory::Capability,
            PartLibError::DuplicatePart { .. } => ErrorCategory::Input,
            PartLibError::Ron(_) => ErrorCategory::Algorithm,
        }
    }
}
