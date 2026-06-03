//! Typed errors for `valenx-subdivision`.

use thiserror::Error;

/// Errors raised by the Wings-style subdivision workbench.
#[derive(Debug, Error)]
pub enum SubdivError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Mesh topology problem — e.g. non-quad face in Catmull-Clark
    /// without the Loop fallback, or non-triangle face in Loop.
    #[error("topology error: {0}")]
    Topology(String),

    /// Index out of range.
    #[error("index out of range: {kind} {idx} (limit {limit})")]
    IndexOutOfRange {
        /// What kind of index — `"face"` / `"edge"` / `"vertex"`.
        kind: &'static str,
        /// Bad index.
        idx: usize,
        /// Upper bound (exclusive).
        limit: usize,
    },
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input / parameter.
    Input,
    /// Topology / mesh structure problem.
    Topology,
}

impl SubdivError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            SubdivError::BadParameter { .. } => "subdiv.bad_parameter",
            SubdivError::Topology(_) => "subdiv.topology",
            SubdivError::IndexOutOfRange { .. } => "subdiv.index_out_of_range",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            SubdivError::BadParameter { .. } => ErrorCategory::Input,
            SubdivError::Topology(_) | SubdivError::IndexOutOfRange { .. } => {
                ErrorCategory::Topology
            }
        }
    }
}
