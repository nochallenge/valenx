//! Typed errors for `valenx-blender-mesh-ops`.

use thiserror::Error;

/// Errors raised by the Blender-style mesh ops.
#[derive(Debug, Error)]
pub enum BlenderOpError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Index out of range.
    #[error("index out of range: {kind} {idx} (limit {limit})")]
    IndexOutOfRange {
        /// What kind — `"face"` / `"edge"` / `"vertex"`.
        kind: &'static str,
        /// Bad index.
        idx: usize,
        /// Upper bound (exclusive).
        limit: usize,
    },

    /// Topology problem.
    #[error("topology: {0}")]
    Topology(String),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Topology / mesh structure.
    Topology,
}

impl BlenderOpError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            BlenderOpError::BadParameter { .. } => "blender-op.bad_parameter",
            BlenderOpError::IndexOutOfRange { .. } => "blender-op.index_out_of_range",
            BlenderOpError::Topology(_) => "blender-op.topology",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BlenderOpError::BadParameter { .. } => ErrorCategory::Input,
            BlenderOpError::IndexOutOfRange { .. } | BlenderOpError::Topology(_) => {
                ErrorCategory::Topology
            }
        }
    }
}
