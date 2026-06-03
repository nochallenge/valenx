//! Typed errors for `valenx-paramhist`.

use thiserror::Error;

/// Errors raised by the NaroCAD-style history.
#[derive(Debug, Error)]
pub enum ParamHistError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Cycle detected in the dependency graph at the given entry.
    #[error("dependency cycle including entry {0}")]
    Cycle(usize),

    /// Index out of range.
    #[error("index out of range: {idx} (limit {limit})")]
    IndexOutOfRange {
        /// Bad index.
        idx: usize,
        /// Upper bound (exclusive).
        limit: usize,
    },

    /// DAG-preserving op refused — the move would create a back-edge.
    #[error("move would violate DAG: {0}")]
    InvalidMove(String),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// DAG / topology.
    Topology,
}

impl ParamHistError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            ParamHistError::BadParameter { .. } => "paramhist.bad_parameter",
            ParamHistError::Cycle(_) => "paramhist.cycle",
            ParamHistError::IndexOutOfRange { .. } => "paramhist.index_out_of_range",
            ParamHistError::InvalidMove(_) => "paramhist.invalid_move",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ParamHistError::BadParameter { .. } => ErrorCategory::Input,
            ParamHistError::Cycle(_)
            | ParamHistError::IndexOutOfRange { .. }
            | ParamHistError::InvalidMove(_) => ErrorCategory::Topology,
        }
    }
}
