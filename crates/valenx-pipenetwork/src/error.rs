//! Error taxonomy for the pipe-network solver.

use thiserror::Error;

/// Errors raised by the pipe-network crate.
///
/// Every fallible constructor in this crate validates its inputs up front
/// and returns one of these variants rather than panicking. Float fields
/// are stored as `f64`; the variants below describe *why* a value was
/// rejected, not the value itself.
#[derive(Debug, Error)]
pub enum NetworkError {
    /// A scalar parameter failed validation (non-finite, or outside the
    /// physically meaningful range documented for that parameter).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// The parameter that was rejected.
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A loop referenced a pipe index that does not exist in the network.
    #[error("loop references unknown pipe index {index} (network has {count} pipes)")]
    UnknownPipe {
        /// The out-of-range index that was referenced.
        index: usize,
        /// The number of pipes actually present in the network.
        count: usize,
    },

    /// A loop was constructed with no pipes. A Hardy-Cross loop must
    /// contain at least one pipe for the correction to be defined.
    #[error("loop `{0}` is empty (a loop needs at least one pipe)")]
    EmptyLoop(String),

    /// The Hardy-Cross iteration did not converge within the configured
    /// iteration budget. The largest remaining absolute loop correction
    /// at the point of giving up is reported.
    #[error(
        "Hardy-Cross did not converge in {iterations} iterations \
         (largest remaining |dQ| = {residual}, tolerance = {tolerance})"
    )]
    NoConvergence {
        /// Number of iterations performed before giving up.
        iterations: usize,
        /// Largest remaining absolute loop correction magnitude.
        residual: f64,
        /// The convergence tolerance that was not met.
        tolerance: f64,
    },
}

/// Coarse category of a [`NetworkError`], useful for UI grouping and
/// for deciding whether a failure is the caller's fault (input) or the
/// solver's (algorithm).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// Bad caller-supplied input (parameters, indices, topology).
    Input,
    /// The numerical algorithm failed to reach a solution.
    Algorithm,
}

impl NetworkError {
    /// Stable, kebab-cased identifier for this error, suitable for logs
    /// and telemetry where the human-readable message may change.
    pub fn code(&self) -> &'static str {
        match self {
            NetworkError::BadParameter { .. } => "pipenetwork.bad_parameter",
            NetworkError::UnknownPipe { .. } => "pipenetwork.unknown_pipe",
            NetworkError::EmptyLoop(_) => "pipenetwork.empty_loop",
            NetworkError::NoConvergence { .. } => "pipenetwork.no_convergence",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            NetworkError::BadParameter { .. }
            | NetworkError::UnknownPipe { .. }
            | NetworkError::EmptyLoop(_) => ErrorCategory::Input,
            NetworkError::NoConvergence { .. } => ErrorCategory::Algorithm,
        }
    }

    /// Construct a [`NetworkError::BadParameter`] from a static name and an
    /// owned reason string. Small helper to keep call sites terse.
    pub(crate) fn bad(name: &'static str, reason: impl Into<String>) -> Self {
        NetworkError::BadParameter {
            name,
            reason: reason.into(),
        }
    }
}
