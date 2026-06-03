//! Error taxonomy shared by every adapter.
//!
//! The UI renders these differently — `ToolNotInstalled` opens a
//! "Install now" dialog; `InvalidCase` highlights the offending
//! parameter inline; `Run` opens the log viewer. Keeping the enum
//! structured (not just strings) is what makes that possible.
//!
//! Spec: RFC 0002 § Error classification.

use std::path::PathBuf;

use semver::Version;
use thiserror::Error;

/// Coarse phase within a solver run, for reporting where a failure
/// happened.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunPhase {
    Startup,
    MeshRead,
    Solve,
    Output,
    Shutdown,
}

/// Error that occurred while translating a canonical case into an
/// adapter's native input deck.
#[derive(Debug, Error)]
pub enum TranslateError {
    #[error("unsupported feature: {feature}")]
    Unsupported { feature: String },
    #[error("required field missing: {field}")]
    MissingField { field: String },
    #[error("value {value} for {field} out of range [{min}, {max}]")]
    OutOfRange {
        field: String,
        value: f64,
        min: f64,
        max: f64,
    },
    #[error("inconsistent case: {0}")]
    Inconsistent(String),
}

/// The structured error type adapters return.
#[derive(Debug, Error)]
pub enum AdapterError {
    /// The required external tool is not installed or not on PATH.
    #[error("tool {name} is not installed ({hint})")]
    ToolNotInstalled { name: &'static str, hint: String },

    /// The installed tool's version is outside the adapter's supported
    /// range.
    #[error("{name} version mismatch: expected {expected}, found {found}")]
    ToolVersionMismatch {
        name: &'static str,
        expected: String,
        found: Version,
    },

    /// The case the user authored is invalid.
    #[error("invalid case at {case_path}: {reason}")]
    InvalidCase { case_path: PathBuf, reason: String },

    /// Case-to-native translation failed.
    #[error("translation failed: {0}")]
    Translate(#[from] TranslateError),

    /// The solver subprocess returned a non-zero exit code.
    #[error("solver exited {exit_code} during {phase:?}: {stderr}")]
    Run {
        exit_code: i32,
        stderr: String,
        phase: RunPhase,
    },

    /// The adapter could not parse one of the tool's output files.
    #[error("failed to parse output {file}: {reason}")]
    ParseOutput { file: PathBuf, reason: String },

    /// The caller asked for cancellation; the adapter honored it.
    #[error("cancelled by caller")]
    Cancelled,

    /// Underlying I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Anything else — kept out of the taxonomy.
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}
