//! The single error type returned by valenx-meshgen.

use std::path::PathBuf;

/// Everything that can go wrong loading or running a mesh-LLM.
#[derive(Debug, thiserror::Error)]
pub enum MeshGenError {
    /// No model file is configured (the user must download + point to one).
    #[error("no mesh-LLM is configured — download a model and set its path")]
    NoModelConfigured,

    /// The model file could not be loaded.
    #[error("failed to load model at {path}: {reason}")]
    ModelLoad {
        /// Path to the model file that failed to load.
        path: PathBuf,
        /// Human-readable reason the load failed.
        reason: String,
    },

    /// Inference failed mid-generation.
    #[error("inference failed: {0}")]
    Inference(String),

    /// The decoded geometry was unusable even after repair.
    #[error("the model produced an invalid mesh: {0}")]
    InvalidGeometry(String),

    /// An I/O error.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}
