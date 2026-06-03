//! CAM workbench error taxonomy.
//!
//! Same stable `code()` / `category()` pattern as `valenx-surface`,
//! `valenx-sketch`, etc. — every variant maps to a kebab-cased string
//! identifier that LLM / scripting layers can branch on without
//! parsing the human-readable message.

use thiserror::Error;

/// Errors raised by tool construction, operation generation,
/// postprocessing, or persistence.
#[derive(Debug, Error)]
pub enum CamError {
    /// Stock block has zero or negative extent in one or more axes.
    #[error("empty stock")]
    EmptyStock,

    /// Tool failed validation (e.g. zero diameter, zero flutes).
    #[error("bad tool: {reason}")]
    BadTool {
        /// Human-readable reason — included in the displayed message.
        reason: String,
    },

    /// Per-operation parameters failed validation (e.g. negative
    /// step-down, step-over > tool diameter * 0.5).
    #[error("bad operation '{name}': {reason}")]
    BadOperation {
        /// Operation name (e.g. `"profile"`, `"pocket"`).
        name: String,
        /// Why the operation could not be generated.
        reason: String,
    },

    /// Postprocessor tried to emit G-code for an empty toolpath.
    #[error("empty toolpath")]
    EmptyToolpath,

    /// Postprocessor refused to format a move (e.g. cut with feed=0).
    #[error("postprocessor failed: {reason}")]
    PostprocessorFailed {
        /// Reason supplied by the failing postprocessor.
        reason: String,
    },

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON serialise / parse error.
    #[error("ron: {0}")]
    Ron(String),
}

impl CamError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            CamError::EmptyStock => "cam.empty_stock",
            CamError::BadTool { .. } => "cam.bad_tool",
            CamError::BadOperation { .. } => "cam.bad_operation",
            CamError::EmptyToolpath => "cam.empty_toolpath",
            CamError::PostprocessorFailed { .. } => "cam.postprocessor_failed",
            CamError::Io(_) => "cam.io",
            CamError::Ron(_) => "cam.ron",
        }
    }

    /// Coarse classification — `"input"`, `"geometry"`, `"post"`, or
    /// `"io"`.
    pub fn category(&self) -> &'static str {
        match self {
            CamError::EmptyStock | CamError::BadTool { .. } | CamError::BadOperation { .. } => {
                "input"
            }
            CamError::EmptyToolpath => "geometry",
            CamError::PostprocessorFailed { .. } => "post",
            CamError::Io(_) | CamError::Ron(_) => "io",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code() {
        let cases: Vec<(CamError, &str, &str)> = vec![
            (CamError::EmptyStock, "cam.empty_stock", "input"),
            (
                CamError::BadTool {
                    reason: "diameter must be > 0".into(),
                },
                "cam.bad_tool",
                "input",
            ),
            (
                CamError::BadOperation {
                    name: "profile".into(),
                    reason: "step_down must be > 0".into(),
                },
                "cam.bad_operation",
                "input",
            ),
            (CamError::EmptyToolpath, "cam.empty_toolpath", "geometry"),
            (
                CamError::PostprocessorFailed {
                    reason: "cut feed = 0".into(),
                },
                "cam.postprocessor_failed",
                "post",
            ),
            (CamError::Ron("malformed".into()), "cam.ron", "io"),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn io_wraps_std() {
        let e: CamError = std::io::Error::other("disk gone").into();
        assert_eq!(e.code(), "cam.io");
        assert_eq!(e.category(), "io");
    }
}
