//! Curves workbench error taxonomy.

use thiserror::Error;

/// Errors raised by curve operations.
#[derive(Debug, Error)]
pub enum CurvesError {
    /// Caller passed a degenerate input — empty point list, zero-length
    /// curve, …
    #[error("degenerate input: {0}")]
    Degenerate(String),

    /// User-supplied parameter out of range (negative distance,
    /// t_start >= t_end, …).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// Wraps an error from `valenx-surface` (typically curve
    /// construction validation).
    #[error("surface: {0}")]
    Surface(String),

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON error.
    #[error("ron: {0}")]
    Ron(String),
}

impl From<valenx_surface::error::SurfaceError> for CurvesError {
    fn from(e: valenx_surface::error::SurfaceError) -> Self {
        CurvesError::Surface(e.to_string())
    }
}

/// Coarse category an LLM / UI can branch on.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input is bad.
    Input,
    /// User knob out of range.
    Config,
    /// Bug in valenx-curves or in a delegated library.
    Internal,
    /// IO / RON transient.
    Runtime,
}

impl CurvesError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            CurvesError::Degenerate(_) => "curves.degenerate",
            CurvesError::BadParameter { .. } => "curves.bad_parameter",
            CurvesError::Surface(_) => "curves.surface",
            CurvesError::Io(_) => "curves.io",
            CurvesError::Ron(_) => "curves.ron",
        }
    }

    /// High-level classification.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CurvesError::Degenerate(_) => ErrorCategory::Input,
            CurvesError::BadParameter { .. } => ErrorCategory::Config,
            CurvesError::Surface(_) => ErrorCategory::Internal,
            CurvesError::Io(_) | CurvesError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(CurvesError, &'static str, ErrorCategory)> = vec![
            (
                CurvesError::Degenerate("zero len".into()),
                "curves.degenerate",
                ErrorCategory::Input,
            ),
            (
                CurvesError::BadParameter {
                    name: "distance",
                    reason: "zero".into(),
                },
                "curves.bad_parameter",
                ErrorCategory::Config,
            ),
            (
                CurvesError::Surface("upstream".into()),
                "curves.surface",
                ErrorCategory::Internal,
            ),
            (
                CurvesError::Ron("bad token".into()),
                "curves.ron",
                ErrorCategory::Runtime,
            ),
        ];
        for (err, code, cat) in cases {
            assert_eq!(err.code(), code);
            assert_eq!(err.category(), cat);
        }
    }
}
