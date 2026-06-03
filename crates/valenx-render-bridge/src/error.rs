//! Render-bridge error taxonomy.

use thiserror::Error;

/// Errors raised by scene-file emission.
#[derive(Debug, Error)]
pub enum RenderError {
    /// Scene has no geometry (no meshes).
    #[error("empty scene")]
    EmptyScene,

    /// Caller passed a bad parameter (negative size, zero fov, …).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Offending parameter.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Emission for the chosen engine isn't supported in v1.
    #[error("engine `{engine}` emission is not implemented yet: {reason}")]
    EngineNotImplemented {
        /// Engine name.
        engine: &'static str,
        /// Why it's deferred.
        reason: String,
    },

    /// The external renderer executable was not found on `PATH`.
    ///
    /// Raised by the subprocess adapters
    /// ([`crate::subprocess::run_cycles`] /
    /// [`crate::subprocess::run_luxcore`]) when the renderer is not
    /// installed. Distinct from [`RenderError::Io`] (the process *was*
    /// found but the spawn failed) — a UI surfaces this as "install
    /// Cycles / LuxCoreRender", not "an I/O error occurred".
    #[error("renderer `{tool}` was not found on PATH — install it to render with this engine")]
    ToolNotAvailable {
        /// The executable name that was searched for.
        tool: &'static str,
    },

    /// The external renderer was found and launched but exited with a
    /// non-zero status.
    #[error("renderer `{tool}` exited with a failure status: {detail}")]
    RendererFailed {
        /// The executable name.
        tool: &'static str,
        /// Exit-status detail.
        detail: String,
    },

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON error.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Not implemented in v1.
    NotImplemented,
    /// Transient / IO / parse.
    Runtime,
}

impl RenderError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            RenderError::EmptyScene => "render.empty_scene",
            RenderError::BadParameter { .. } => "render.bad_parameter",
            RenderError::EngineNotImplemented { .. } => "render.engine_not_implemented",
            RenderError::ToolNotAvailable { .. } => "render.tool_not_available",
            RenderError::RendererFailed { .. } => "render.renderer_failed",
            RenderError::Io(_) => "render.io",
            RenderError::Ron(_) => "render.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RenderError::EmptyScene => ErrorCategory::Input,
            RenderError::BadParameter { .. } => ErrorCategory::Config,
            RenderError::EngineNotImplemented { .. }
            | RenderError::ToolNotAvailable { .. } => ErrorCategory::NotImplemented,
            RenderError::RendererFailed { .. }
            | RenderError::Io(_)
            | RenderError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One value of every `RenderError` variant.
    fn one_of_each() -> Vec<RenderError> {
        vec![
            RenderError::EmptyScene,
            RenderError::BadParameter {
                name: "fov",
                reason: "must be positive".into(),
            },
            RenderError::EngineNotImplemented {
                engine: "Octane",
                reason: "no emitter yet".into(),
            },
            RenderError::ToolNotAvailable { tool: "cycles" },
            RenderError::RendererFailed {
                tool: "luxcoreui",
                detail: "exit 1".into(),
            },
            RenderError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "missing",
            )),
            RenderError::Ron("bad syntax".into()),
        ]
    }

    #[test]
    fn every_variant_has_a_stable_code() {
        // Drives all seven arms of `code()`; codes are unique +
        // namespaced.
        let codes: Vec<&str> = one_of_each().iter().map(|e| e.code()).collect();
        assert_eq!(
            codes,
            [
                "render.empty_scene",
                "render.bad_parameter",
                "render.engine_not_implemented",
                "render.tool_not_available",
                "render.renderer_failed",
                "render.io",
                "render.ron",
            ]
        );
        for c in &codes {
            assert!(c.starts_with("render."), "code `{c}` must be namespaced");
        }
    }

    #[test]
    fn category_classifies_every_variant() {
        // Drives all arms of `category()`.
        let cats: Vec<ErrorCategory> =
            one_of_each().iter().map(|e| e.category()).collect();
        assert_eq!(cats[0], ErrorCategory::Input); // EmptyScene
        assert_eq!(cats[1], ErrorCategory::Config); // BadParameter
        assert_eq!(cats[2], ErrorCategory::NotImplemented); // EngineNotImpl
        assert_eq!(cats[3], ErrorCategory::NotImplemented); // ToolNotAvailable
        assert_eq!(cats[4], ErrorCategory::Runtime); // RendererFailed
        assert_eq!(cats[5], ErrorCategory::Runtime); // Io
        assert_eq!(cats[6], ErrorCategory::Runtime); // Ron
    }

    #[test]
    fn display_messages_are_non_empty_and_informative() {
        for e in one_of_each() {
            let msg = e.to_string();
            assert!(!msg.is_empty(), "every error needs a Display message");
        }
        // The ToolNotAvailable message names the missing tool — a UI
        // surfaces it as an "install X" hint.
        let tna = RenderError::ToolNotAvailable { tool: "cycles" };
        assert!(tna.to_string().contains("cycles"));
    }

    #[test]
    fn io_error_converts_via_from() {
        // The `#[from] std::io::Error` conversion.
        let io = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "no");
        let err: RenderError = io.into();
        assert_eq!(err.code(), "render.io");
        assert_eq!(err.category(), ErrorCategory::Runtime);
    }
}
