//! Error type for the path tracer.
//!
//! The renderer itself is total — `render` always produces a
//! framebuffer — so most of the surface area never returns a `Result`.
//! This type covers the few genuine *input* failure modes (a malformed
//! HDR environment file, a nonsensical render size) so a caller wiring
//! the tracer into a UI or a job runner has a typed error to branch on.

use thiserror::Error;

/// Errors raised when constructing a path-trace job from external
/// input.
#[derive(Debug, Error)]
pub enum PathTraceError {
    /// The supplied HDR environment file could not be decoded.
    ///
    /// Wraps the message from
    /// [`valenx_render_bridge::environment::EnvironmentMap::from_radiance_hdr`]
    /// — a bad Radiance magic line, an unsupported orientation, a
    /// truncated pixel stream, and so on.
    #[error("failed to decode the HDR environment: {0}")]
    BadEnvironment(String),

    /// A render parameter was out of range — a zero image dimension, a
    /// zero sample count, a non-positive field of view.
    #[error("invalid render parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// The scene has no geometry *and* no environment light, so a
    /// render would produce a uniformly black image. Surfaced as an
    /// error so a caller does not silently ship an empty frame.
    #[error("scene has no geometry and no environment light — nothing to render")]
    EmptyScene,
}

impl PathTraceError {
    /// A stable, kebab-cased identifier for the error — never changes
    /// across versions, so a UI can key off it.
    pub fn code(&self) -> &'static str {
        match self {
            PathTraceError::BadEnvironment(_) => "pathtrace.bad_environment",
            PathTraceError::BadParameter { .. } => "pathtrace.bad_parameter",
            PathTraceError::EmptyScene => "pathtrace.empty_scene",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_a_stable_code() {
        assert_eq!(
            PathTraceError::BadEnvironment("x".into()).code(),
            "pathtrace.bad_environment"
        );
        assert_eq!(
            PathTraceError::BadParameter {
                name: "width",
                reason: "zero".into()
            }
            .code(),
            "pathtrace.bad_parameter"
        );
        assert_eq!(PathTraceError::EmptyScene.code(), "pathtrace.empty_scene");
    }

    #[test]
    fn display_carries_context() {
        let e = PathTraceError::BadParameter {
            name: "samples_per_pixel",
            reason: "must be positive".into(),
        };
        let msg = e.to_string();
        assert!(msg.contains("samples_per_pixel"));
        assert!(msg.contains("must be positive"));
    }
}
