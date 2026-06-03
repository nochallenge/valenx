//! Error taxonomy for the 3-D external-aerodynamics solver.
//!
//! The steady driver itself is *total* — [`crate::run_windtunnel`]
//! always returns an [`crate::AeroResult`] with a `converged` flag
//! rather than an `Err` on an iteration-cap miss, because a
//! non-converged field is still useful (it can be inspected, or the
//! run resumed with a larger iteration budget). This type covers the
//! genuine *input* failure modes a caller building a wind-tunnel job
//! from external data must handle.
//!
//! Same shape as the sibling crates' error types: a stable
//! kebab-cased [`AeroError::code`] an LLM / UI can branch on, plus a
//! coarse [`ErrorCategory`] for escalation routing.

use thiserror::Error;

/// Errors raised when setting up or running a virtual-wind-tunnel case.
#[derive(Debug, Error)]
pub enum AeroError {
    /// A grid, fluid or wind parameter was out of range — a zero cell
    /// count, a non-positive domain length, a non-physical density /
    /// viscosity, a negative free-stream speed.
    #[error("invalid aero parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// The supplied body geometry cannot be voxelized into the tunnel
    /// — an empty triangle list, degenerate triangles only, or a body
    /// whose bounding box has zero extent.
    #[error("invalid body geometry: {0}")]
    BadGeometry(String),

    /// The body does not fit inside the requested domain, or the
    /// auto-sized domain came out degenerate — for example a body
    /// larger than a user-pinned tunnel box.
    #[error("body does not fit the wind-tunnel domain: {0}")]
    DomainTooSmall(String),

    /// The wind specification cannot drive a well-posed external flow
    /// — a zero-magnitude free-stream direction vector, or a
    /// turbulence intensity outside `[0, 1]`.
    #[error("ill-posed wind specification: {0}")]
    IllPosedWind(String),

    /// A requested feature is a documented stub, not yet a real
    /// implementation in this v1.
    #[error("feature not yet implemented: {0}")]
    NotYetImplemented(&'static str),
}

/// Coarse classification an LLM / UI can branch on to decide who to
/// escalate the error to. Mirrors `valenx_arch::ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied bad data — fix the geometry / wind vector.
    Input,
    /// A user-tunable knob is out of range — fix the dimension / count.
    Config,
    /// A transient or environmental failure — retry may help.
    Runtime,
    /// A bug in `valenx-aero` — file a report.
    Internal,
}

impl AeroError {
    /// A stable, kebab-cased identifier for the error. Never changes
    /// across versions, so an LLM agent can branch on it safely.
    pub fn code(&self) -> &'static str {
        match self {
            AeroError::BadParameter { .. } => "aero.bad_parameter",
            AeroError::BadGeometry(_) => "aero.bad_geometry",
            AeroError::DomainTooSmall(_) => "aero.domain_too_small",
            AeroError::IllPosedWind(_) => "aero.ill_posed_wind",
            AeroError::NotYetImplemented(_) => "aero.not_yet_implemented",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            AeroError::BadGeometry(_) => ErrorCategory::Input,
            AeroError::IllPosedWind(_) => ErrorCategory::Input,
            AeroError::BadParameter { .. } => ErrorCategory::Config,
            AeroError::DomainTooSmall(_) => ErrorCategory::Config,
            AeroError::NotYetImplemented(_) => ErrorCategory::Internal,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(AeroError, &'static str, ErrorCategory)> = vec![
            (
                AeroError::BadParameter {
                    name: "density",
                    reason: "must be positive".into(),
                },
                "aero.bad_parameter",
                ErrorCategory::Config,
            ),
            (
                AeroError::BadGeometry("empty mesh".into()),
                "aero.bad_geometry",
                ErrorCategory::Input,
            ),
            (
                AeroError::DomainTooSmall("body 5 m, box 2 m".into()),
                "aero.domain_too_small",
                ErrorCategory::Config,
            ),
            (
                AeroError::IllPosedWind("zero direction".into()),
                "aero.ill_posed_wind",
                ErrorCategory::Input,
            ),
            (
                AeroError::NotYetImplemented("LES"),
                "aero.not_yet_implemented",
                ErrorCategory::Internal,
            ),
        ];
        for (err, code, cat) in cases {
            assert_eq!(err.code(), code, "wrong code for {err:?}");
            assert_eq!(err.category(), cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn display_carries_context() {
        let e = AeroError::BadParameter {
            name: "free_stream_speed",
            reason: "must be non-negative".into(),
        };
        let s = e.to_string();
        assert!(s.contains("free_stream_speed"));
        assert!(s.contains("non-negative"));
    }
}
