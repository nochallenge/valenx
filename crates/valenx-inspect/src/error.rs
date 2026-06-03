//! Inspection workbench error taxonomy.
//!
//! Same shape as `valenx-draft::DraftError`: stable `code()` strings
//! an LLM can branch on plus a coarse [`ErrorCategory`] for UI /
//! escalation routing.

use thiserror::Error;

/// Errors raised by inspection-workbench operations.
#[derive(Debug, Error)]
pub enum InspectError {
    /// A measurement requires more inputs than were supplied (e.g.
    /// a polyline with fewer than 2 points, an empty mesh for volume).
    #[error("not enough geometry for `{kind}`: {reason}")]
    NotEnoughGeometry {
        /// Measurement kind label (e.g. `"PolylineLength"`).
        kind: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// User-supplied parameter was out of range (negative tolerance,
    /// degenerate triangle, …).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"upper_dev"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// Asked to verify a GD&T characteristic that v1 doesn't know how
    /// to evaluate from a single scalar (e.g. Cylindricity needs many
    /// surface points, not a single radius).
    #[error("GD&T characteristic `{characteristic}` is not implemented yet: {reason}")]
    GdtNotImplemented {
        /// Characteristic name (e.g. `"Cylindricity"`).
        characteristic: &'static str,
        /// Why it's deferred and what's needed.
        reason: String,
    },

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON parse / serialise error.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse category an LLM / UI can branch on to decide who to escalate
/// the error to.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad geometry (fix the picks).
    Input,
    /// User-tunable knob out of range.
    Config,
    /// Functionality not yet implemented in v1.
    NotImplemented,
    /// Transient or environmental failure — IO, parsing.
    Runtime,
}

impl InspectError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            InspectError::NotEnoughGeometry { .. } => "inspect.not_enough_geometry",
            InspectError::BadParameter { .. } => "inspect.bad_parameter",
            InspectError::GdtNotImplemented { .. } => "inspect.gdt_not_implemented",
            InspectError::Io(_) => "inspect.io",
            InspectError::Ron(_) => "inspect.ron",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            InspectError::NotEnoughGeometry { .. } => ErrorCategory::Input,
            InspectError::BadParameter { .. } => ErrorCategory::Config,
            InspectError::GdtNotImplemented { .. } => ErrorCategory::NotImplemented,
            InspectError::Io(_) | InspectError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(InspectError, &'static str, ErrorCategory)> = vec![
            (
                InspectError::NotEnoughGeometry {
                    kind: "PolylineLength",
                    reason: "needs 2 points".into(),
                },
                "inspect.not_enough_geometry",
                ErrorCategory::Input,
            ),
            (
                InspectError::BadParameter {
                    name: "upper_dev",
                    reason: "must be >= lower_dev".into(),
                },
                "inspect.bad_parameter",
                ErrorCategory::Config,
            ),
            (
                InspectError::GdtNotImplemented {
                    characteristic: "Cylindricity",
                    reason: "needs surface samples".into(),
                },
                "inspect.gdt_not_implemented",
                ErrorCategory::NotImplemented,
            ),
            (
                InspectError::Ron("bad token".into()),
                "inspect.ron",
                ErrorCategory::Runtime,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code);
            assert_eq!(err.category(), expected_cat);
        }
    }
}
