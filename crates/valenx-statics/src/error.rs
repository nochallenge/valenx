//! Error taxonomy for `valenx-statics`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, StaticsError>`]. The variants are intentionally coarse:
//! a statics caller usually only cares about three things.
//!
//! 1. Did the caller pass a nonsense argument — a non-finite
//!    coordinate, a negative or non-positive span
//!    ([`StaticsError::Invalid`])?
//! 2. Is the structure geometrically degenerate — the two supports of
//!    a beam placed at the same location, so the moment equation cannot
//!    be solved ([`StaticsError::Degenerate`])?
//! 3. Is the structure statically *indeterminate* or a *mechanism* —
//!    not solvable by the three planar equilibrium equations alone
//!    ([`StaticsError::Indeterminate`])?
//!
//! Use [`StaticsError::code`] for stable log / telemetry tagging and
//! [`StaticsError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-gears`'s `GearsError`.

use thiserror::Error;

/// Errors produced by `valenx-statics`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum StaticsError {
    /// Caller passed an argument the model cannot accept: a non-finite
    /// (`NaN` / infinite) coordinate or force component, a non-positive
    /// beam span, a support / load placed off the beam, and so on. A
    /// property of the *call*, not of the structure's topology.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"span"`, `"load.position"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The structure is geometrically degenerate: the two supports of a
    /// simply-supported beam coincide (zero lever arm between them), so
    /// the moment-balance equation has no unique solution.
    #[error("degenerate geometry: {0}")]
    Degenerate(String),

    /// The structure cannot be solved by the three planar equilibrium
    /// equations alone — it is statically indeterminate (too many
    /// reactions) or a mechanism (too few). The string explains which.
    #[error("statically indeterminate or a mechanism: {0}")]
    Indeterminate(String),
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the full error variant set.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A supplied argument is wrong (bad number, off-beam position).
    Input,
    /// The structure's geometry / topology is unsolvable
    /// (degenerate supports, indeterminate, or a mechanism).
    Geometry,
}

impl StaticsError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"statics.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            StaticsError::Invalid { .. } => "statics.invalid",
            StaticsError::Degenerate(_) => "statics.degenerate",
            StaticsError::Indeterminate(_) => "statics.indeterminate",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            StaticsError::Invalid { .. } => "input",
            StaticsError::Degenerate(_) | StaticsError::Indeterminate(_) => "geometry",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            StaticsError::Invalid { .. } => ErrorCategory::Input,
            StaticsError::Degenerate(_) | StaticsError::Indeterminate(_) => ErrorCategory::Geometry,
        }
    }

    /// Convenience constructor for [`StaticsError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        StaticsError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`StaticsError::Degenerate`].
    pub fn degenerate(reason: impl Into<String>) -> Self {
        StaticsError::Degenerate(reason.into())
    }

    /// Convenience constructor for [`StaticsError::Indeterminate`].
    pub fn indeterminate(reason: impl Into<String>) -> Self {
        StaticsError::Indeterminate(reason.into())
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, StaticsError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = StaticsError::invalid("span", "must be positive");
        assert_eq!(err.code(), "statics.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = StaticsError::degenerate("supports coincide");
        assert_eq!(err.code(), "statics.degenerate");
        assert_eq!(err.category(), "geometry");
        assert_eq!(err.category_enum(), ErrorCategory::Geometry);

        let err = StaticsError::indeterminate("3 reactions, 3 equations, 1 redundant");
        assert_eq!(err.code(), "statics.indeterminate");
        assert_eq!(err.category(), "geometry");
        assert_eq!(err.category_enum(), ErrorCategory::Geometry);
    }

    #[test]
    fn display_is_informative() {
        let msg = StaticsError::invalid("load.position", "off the beam").to_string();
        assert!(msg.contains("load.position"), "got: {msg}");
        assert!(msg.contains("off the beam"), "got: {msg}");

        let msg = StaticsError::degenerate("zero lever arm").to_string();
        assert!(msg.contains("zero lever arm"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(StaticsError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
