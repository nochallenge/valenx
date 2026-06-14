//! Error taxonomy for `valenx-bonemech`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, BoneError>`](crate::Result). The variants are deliberately
//! coarse — a structural-mechanics caller usually only needs to know
//! whether an input was out of its physical domain (a non-positive area,
//! a negative modulus, a hollow cross-section whose inner diameter is not
//! smaller than its outer diameter) or whether two inputs were
//! geometrically inconsistent.
//!
//! The enum derives [`thiserror::Error`]; each variant carries a
//! human-readable `Display` message via its `#[error(...)]` attribute,
//! mirroring `valenx-springs`'s `SpringsError`. [`BoneError::code`]
//! returns a stable kebab-cased identifier for log / telemetry tagging
//! and [`BoneError::category`] buckets failures without matching every
//! variant.

use thiserror::Error;

/// Errors produced by `valenx-bonemech`.
///
/// Marked `#[non_exhaustive]`: more failure modes may be added as the
/// crate grows, so downstream matches must include a wildcard arm.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum BoneError {
    /// A scalar parameter fell outside its physical domain: a
    /// non-positive cross-sectional area, a non-positive or non-finite
    /// elastic modulus, a negative stress, a negative apparent density,
    /// and so on. `name` is the logical parameter name; `reason` is a
    /// human-readable explanation surfaced verbatim in the UI.
    #[error("invalid `{name}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"area_mm2"`, `"elastic_modulus_gpa"`).
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// Two inputs describing one geometry disagree: the inner diameter of
    /// a hollow circular cross-section is not strictly smaller than the
    /// outer diameter, the extreme-fibre distance `c` exceeds the section
    /// radius, etc. A property of the *pair* of inputs, not of either one
    /// alone.
    #[error("inconsistent geometry: {0}")]
    Geometry(String),
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on the full variant set.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A single scalar input was out of its physical domain.
    Input,
    /// Two inputs were geometrically inconsistent.
    Geometry,
}

impl BoneError {
    /// Convenience constructor for [`BoneError::Invalid`].
    pub fn invalid(name: &'static str, reason: impl Into<String>) -> Self {
        BoneError::Invalid {
            name,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`BoneError::Geometry`].
    pub fn geometry(reason: impl Into<String>) -> Self {
        BoneError::Geometry(reason.into())
    }

    /// Stable kebab-cased error code suitable for log / telemetry
    /// tagging. Format: `"bonemech.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            BoneError::Invalid { .. } => "bonemech.invalid",
            BoneError::Geometry(_) => "bonemech.geometry",
        }
    }

    /// Coarse [`ErrorCategory`] for callers that want to `match` on a
    /// small stable enum instead of the full variant set.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BoneError::Invalid { .. } => ErrorCategory::Input,
            BoneError::Geometry(_) => ErrorCategory::Geometry,
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, BoneError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let e = BoneError::invalid("area_mm2", "must be positive");
        assert_eq!(e.code(), "bonemech.invalid");
        assert_eq!(e.category(), ErrorCategory::Input);

        let e = BoneError::geometry("inner diameter >= outer diameter");
        assert_eq!(e.code(), "bonemech.geometry");
        assert_eq!(e.category(), ErrorCategory::Geometry);
    }

    #[test]
    fn display_is_informative() {
        let msg = BoneError::invalid("elastic_modulus_gpa", "must be finite").to_string();
        assert!(msg.contains("elastic_modulus_gpa"), "got: {msg}");
        assert!(msg.contains("must be finite"), "got: {msg}");

        let msg = BoneError::geometry("c exceeds outer radius").to_string();
        assert!(msg.contains("c exceeds outer radius"), "got: {msg}");
    }

    #[test]
    fn error_is_a_std_error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(BoneError::invalid("force_n", "negative"));
        assert!(err.to_string().contains("force_n"));
    }
}
