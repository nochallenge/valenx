//! Error taxonomy for `valenx-structpredict`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, StructPredictError>`]. The variants are intentionally
//! coarse — a structure-prediction caller usually only cares about
//! five things:
//!
//! 1. Did the caller pass nonsense arguments — an empty sequence, a
//!    backbone with no Cα atoms, a negative grid size
//!    ([`StructPredictError::Invalid`])?
//! 2. Did a search fail to find a usable result — no template scored
//!    above threshold, no fragment matched a window, no particle
//!    cleared the picker ([`StructPredictError::NotFound`])?
//! 3. Did an iterative optimiser run out of steps without reaching its
//!    tolerance ([`StructPredictError::NotConverged`])?
//! 4. Did a file / data buffer fail to parse — a malformed MRC header,
//!    a truncated particle stack ([`StructPredictError::Parse`])?
//! 5. Is this a documented stub awaiting deeper work
//!    ([`StructPredictError::NotYetImplemented`])?
//!
//! Use [`StructPredictError::code`] for stable log / telemetry tagging
//! and [`StructPredictError::category`] to bucket failures into
//! Input / Search / Convergence / Parse / Capability without matching
//! every variant. The pattern mirrors `valenx-cheminf`'s
//! `CheminfError` and `valenx-biostruct`'s `BiostructError`.

use std::fmt;

/// Errors produced by `valenx-structpredict`.
#[derive(Debug, Clone, PartialEq)]
pub enum StructPredictError {
    /// Caller passed an argument the algorithm cannot accept: an empty
    /// sequence, a backbone missing its Cα atoms, a non-positive count
    /// or grid size, a mismatched array length, etc. A property of the
    /// *call*.
    Invalid {
        /// Logical parameter name (e.g. `"sequence"`, `"box_size"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// A search step found nothing usable — no template scored above
    /// the acceptance threshold, no backbone fragment matched a
    /// sequence window, no particle cleared the picker. Distinct from
    /// [`Invalid`](Self::Invalid): the inputs were well-formed, the
    /// search simply had no hit.
    NotFound {
        /// What was being searched for (`"template"`, `"fragment"`,
        /// `"particle"`, `"rotamer"`).
        what: &'static str,
        /// Human-readable context.
        detail: String,
    },

    /// An iterative optimiser exhausted its step budget without
    /// reaching the requested tolerance. `iterations` is how many
    /// steps ran; `detail` names the optimiser.
    NotConverged {
        /// The optimiser / protocol that did not converge.
        detail: String,
        /// Iterations performed before giving up.
        iterations: usize,
    },

    /// A file or data buffer failed to parse. `format` names the
    /// expected format (`"mrc"`, `"star"`); `detail` is a
    /// human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Format being parsed (`"mrc"`, `"star"`, …).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"de_novo_helix_bundle"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 5+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (bad argument or malformed model).
    Input,
    /// A search produced no usable result.
    Search,
    /// An iterative optimiser did not converge.
    Convergence,
    /// A file / data buffer failed to parse.
    Parse,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl StructPredictError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"structpredict.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            StructPredictError::Invalid { .. } => "structpredict.invalid",
            StructPredictError::NotFound { .. } => "structpredict.not_found",
            StructPredictError::NotConverged { .. } => "structpredict.not_converged",
            StructPredictError::Parse { .. } => "structpredict.parse",
            StructPredictError::NotYetImplemented { .. } => "structpredict.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            StructPredictError::Invalid { .. } => "input",
            StructPredictError::NotFound { .. } => "search",
            StructPredictError::NotConverged { .. } => "convergence",
            StructPredictError::Parse { .. } => "parse",
            StructPredictError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            StructPredictError::Invalid { .. } => ErrorCategory::Input,
            StructPredictError::NotFound { .. } => ErrorCategory::Search,
            StructPredictError::NotConverged { .. } => ErrorCategory::Convergence,
            StructPredictError::Parse { .. } => ErrorCategory::Parse,
            StructPredictError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`StructPredictError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        StructPredictError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`StructPredictError::NotFound`].
    pub fn not_found(what: &'static str, detail: impl Into<String>) -> Self {
        StructPredictError::NotFound {
            what,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`StructPredictError::NotConverged`].
    pub fn not_converged(detail: impl Into<String>, iterations: usize) -> Self {
        StructPredictError::NotConverged {
            detail: detail.into(),
            iterations,
        }
    }

    /// Convenience constructor for [`StructPredictError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        StructPredictError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`StructPredictError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        StructPredictError::NotYetImplemented { feature }
    }
}

impl fmt::Display for StructPredictError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructPredictError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            StructPredictError::NotFound { what, detail } => {
                write!(f, "no {what} found: {detail}")
            }
            StructPredictError::NotConverged { detail, iterations } => {
                write!(
                    f,
                    "{detail} did not converge after {iterations} iterations"
                )
            }
            StructPredictError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            StructPredictError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "structpredict feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for StructPredictError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, StructPredictError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = StructPredictError::invalid("sequence", "empty");
        assert_eq!(err.code(), "structpredict.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = StructPredictError::not_found("template", "no hit above 30%");
        assert_eq!(err.code(), "structpredict.not_found");
        assert_eq!(err.category(), "search");
        assert_eq!(err.category_enum(), ErrorCategory::Search);

        let err = StructPredictError::not_converged("CCD loop closure", 200);
        assert_eq!(err.code(), "structpredict.not_converged");
        assert_eq!(err.category(), "convergence");
        assert_eq!(err.category_enum(), ErrorCategory::Convergence);

        let err = StructPredictError::parse("mrc", "bad magic");
        assert_eq!(err.code(), "structpredict.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = StructPredictError::not_yet("de_novo");
        assert_eq!(err.code(), "structpredict.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = StructPredictError::not_converged("annealer", 5000).to_string();
        assert!(msg.contains("annealer"), "got: {msg}");
        assert!(msg.contains("5000"), "got: {msg}");

        let msg = StructPredictError::parse("star", "missing column").to_string();
        assert!(msg.contains("star"), "got: {msg}");
        assert!(msg.contains("missing column"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(StructPredictError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
