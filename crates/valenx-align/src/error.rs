//! Error taxonomy for `valenx-align`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, AlignError>`]. The variants are intentionally coarse —
//! an alignment-toolkit caller usually only cares about four things:
//!
//! 1. Did a file / format fail to parse ([`AlignError::Parse`])?
//! 2. Did the caller pass a nonsense argument — empty input, an
//!    out-of-range index, a non-positive `k` ([`AlignError::Invalid`])?
//! 3. Did two inputs disagree in shape — mismatched profile widths,
//!    a query that does not fit a band, unequal alignment-row lengths
//!    ([`AlignError::Dimension`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`AlignError::NotYetImplemented`])?
//!
//! Use [`AlignError::code`] for stable log / telemetry tagging and
//! [`AlignError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError` and `valenx-occt-*`'s error enums.

use std::fmt;

/// Errors produced by `valenx-align`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlignError {
    /// A file or in-memory string failed to parse. `format` names the
    /// expected format (`"clustal"`, `"stockholm"`, …); `detail` is a
    /// human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Format being parsed (e.g. `"clustal"`, `"phylip"`, `"sam"`).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive length, an out-of-range coordinate, etc.
    /// A property of the *call*, not of a file being parsed.
    Invalid {
        /// Logical parameter name (e.g. `"k"`, `"band"`, `"window"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Two inputs disagree in shape — e.g. profile columns of unequal
    /// width, an aligned MSA whose rows have different lengths, or a
    /// k-band too narrow to admit the diagonal a global alignment
    /// needs. Distinct from [`AlignError::Invalid`] because the
    /// individual inputs are each well-formed; only their combination
    /// is wrong.
    Dimension {
        /// Short description of the mismatch.
        detail: String,
    },

    /// A dynamic-programming matrix would be too large to allocate
    /// safely: the requested cell count exceeds
    /// [`MAX_DP_CELLS`](crate::limits::MAX_DP_CELLS), or the product
    /// `(n+1)·(m+1)` overflowed `usize`. Guards against an
    /// out-of-memory / denial-of-service from two very long sequences
    /// (a full O(n·m) matrix of, say, two 50 kb sequences is ~10 GB).
    /// The caller should use a linear-space routine
    /// ([`hirschberg`](crate::pairwise::hirschberg::hirschberg)) or a
    /// banded one instead.
    TooLarge {
        /// Number of DP cells requested (`usize::MAX` reported when the
        /// dimension product itself overflowed).
        cells: usize,
        /// The configured ceiling that was exceeded.
        max: usize,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"sw_traceback_gpu"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 4+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A file / string failed to parse.
    Parse,
    /// User-supplied input is wrong (bad argument or shape mismatch).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl AlignError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"align.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            AlignError::Parse { .. } => "align.parse",
            AlignError::Invalid { .. } => "align.invalid",
            AlignError::Dimension { .. } => "align.dimension",
            AlignError::TooLarge { .. } => "align.too_large",
            AlignError::NotYetImplemented { .. } => "align.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            AlignError::Parse { .. } => "parse",
            AlignError::Invalid { .. }
            | AlignError::Dimension { .. }
            | AlignError::TooLarge { .. } => "input",
            AlignError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            AlignError::Parse { .. } => ErrorCategory::Parse,
            AlignError::Invalid { .. }
            | AlignError::Dimension { .. }
            | AlignError::TooLarge { .. } => ErrorCategory::Input,
            AlignError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`AlignError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        AlignError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`AlignError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        AlignError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`AlignError::Dimension`].
    pub fn dimension(detail: impl Into<String>) -> Self {
        AlignError::Dimension {
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`AlignError::TooLarge`].
    pub fn too_large(cells: usize, max: usize) -> Self {
        AlignError::TooLarge { cells, max }
    }

    /// Convenience constructor for [`AlignError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        AlignError::NotYetImplemented { feature }
    }
}

impl fmt::Display for AlignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            AlignError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            AlignError::Dimension { detail } => {
                write!(f, "dimension mismatch: {detail}")
            }
            AlignError::TooLarge { cells, max } => {
                write!(
                    f,
                    "dynamic-programming matrix too large: {cells} cells \
                     exceeds the {max}-cell limit; use a linear-space \
                     (Hirschberg) or banded routine"
                )
            }
            AlignError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "align feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for AlignError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, AlignError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = AlignError::parse("clustal", "bad header");
        assert_eq!(err.code(), "align.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = AlignError::invalid("k", "must be positive");
        assert_eq!(err.code(), "align.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = AlignError::dimension("rows differ");
        assert_eq!(err.code(), "align.dimension");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = AlignError::too_large(1_000, 100);
        assert_eq!(err.code(), "align.too_large");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = AlignError::not_yet("gpu_sw");
        assert_eq!(err.code(), "align.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn too_large_display_mentions_counts() {
        let msg = AlignError::too_large(5_000, 100).to_string();
        assert!(msg.contains("5000"), "got: {msg}");
        assert!(msg.contains("100"), "got: {msg}");
    }

    #[test]
    fn display_is_informative() {
        let msg = AlignError::dimension("widths 3 vs 4").to_string();
        assert!(msg.contains('3'), "got: {msg}");
        assert!(msg.contains('4'), "got: {msg}");

        let msg = AlignError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(AlignError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
