//! Error taxonomy for `valenx-md`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, MdError>`]. The variants are intentionally coarse — a
//! molecular-dynamics caller usually only cares about five things:
//!
//! 1. Did a structure / trajectory / topology file fail to parse
//!    ([`MdError::Parse`])?
//! 2. Did the caller hand over a malformed system or topology
//!    ([`MdError::Invalid`]) — an atom referenced by a bond that does
//!    not exist, a negative mass, an empty system?
//! 3. Did two arrays / a system and a force field disagree on size
//!    ([`MdError::DimensionMismatch`])?
//! 4. Did an iterative routine — a minimiser, an SCF-style self-
//!    consistent constraint solve, a PME convergence loop — fail to
//!    converge ([`MdError::NotConverged`])?
//! 5. Is this a documented stub awaiting deeper work
//!    ([`MdError::NotYetImplemented`])?
//!
//! Use [`MdError::code`] for stable log / telemetry tagging and
//! [`MdError::category`] to bucket failures into Parse / Input /
//! Numerics / Capability without matching every variant. The pattern
//! mirrors `valenx-bioseq`'s `BioseqError` and the other Round 6
//! crates.

use std::fmt;

/// Errors produced by `valenx-md`.
#[derive(Debug, Clone, PartialEq)]
pub enum MdError {
    /// A structure / trajectory / topology file or string failed to
    /// parse. `format` names the expected format; `detail` is a
    /// human-readable reason surfaced verbatim.
    Parse {
        /// Format being parsed (e.g. `"pdb"`, `"xyz"`, `"gro"`, `"dcd"`).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// The caller supplied a system, topology or force field that is
    /// internally inconsistent — a bond referencing a nonexistent
    /// atom, a non-finite or non-positive mass, an empty system, a
    /// degenerate simulation box. A property of the *input*.
    Invalid {
        /// Logical parameter / object name (e.g. `"mass"`, `"bond"`,
        /// `"box"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Two quantities that must agree on a length / dimension do not —
    /// a force-field that types fewer atoms than the system has, a
    /// trajectory frame of the wrong width, mismatched position and
    /// velocity arrays.
    DimensionMismatch {
        /// Human-readable reason, naming both sides and their sizes.
        reason: String,
    },

    /// An iterative algorithm — energy minimisation, the SHAKE /
    /// RATTLE constraint solve, a PME parameter search — ran out of
    /// iterations before reaching its tolerance.
    NotConverged {
        /// Which algorithm failed to converge.
        algorithm: &'static str,
        /// Iterations performed before giving up.
        iterations: usize,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier.
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A structure / trajectory file failed to parse.
    Parse,
    /// User-supplied input is wrong (bad system, bad topology, size
    /// mismatch).
    Input,
    /// A numerical routine failed to converge.
    Numerics,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl MdError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"md.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            MdError::Parse { .. } => "md.parse",
            MdError::Invalid { .. } => "md.invalid",
            MdError::DimensionMismatch { .. } => "md.dimension_mismatch",
            MdError::NotConverged { .. } => "md.not_converged",
            MdError::NotYetImplemented { .. } => "md.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            MdError::Parse { .. } => "parse",
            MdError::Invalid { .. } | MdError::DimensionMismatch { .. } => "input",
            MdError::NotConverged { .. } => "numerics",
            MdError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            MdError::Parse { .. } => ErrorCategory::Parse,
            MdError::Invalid { .. } | MdError::DimensionMismatch { .. } => ErrorCategory::Input,
            MdError::NotConverged { .. } => ErrorCategory::Numerics,
            MdError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`MdError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        MdError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`MdError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        MdError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`MdError::DimensionMismatch`].
    pub fn dimension(reason: impl Into<String>) -> Self {
        MdError::DimensionMismatch {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`MdError::NotConverged`].
    pub fn not_converged(algorithm: &'static str, iterations: usize) -> Self {
        MdError::NotConverged {
            algorithm,
            iterations,
        }
    }

    /// Convenience constructor for [`MdError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        MdError::NotYetImplemented { feature }
    }
}

impl fmt::Display for MdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MdError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            MdError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            MdError::DimensionMismatch { reason } => {
                write!(f, "dimension mismatch: {reason}")
            }
            MdError::NotConverged {
                algorithm,
                iterations,
            } => {
                write!(
                    f,
                    "`{algorithm}` did not converge after {iterations} iterations"
                )
            }
            MdError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "md feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for MdError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, MdError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = MdError::parse("pdb", "bad ATOM record");
        assert_eq!(err.code(), "md.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = MdError::invalid("mass", "must be positive");
        assert_eq!(err.code(), "md.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = MdError::dimension("3 charges for 4 atoms");
        assert_eq!(err.code(), "md.dimension_mismatch");
        assert_eq!(err.category(), "input");

        let err = MdError::not_converged("steepest-descent", 500);
        assert_eq!(err.code(), "md.not_converged");
        assert_eq!(err.category(), "numerics");
        assert_eq!(err.category_enum(), ErrorCategory::Numerics);

        let err = MdError::not_yet("free-energy-perturbation");
        assert_eq!(err.code(), "md.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = MdError::invalid("box", "zero volume").to_string();
        assert!(msg.contains("zero volume"), "got: {msg}");

        let msg = MdError::not_converged("l-bfgs", 1000).to_string();
        assert!(msg.contains("1000"), "got: {msg}");

        let msg = MdError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(MdError::dimension("x"));
        assert!(err.to_string().contains('x'));
    }
}
