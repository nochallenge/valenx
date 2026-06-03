//! Error taxonomy for `valenx-sysbio`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, SysbioError>`]. The variants are intentionally coarse
//! — a systems-biology caller usually only cares about five things:
//!
//! 1. Did an SBML / SBOL document fail to parse
//!    ([`SysbioError::Parse`])?
//! 2. Is a constructed model itself structurally invalid — a reaction
//!    referencing an unknown species, a stoichiometry matrix with a
//!    dimension mismatch, a circuit gate with no assigned part
//!    ([`SysbioError::InvalidModel`])?
//! 3. Did the caller pass nonsense arguments — an empty species set, a
//!    negative time step, a parameter-scan range with `lo > hi`
//!    ([`SysbioError::Invalid`])?
//! 4. Did a numerical method fail to make progress — a Newton solve
//!    that did not converge, an LP that is infeasible or unbounded, an
//!    integrator that hit its step-count ceiling
//!    ([`SysbioError::NotConverged`])?
//! 5. Is this a documented stub awaiting deeper work
//!    ([`SysbioError::NotYetImplemented`])?
//!
//! Use [`SysbioError::code`] for stable log / telemetry tagging and
//! [`SysbioError::category`] to bucket failures into Parse / Input /
//! Numeric / Capability without matching every variant. The pattern
//! mirrors `valenx-genomics`'s `GenomicsError` and `valenx-cheminf`'s
//! `CheminfError`.

use std::fmt;

/// Errors produced by `valenx-sysbio`.
#[derive(Debug, Clone, PartialEq)]
pub enum SysbioError {
    /// An SBML / SBOL document or a kinetic-law expression failed to
    /// parse. `format` names the expected format (`"sbml"`, `"sbol"`,
    /// `"ratelaw"`); `detail` is a human-readable reason surfaced
    /// verbatim in the UI.
    Parse {
        /// Format being parsed (`"sbml"`, `"sbol"`, `"ratelaw"`, …).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// A constructed model is structurally inconsistent: a reaction
    /// referencing a species that is not in the species table, a
    /// stoichiometry matrix whose row count disagrees with the species
    /// list, a genetic circuit with a dangling wire. A property of the
    /// *model*, not of a parse position or a caller argument.
    InvalidModel {
        /// Model kind (`"reaction_network"`, `"circuit"`, `"sbol"`, …).
        kind: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive step size, a probability outside `[0, 1]`,
    /// a scan range with `lo > hi`, a percentile outside `[0, 100]`. A
    /// property of the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"dt"`, `"n_samples"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// A numerical method failed to make progress: a Newton iteration
    /// that exceeded its iteration budget without satisfying the
    /// tolerance, an LP that is infeasible or unbounded, an adaptive
    /// integrator that could not shrink its step below the floor.
    NotConverged {
        /// Method name (`"newton"`, `"simplex"`, `"rk45"`, `"bdf"`).
        method: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"sbml_l3_packages"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 5+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A document / expression failed to parse.
    Parse,
    /// User-supplied input is wrong (bad model or bad argument).
    Input,
    /// A numerical method did not converge / a problem was infeasible.
    Numeric,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl SysbioError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"sysbio.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            SysbioError::Parse { .. } => "sysbio.parse",
            SysbioError::InvalidModel { .. } => "sysbio.invalid_model",
            SysbioError::Invalid { .. } => "sysbio.invalid",
            SysbioError::NotConverged { .. } => "sysbio.not_converged",
            SysbioError::NotYetImplemented { .. } => "sysbio.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            SysbioError::Parse { .. } => "parse",
            SysbioError::InvalidModel { .. } | SysbioError::Invalid { .. } => "input",
            SysbioError::NotConverged { .. } => "numeric",
            SysbioError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            SysbioError::Parse { .. } => ErrorCategory::Parse,
            SysbioError::InvalidModel { .. } | SysbioError::Invalid { .. } => ErrorCategory::Input,
            SysbioError::NotConverged { .. } => ErrorCategory::Numeric,
            SysbioError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`SysbioError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        SysbioError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`SysbioError::InvalidModel`].
    pub fn invalid_model(kind: &'static str, reason: impl Into<String>) -> Self {
        SysbioError::InvalidModel {
            kind,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`SysbioError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        SysbioError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`SysbioError::NotConverged`].
    pub fn not_converged(method: &'static str, reason: impl Into<String>) -> Self {
        SysbioError::NotConverged {
            method,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`SysbioError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        SysbioError::NotYetImplemented { feature }
    }
}

impl fmt::Display for SysbioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SysbioError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            SysbioError::InvalidModel { kind, reason } => {
                write!(f, "invalid {kind} model: {reason}")
            }
            SysbioError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            SysbioError::NotConverged { method, reason } => {
                write!(f, "{method} did not converge: {reason}")
            }
            SysbioError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "sysbio feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for SysbioError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, SysbioError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = SysbioError::parse("sbml", "missing listOfSpecies");
        assert_eq!(err.code(), "sysbio.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = SysbioError::invalid_model("reaction_network", "unknown species");
        assert_eq!(err.code(), "sysbio.invalid_model");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = SysbioError::invalid("dt", "must be positive");
        assert_eq!(err.code(), "sysbio.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = SysbioError::not_converged("simplex", "infeasible");
        assert_eq!(err.code(), "sysbio.not_converged");
        assert_eq!(err.category(), "numeric");
        assert_eq!(err.category_enum(), ErrorCategory::Numeric);

        let err = SysbioError::not_yet("sbml_l3_packages");
        assert_eq!(err.code(), "sysbio.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = SysbioError::parse("sbml", "bad xml").to_string();
        assert!(msg.contains("sbml"), "got: {msg}");
        assert!(msg.contains("bad xml"), "got: {msg}");

        let msg = SysbioError::not_converged("newton", "stalled").to_string();
        assert!(msg.contains("newton"), "got: {msg}");
        assert!(msg.contains("stalled"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(SysbioError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
