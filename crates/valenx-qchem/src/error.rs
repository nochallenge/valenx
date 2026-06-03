//! Error taxonomy for `valenx-qchem`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, QchemError>`]. The variants are intentionally coarse — a
//! quantum-chemistry caller usually only cares about five things:
//!
//! 1. Did a notation string / file fail to parse
//!    ([`QchemError::Parse`])?
//! 2. Is the molecular input itself nonsensical — a non-integer
//!    electron count, a charge / multiplicity pair with no valid
//!    determinant, two atoms on top of each other
//!    ([`QchemError::InvalidInput`])?
//! 3. Was a requested basis set unavailable for one of the elements in
//!    the molecule ([`QchemError::BasisNotFound`])?
//! 4. Did the SCF iteration fail to converge
//!    ([`QchemError::ScfNotConverged`])?
//! 5. Is this a documented stub awaiting a deeper subsystem
//!    ([`QchemError::NotYetImplemented`])?
//!
//! Use [`QchemError::code`] for stable log / telemetry tagging and
//! [`QchemError::category`] to bucket failures into Parse / Input /
//! Convergence / Capability without matching every variant. The pattern
//! mirrors `valenx-cheminf`'s `CheminfError`.

use std::fmt;

/// Errors produced by `valenx-qchem`.
#[derive(Debug, Clone, PartialEq)]
pub enum QchemError {
    /// A notation string or file failed to parse. `format` names the
    /// expected notation (`"xyz"`, `"basis"`); `detail` is a
    /// human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Notation being parsed (`"xyz"`, `"basis"`, …).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// The molecular specification is internally inconsistent: a
    /// charge / multiplicity pair that implies a negative or fractional
    /// electron count, an unphysical multiplicity, atoms placed at the
    /// same point, an empty molecule. A property of the *input*, not of
    /// a parse or a numerical failure.
    InvalidInput {
        /// Human-readable reason.
        reason: String,
    },

    /// A built-in basis set has no shell definition for one of the
    /// elements present in the molecule. `basis` is the requested set
    /// name; `element` is the offending atomic symbol.
    BasisNotFound {
        /// Requested basis-set name (e.g. `"sto-3g"`).
        basis: &'static str,
        /// Atomic symbol with no definition (e.g. `"Fe"`).
        element: String,
    },

    /// The self-consistent-field iteration ran out of cycles without
    /// satisfying the energy / density convergence thresholds.
    ScfNotConverged {
        /// Number of cycles attempted.
        iterations: usize,
        /// Final energy change between the last two cycles (Hartree).
        last_delta_energy: f64,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which subsystem
    /// the caller asked for; `detail` explains why it is a separate
    /// piece of work.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"dft"`).
        feature: &'static str,
        /// Human-readable explanation of the missing subsystem.
        detail: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 5+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A notation string / file failed to parse.
    Parse,
    /// User-supplied input is wrong (bad geometry or missing basis).
    Input,
    /// A numerical procedure failed to converge.
    Convergence,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl QchemError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"qchem.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            QchemError::Parse { .. } => "qchem.parse",
            QchemError::InvalidInput { .. } => "qchem.invalid_input",
            QchemError::BasisNotFound { .. } => "qchem.basis_not_found",
            QchemError::ScfNotConverged { .. } => "qchem.scf_not_converged",
            QchemError::NotYetImplemented { .. } => "qchem.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            QchemError::Parse { .. } => "parse",
            QchemError::InvalidInput { .. } | QchemError::BasisNotFound { .. } => "input",
            QchemError::ScfNotConverged { .. } => "convergence",
            QchemError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            QchemError::Parse { .. } => ErrorCategory::Parse,
            QchemError::InvalidInput { .. } | QchemError::BasisNotFound { .. } => {
                ErrorCategory::Input
            }
            QchemError::ScfNotConverged { .. } => ErrorCategory::Convergence,
            QchemError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`QchemError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        QchemError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`QchemError::InvalidInput`].
    pub fn invalid(reason: impl Into<String>) -> Self {
        QchemError::InvalidInput {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`QchemError::BasisNotFound`].
    pub fn basis_not_found(basis: &'static str, element: impl Into<String>) -> Self {
        QchemError::BasisNotFound {
            basis,
            element: element.into(),
        }
    }

    /// Convenience constructor for [`QchemError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str, detail: &'static str) -> Self {
        QchemError::NotYetImplemented { feature, detail }
    }
}

impl fmt::Display for QchemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QchemError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            QchemError::InvalidInput { reason } => {
                write!(f, "invalid input: {reason}")
            }
            QchemError::BasisNotFound { basis, element } => {
                write!(
                    f,
                    "basis set `{basis}` has no definition for element `{element}`"
                )
            }
            QchemError::ScfNotConverged {
                iterations,
                last_delta_energy,
            } => {
                write!(
                    f,
                    "SCF did not converge in {iterations} iterations \
                     (last |dE| = {last_delta_energy:.3e} Hartree)"
                )
            }
            QchemError::NotYetImplemented { feature, detail } => {
                write!(
                    f,
                    "qchem feature `{feature}` is not yet implemented: {detail}"
                )
            }
        }
    }
}

impl std::error::Error for QchemError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, QchemError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = QchemError::parse("xyz", "bad atom line");
        assert_eq!(err.code(), "qchem.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = QchemError::invalid("fractional electron count");
        assert_eq!(err.code(), "qchem.invalid_input");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = QchemError::basis_not_found("sto-3g", "Fe");
        assert_eq!(err.code(), "qchem.basis_not_found");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = QchemError::ScfNotConverged {
            iterations: 128,
            last_delta_energy: 1.0e-3,
        };
        assert_eq!(err.code(), "qchem.scf_not_converged");
        assert_eq!(err.category(), "convergence");
        assert_eq!(err.category_enum(), ErrorCategory::Convergence);

        let err = QchemError::not_yet("dft", "needs an XC-functional library");
        assert_eq!(err.code(), "qchem.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = QchemError::parse("xyz", "missing element").to_string();
        assert!(msg.contains("xyz"), "got: {msg}");
        assert!(msg.contains("missing element"), "got: {msg}");

        let msg = QchemError::basis_not_found("6-31g", "Br").to_string();
        assert!(msg.contains("6-31g"), "got: {msg}");
        assert!(msg.contains("Br"), "got: {msg}");

        let msg = QchemError::not_yet("geomopt", "needs analytic gradients").to_string();
        assert!(msg.contains("geomopt"), "got: {msg}");
        assert!(msg.contains("analytic gradients"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(QchemError::invalid("empty molecule"));
        assert!(err.to_string().contains("empty"));
    }
}
