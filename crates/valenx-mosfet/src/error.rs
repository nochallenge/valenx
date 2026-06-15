//! Error taxonomy for `valenx-mosfet`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, MosfetError>`]. The variants are intentionally coarse —
//! a device-model caller usually only cares about two things:
//!
//! 1. Did the caller pass a nonsensical parameter — a negative or zero
//!    transconductance parameter, a non-finite voltage
//!    ([`MosfetError::Invalid`])?
//! 2. Was a parameter outside the physical domain the square-law model
//!    accepts — e.g. a NaN supplied where a real bias is required
//!    ([`MosfetError::Domain`])?
//!
//! Use [`MosfetError::code`] for stable log / telemetry tagging and
//! [`MosfetError::category`] to bucket failures into Input / Domain
//! without matching every variant. The pattern mirrors the
//! `Invalid { what, reason }` shape used across the Valenx physics
//! crates (e.g. `valenx-popgen`'s `PopgenError`).

use thiserror::Error;

/// Errors produced by `valenx-mosfet`.
///
/// All constructors that build a [`crate::device::Mosfet`] or evaluate
/// its IV / transconductance equations funnel their rejections through
/// these two variants.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum MosfetError {
    /// Caller passed an argument the model cannot accept: a
    /// non-positive transconductance parameter `k`, a non-positive
    /// gate-oxide quantity, or any device parameter that must be
    /// strictly positive but was supplied as zero or negative. A
    /// property of the *call*, not of a file being parsed.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"k"`, `"vth"`, `"vgs"`).
        what: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A bias voltage or model parameter was outside the domain the
    /// square-law equations are defined on — most commonly a non-finite
    /// value (`NaN` / `±∞`) supplied where a real number is required.
    #[error("`{what}` out of domain: {reason}")]
    Domain {
        /// Logical parameter name (e.g. `"vgs"`, `"vds"`).
        what: &'static str,
        /// Human-readable reason the value is out of domain.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the full set of error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (bad argument value).
    Input,
    /// A value lies outside the model's mathematical domain.
    Domain,
}

impl MosfetError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"mosfet.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            MosfetError::Invalid { .. } => "mosfet.invalid",
            MosfetError::Domain { .. } => "mosfet.domain",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            MosfetError::Invalid { .. } => ErrorCategory::Input,
            MosfetError::Domain { .. } => ErrorCategory::Domain,
        }
    }

    /// Convenience constructor for [`MosfetError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        MosfetError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`MosfetError::Domain`].
    pub fn domain(what: &'static str, reason: impl Into<String>) -> Self {
        MosfetError::Domain {
            what,
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, MosfetError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = MosfetError::invalid("k", "must be > 0");
        assert_eq!(err.code(), "mosfet.invalid");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = MosfetError::domain("vgs", "must be finite");
        assert_eq!(err.code(), "mosfet.domain");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = MosfetError::invalid("k", "must be > 0").to_string();
        assert!(msg.contains('k'), "got: {msg}");
        assert!(msg.contains("must be > 0"), "got: {msg}");

        let msg = MosfetError::domain("vds", "got NaN").to_string();
        assert!(msg.contains("vds"), "got: {msg}");
        assert!(msg.contains("NaN"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(MosfetError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
