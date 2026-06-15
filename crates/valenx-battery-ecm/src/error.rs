//! Error taxonomy for `valenx-battery-ecm`.
//!
//! Every fallible public function returns [`Result<_, EcmError>`]. The
//! variants are deliberately coarse — a caller of an equivalent-circuit
//! model usually only needs to know whether a *parameter* was rejected
//! (a non-positive resistance, capacitance or capacity; a state-of-charge
//! outside `[0, 1]`), whether an *OCV-SoC table* was malformed (too few
//! points, out-of-order or non-monotonic breakpoints), or whether a
//! *simulation step* was nonsensical (a negative time step).
//!
//! Use [`EcmError::code`] for stable log / telemetry tagging and
//! [`EcmError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-gears`' `GearsError` and
//! `valenx-popgen`'s `PopgenError`.

use thiserror::Error;

/// Errors produced by `valenx-battery-ecm`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum EcmError {
    /// A scalar parameter is outside its physical domain — a
    /// non-positive resistance, capacitance, capacity or time constant,
    /// a non-finite value, etc. A property of the *call's argument*,
    /// not of a table being parsed.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"r0"`, `"capacity_ah"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// A state-of-charge value (or breakpoint) is outside the closed
    /// unit interval `[0, 1]`. SoC is a fraction of full charge, so any
    /// value below `0` or above `1` is unphysical.
    #[error("state-of-charge {soc} out of range [0, 1] ({context})")]
    SocOutOfRange {
        /// The offending value.
        soc: f64,
        /// Short context label (e.g. `"initial SoC"`, `"OCV breakpoint"`).
        context: &'static str,
    },

    /// An open-circuit-voltage / state-of-charge lookup table is
    /// malformed: fewer than two points, a SoC and OCV column of
    /// different lengths, SoC breakpoints that are not strictly
    /// increasing, or OCV values that are not monotonically
    /// non-decreasing with SoC (a real cell's rest voltage never falls
    /// as it charges).
    #[error("invalid OCV-SoC table: {reason}")]
    Table {
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A scalar argument is out of its physical domain.
    Parameter,
    /// An open-circuit-voltage lookup table is malformed.
    Table,
}

impl EcmError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"battery_ecm.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            EcmError::Invalid { .. } => "battery_ecm.invalid",
            EcmError::SocOutOfRange { .. } => "battery_ecm.soc_out_of_range",
            EcmError::Table { .. } => "battery_ecm.table",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            EcmError::Invalid { .. } | EcmError::SocOutOfRange { .. } => ErrorCategory::Parameter,
            EcmError::Table { .. } => ErrorCategory::Table,
        }
    }

    /// Convenience constructor for [`EcmError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        EcmError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`EcmError::Table`].
    pub fn table(reason: impl Into<String>) -> Self {
        EcmError::Table {
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, EcmError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = EcmError::invalid("r0", "must be positive");
        assert_eq!(err.code(), "battery_ecm.invalid");
        assert_eq!(err.category(), ErrorCategory::Parameter);

        let err = EcmError::SocOutOfRange {
            soc: 1.5,
            context: "initial SoC",
        };
        assert_eq!(err.code(), "battery_ecm.soc_out_of_range");
        assert_eq!(err.category(), ErrorCategory::Parameter);

        let err = EcmError::table("breakpoints not increasing");
        assert_eq!(err.code(), "battery_ecm.table");
        assert_eq!(err.category(), ErrorCategory::Table);
    }

    #[test]
    fn display_is_informative() {
        let msg = EcmError::invalid("capacity_ah", "must be positive").to_string();
        assert!(msg.contains("capacity_ah"), "got: {msg}");
        assert!(msg.contains("positive"), "got: {msg}");

        let msg = EcmError::SocOutOfRange {
            soc: 1.5,
            context: "initial SoC",
        }
        .to_string();
        assert!(msg.contains("1.5"), "got: {msg}");
        assert!(msg.contains("initial SoC"), "got: {msg}");

        let msg = EcmError::table("non-monotonic OCV").to_string();
        assert!(msg.contains("non-monotonic OCV"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(EcmError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
