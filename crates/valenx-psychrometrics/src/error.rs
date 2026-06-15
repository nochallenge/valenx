//! Error taxonomy for `valenx-psychrometrics`.
//!
//! Every fallible public function returns
//! [`Result<_, PsychroError>`]. The variants are intentionally coarse —
//! a psychrometrics caller usually only cares about two things:
//! whether an argument was outside its physical domain
//! ([`PsychroError::BadParameter`]), or whether the requested state is
//! itself physically impossible — for example a vapour pressure that
//! meets or exceeds the total pressure, which would drive the humidity
//! ratio to infinity ([`PsychroError::Unphysical`]).
//!
//! Use [`PsychroError::code`] for stable log / telemetry tagging and
//! [`PsychroError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-hvac`'s `HvacError` and
//! `valenx-springs`'s `SpringsError`.

use thiserror::Error;

/// Errors raised by the psychrometrics crate.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum PsychroError {
    /// Caller passed an argument outside its physical domain: a
    /// non-positive absolute pressure, a relative humidity outside
    /// `[0, 1]`, a negative humidity ratio, or a temperature below
    /// absolute zero. A property of the *call*, not of a derived state.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Logical parameter name (e.g. `"pressure_pa"`, `"rh"`).
        name: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The requested moist-air state is physically impossible — most
    /// commonly a partial vapour pressure that meets or exceeds the
    /// total barometric pressure, which would make the humidity ratio
    /// `w = 0.622 pv / (p - pv)` diverge.
    #[error("unphysical state: {0}")]
    Unphysical(String),
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is outside its physical domain.
    Input,
    /// The derived moist-air state is physically impossible.
    Domain,
}

impl PsychroError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"psychro.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            PsychroError::BadParameter { .. } => "psychro.bad_parameter",
            PsychroError::Unphysical(_) => "psychro.unphysical",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            PsychroError::BadParameter { .. } => ErrorCategory::Input,
            PsychroError::Unphysical(_) => ErrorCategory::Domain,
        }
    }

    /// Convenience constructor for [`PsychroError::BadParameter`].
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        PsychroError::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PsychroError::Unphysical`].
    pub fn unphysical(reason: impl Into<String>) -> Self {
        PsychroError::Unphysical(reason.into())
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, PsychroError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = PsychroError::bad_parameter("rh", "must be in [0, 1]");
        assert_eq!(err.code(), "psychro.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = PsychroError::unphysical("pv >= p");
        assert_eq!(err.code(), "psychro.unphysical");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = PsychroError::bad_parameter("pressure_pa", "must be positive").to_string();
        assert!(msg.contains("pressure_pa"), "got: {msg}");
        assert!(msg.contains("must be positive"), "got: {msg}");

        let msg = PsychroError::unphysical("pv >= p").to_string();
        assert!(msg.contains("pv >= p"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(PsychroError::bad_parameter("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
