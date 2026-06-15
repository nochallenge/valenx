//! Error taxonomy for `valenx-heatpump`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, HeatPumpError>`]. The variants are intentionally coarse —
//! a heat-pump caller usually only cares about three things:
//!
//! 1. Did the caller pass a physically meaningless argument — a
//!    non-positive absolute temperature, a negative efficiency, a
//!    fraction outside `[0, 1]` ([`HeatPumpError::Invalid`])?
//! 2. Is the temperature lift degenerate — the hot and cold reservoirs
//!    are equal (or inverted), so the Carnot COP would divide by zero
//!    or go negative ([`HeatPumpError::DegenerateLift`])?
//! 3. Did a numerical solver fail to bracket or converge on a root —
//!    the balance-point search found no sign change in the supplied
//!    interval ([`HeatPumpError::NoConvergence`])?
//!
//! Use [`HeatPumpError::code`] for stable log / telemetry tagging and
//! [`HeatPumpError::category`] to bucket failures without matching every
//! variant. The pattern mirrors the other Valenx physics crates'
//! `*Error` enums (e.g. `valenx-springs`'s `SpringsError`).

use thiserror::Error;

/// Errors produced by `valenx-heatpump`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum HeatPumpError {
    /// Caller passed an argument outside its physical domain: a
    /// non-positive absolute temperature, a negative or super-unity
    /// Carnot fraction, a non-positive slope, an empty interval, etc. A
    /// property of the *call*, not of a numerical search.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"t_hot_k"`, `"carnot_fraction"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The temperature lift `T_h - T_c` is zero or negative, so a Carnot
    /// coefficient of performance is undefined (division by zero) or
    /// physically meaningless. Carries both reservoir temperatures in
    /// kelvin so the caller can report them verbatim.
    #[error(
        "degenerate temperature lift: T_hot = {t_hot_k} K must be strictly greater \
         than T_cold = {t_cold_k} K"
    )]
    DegenerateLift {
        /// Hot-reservoir absolute temperature supplied, in kelvin.
        t_hot_k: f64,
        /// Cold-reservoir absolute temperature supplied, in kelvin.
        t_cold_k: f64,
    },

    /// A numerical root-finder (the balance-point solver) could not find
    /// a solution: the load and capacity curves do not cross inside the
    /// bracketing interval, or the bracket endpoints share a sign.
    #[error("balance-point solver failed to converge: {reason}")]
    NoConvergence {
        /// Human-readable reason (e.g. `"load exceeds capacity across the
        /// whole interval"`).
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is outside its physical domain.
    Input,
    /// A numerical solver failed to converge.
    Numerical,
}

impl HeatPumpError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"heatpump.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            HeatPumpError::Invalid { .. } => "heatpump.invalid",
            HeatPumpError::DegenerateLift { .. } => "heatpump.degenerate_lift",
            HeatPumpError::NoConvergence { .. } => "heatpump.no_convergence",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            HeatPumpError::Invalid { .. } | HeatPumpError::DegenerateLift { .. } => "input",
            HeatPumpError::NoConvergence { .. } => "numerical",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            HeatPumpError::Invalid { .. } | HeatPumpError::DegenerateLift { .. } => {
                ErrorCategory::Input
            }
            HeatPumpError::NoConvergence { .. } => ErrorCategory::Numerical,
        }
    }

    /// Convenience constructor for [`HeatPumpError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        HeatPumpError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`HeatPumpError::DegenerateLift`].
    pub fn degenerate_lift(t_hot_k: f64, t_cold_k: f64) -> Self {
        HeatPumpError::DegenerateLift { t_hot_k, t_cold_k }
    }

    /// Convenience constructor for [`HeatPumpError::NoConvergence`].
    pub fn no_convergence(reason: impl Into<String>) -> Self {
        HeatPumpError::NoConvergence {
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, HeatPumpError>;

/// Validate that an absolute temperature is finite and strictly above
/// absolute zero.
///
/// Returns [`HeatPumpError::Invalid`] otherwise. Shared by every
/// constructor that accepts a kelvin temperature so the domain check is
/// written once.
///
/// # Examples
///
/// ```
/// use valenx_heatpump::error::check_temperature_k;
/// assert!(check_temperature_k("t_hot_k", 300.0).is_ok());
/// assert!(check_temperature_k("t_hot_k", 0.0).is_err());
/// assert!(check_temperature_k("t_hot_k", f64::NAN).is_err());
/// ```
pub fn check_temperature_k(what: &'static str, t_k: f64) -> Result<f64> {
    if !t_k.is_finite() {
        return Err(HeatPumpError::invalid(what, "must be a finite number"));
    }
    if t_k <= 0.0 {
        return Err(HeatPumpError::invalid(
            what,
            format!("must be a positive absolute temperature in kelvin, got {t_k}"),
        ));
    }
    Ok(t_k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = HeatPumpError::invalid("t_hot_k", "must be positive");
        assert_eq!(err.code(), "heatpump.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = HeatPumpError::degenerate_lift(290.0, 300.0);
        assert_eq!(err.code(), "heatpump.degenerate_lift");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = HeatPumpError::no_convergence("no sign change");
        assert_eq!(err.code(), "heatpump.no_convergence");
        assert_eq!(err.category(), "numerical");
        assert_eq!(err.category_enum(), ErrorCategory::Numerical);
    }

    #[test]
    fn display_is_informative() {
        let msg = HeatPumpError::degenerate_lift(290.0, 300.0).to_string();
        assert!(msg.contains("290"), "got: {msg}");
        assert!(msg.contains("300"), "got: {msg}");

        let msg = HeatPumpError::invalid("carnot_fraction", "out of range").to_string();
        assert!(msg.contains("carnot_fraction"), "got: {msg}");
        assert!(msg.contains("out of range"), "got: {msg}");

        let msg = HeatPumpError::no_convergence("flat residual").to_string();
        assert!(msg.contains("flat residual"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(HeatPumpError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }

    #[test]
    fn temperature_check_accepts_positive_finite() {
        assert_eq!(check_temperature_k("t", 273.15).unwrap(), 273.15);
    }

    #[test]
    fn temperature_check_rejects_non_positive_and_nonfinite() {
        assert!(check_temperature_k("t", 0.0).is_err());
        assert!(check_temperature_k("t", -5.0).is_err());
        assert!(check_temperature_k("t", f64::NAN).is_err());
        assert!(check_temperature_k("t", f64::INFINITY).is_err());
    }
}
