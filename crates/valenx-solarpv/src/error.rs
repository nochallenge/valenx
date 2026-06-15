//! Error taxonomy for `valenx-solarpv`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, SolarPvError>`]. The variants are intentionally coarse — a
//! photovoltaic caller usually only cares about three things:
//!
//! 1. Did the caller pass a non-physical parameter — a non-positive
//!    temperature, a negative saturation current, an ideality factor
//!    outside its sensible band ([`SolarPvError::Invalid`])?
//! 2. Did a numerical root-find or maximum-power search fail to converge
//!    ([`SolarPvError::NoConvergence`])?
//! 3. Is a derived quantity undefined for the supplied operating point —
//!    a fill factor when `Voc*Isc` is zero, an efficiency when the
//!    incident power is zero ([`SolarPvError::Undefined`])?
//!
//! Use [`SolarPvError::code`] for stable log / telemetry tagging and
//! [`SolarPvError::category`] to bucket failures into Input / Numeric /
//! Domain without matching every variant. The shape mirrors the rest of
//! the workspace (`valenx-gears`'s `GearsError`, `valenx-popgen`'s
//! `PopgenError`).

use thiserror::Error;

/// Errors produced by `valenx-solarpv`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum SolarPvError {
    /// Caller passed a parameter the model cannot accept: a non-positive
    /// cell temperature, a negative photocurrent or saturation current,
    /// an ideality factor outside `(0, 10]`, a non-positive cell area,
    /// and so on. A property of the *call*, not of any numerical state.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"temperature_k"`, `"i0"`).
        what: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A numerical procedure — the implicit current root-find or the
    /// maximum-power-point search — exhausted its iteration budget
    /// without meeting the requested tolerance.
    #[error(
        "`{routine}` failed to converge after {iterations} iterations (residual {residual:e})"
    )]
    NoConvergence {
        /// Short name of the routine that gave up (e.g. `"solve_current"`).
        routine: &'static str,
        /// Number of iterations attempted before giving up.
        iterations: u32,
        /// Last residual magnitude reached.
        residual: f64,
    },

    /// A derived quantity is mathematically undefined for the supplied
    /// operating point: a fill factor when the `Voc * Isc` denominator is
    /// zero, an efficiency when the incident optical power
    /// (`irradiance * area`) is zero.
    #[error("undefined quantity `{quantity}`: {reason}")]
    Undefined {
        /// Quantity that could not be formed (e.g. `"fill_factor"`).
        quantity: &'static str,
        /// Why it is undefined here.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// Caller-supplied parameter is non-physical.
    Input,
    /// A numerical solver failed to converge.
    Numeric,
    /// A derived quantity is undefined for this operating point.
    Domain,
}

impl SolarPvError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"solarpv.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            SolarPvError::Invalid { .. } => "solarpv.invalid",
            SolarPvError::NoConvergence { .. } => "solarpv.no_convergence",
            SolarPvError::Undefined { .. } => "solarpv.undefined",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            SolarPvError::Invalid { .. } => "input",
            SolarPvError::NoConvergence { .. } => "numeric",
            SolarPvError::Undefined { .. } => "domain",
        }
    }

    /// Typed category enum (for callers that prefer to `match` rather
    /// than compare the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            SolarPvError::Invalid { .. } => ErrorCategory::Input,
            SolarPvError::NoConvergence { .. } => ErrorCategory::Numeric,
            SolarPvError::Undefined { .. } => ErrorCategory::Domain,
        }
    }

    /// Convenience constructor for [`SolarPvError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        SolarPvError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`SolarPvError::NoConvergence`].
    pub fn no_convergence(routine: &'static str, iterations: u32, residual: f64) -> Self {
        SolarPvError::NoConvergence {
            routine,
            iterations,
            residual,
        }
    }

    /// Convenience constructor for [`SolarPvError::Undefined`].
    pub fn undefined(quantity: &'static str, reason: impl Into<String>) -> Self {
        SolarPvError::Undefined {
            quantity,
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, SolarPvError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = SolarPvError::invalid("temperature_k", "must be > 0");
        assert_eq!(err.code(), "solarpv.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = SolarPvError::no_convergence("solve_current", 100, 1.0e-3);
        assert_eq!(err.code(), "solarpv.no_convergence");
        assert_eq!(err.category(), "numeric");
        assert_eq!(err.category_enum(), ErrorCategory::Numeric);

        let err = SolarPvError::undefined("fill_factor", "Voc*Isc is zero");
        assert_eq!(err.code(), "solarpv.undefined");
        assert_eq!(err.category(), "domain");
        assert_eq!(err.category_enum(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = SolarPvError::invalid("i0", "must be non-negative").to_string();
        assert!(msg.contains("i0"), "got: {msg}");
        assert!(msg.contains("non-negative"), "got: {msg}");

        let msg = SolarPvError::no_convergence("mpp", 200, 2.5e-6).to_string();
        assert!(msg.contains("mpp"), "got: {msg}");
        assert!(msg.contains("200"), "got: {msg}");

        let msg = SolarPvError::undefined("efficiency", "zero irradiance").to_string();
        assert!(msg.contains("efficiency"), "got: {msg}");
        assert!(msg.contains("zero irradiance"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(SolarPvError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
