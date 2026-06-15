//! Error taxonomy for the refrigeration crate.
//!
//! Every fallible entry point returns [`Result<T>`], an alias for
//! [`core::result::Result<T, RefrigError>`]. The error type exposes a
//! stable kebab-cased [`code`](RefrigError::code) and a coarse
//! [`category`](RefrigError::category) for telemetry, plus a small set
//! of validated constructors used throughout the crate so that domain
//! invariants (positive temperatures, a thermodynamically ordered
//! temperature lift, non-negative heat flows, a strictly positive
//! compressor work) are checked in exactly one place.

use thiserror::Error;

/// Convenience alias for results produced by this crate.
pub type Result<T> = core::result::Result<T, RefrigError>;

/// Errors raised by the refrigeration thermodynamics routines.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RefrigError {
    /// A numeric parameter was outside its physically valid domain
    /// (for example a non-positive absolute temperature, a negative
    /// heat flow, or a non-finite enthalpy).
    #[error("bad parameter `{name}` = {value}: {reason}")]
    BadParameter {
        /// The offending parameter name.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
        /// Why the value is rejected.
        reason: &'static str,
    },

    /// The hot- and cold-side absolute temperatures do not form a valid
    /// refrigeration temperature lift: a refrigerator or heat pump moves
    /// heat from cold to hot, so the condensing (hot) temperature must be
    /// strictly greater than the evaporating (cold) temperature.
    #[error(
        "invalid temperature lift: hot side {t_hot} K must exceed cold side {t_cold} K \
         (Th - Tc must be > 0)"
    )]
    InvalidLift {
        /// Hot-side (condenser / heat-rejection) absolute temperature, in kelvin.
        t_hot: f64,
        /// Cold-side (evaporator / refrigerated-space) absolute temperature, in kelvin.
        t_cold: f64,
    },

    /// The compressor work is zero or negative, so a coefficient of
    /// performance is undefined (it would divide by zero or yield a
    /// nonsensical negative value).
    #[error("compressor work {0} kJ/kg must be strictly positive to define a COP")]
    NonPositiveWork(f64),

    /// The cycle enthalpies are not thermodynamically ordered for a
    /// vapor-compression refrigeration loop, so the derived effects would
    /// be negative. The reason string names the specific violation.
    #[error("inconsistent cycle enthalpies: {0}")]
    InconsistentCycle(&'static str),
}

impl RefrigError {
    /// A stable kebab-cased identifier suitable for logs and metrics.
    ///
    /// The string is part of the crate's public contract and will not
    /// change for an existing variant.
    pub fn code(&self) -> &'static str {
        match self {
            RefrigError::BadParameter { .. } => "refrigeration.bad_parameter",
            RefrigError::InvalidLift { .. } => "refrigeration.invalid_lift",
            RefrigError::NonPositiveWork(_) => "refrigeration.non_positive_work",
            RefrigError::InconsistentCycle(_) => "refrigeration.inconsistent_cycle",
        }
    }

    /// A coarse category for grouping errors in dashboards.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RefrigError::BadParameter { .. } => ErrorCategory::Input,
            RefrigError::InvalidLift { .. } => ErrorCategory::Input,
            RefrigError::NonPositiveWork(_) => ErrorCategory::Input,
            RefrigError::InconsistentCycle(_) => ErrorCategory::Thermodynamics,
        }
    }

    /// Validate that `value` is a strictly positive, finite number.
    ///
    /// Used for absolute temperatures and the compressor work as a raw
    /// scalar. Returns the value unchanged on success.
    pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(RefrigError::BadParameter {
                name,
                value,
                reason: "must be a finite, strictly positive value",
            })
        }
    }

    /// Validate that `value` is a finite, non-negative number.
    ///
    /// Used for heat flows, which may legitimately be zero. Returns the
    /// value unchanged on success.
    pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(RefrigError::BadParameter {
                name,
                value,
                reason: "must be a finite, non-negative value",
            })
        }
    }

    /// Validate that `value` is simply finite (no NaN / infinity).
    ///
    /// Used for specific enthalpies, where only the differences carry
    /// physical meaning so any finite datum is admissible. Returns the
    /// value unchanged on success.
    pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(RefrigError::BadParameter {
                name,
                value,
                reason: "must be a finite value (not NaN or infinity)",
            })
        }
    }
}

/// Coarse bucket used to triage [`RefrigError`] values.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-domain input.
    Input,
    /// The supplied thermodynamic state is internally inconsistent.
    Thermodynamics,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_guard_accepts_and_rejects() {
        assert_eq!(RefrigError::require_positive("t", 300.0).unwrap(), 300.0);
        assert!(RefrigError::require_positive("t", 0.0).is_err());
        assert!(RefrigError::require_positive("t", -1.0).is_err());
        assert!(RefrigError::require_positive("t", f64::NAN).is_err());
        assert!(RefrigError::require_positive("t", f64::INFINITY).is_err());
    }

    #[test]
    fn non_negative_guard_allows_zero() {
        assert_eq!(RefrigError::require_non_negative("q", 0.0).unwrap(), 0.0);
        assert!(RefrigError::require_non_negative("q", -0.001).is_err());
        assert!(RefrigError::require_non_negative("q", f64::NAN).is_err());
    }

    #[test]
    fn finite_guard_allows_negative() {
        assert_eq!(RefrigError::require_finite("h", -42.0).unwrap(), -42.0);
        assert!(RefrigError::require_finite("h", f64::INFINITY).is_err());
    }

    #[test]
    fn codes_and_categories_are_stable() {
        let bad = RefrigError::BadParameter {
            name: "x",
            value: -1.0,
            reason: "r",
        };
        assert_eq!(bad.code(), "refrigeration.bad_parameter");
        assert_eq!(bad.category(), ErrorCategory::Input);

        let lift = RefrigError::InvalidLift {
            t_hot: 280.0,
            t_cold: 300.0,
        };
        assert_eq!(lift.code(), "refrigeration.invalid_lift");
        assert_eq!(lift.category(), ErrorCategory::Input);

        let work = RefrigError::NonPositiveWork(0.0);
        assert_eq!(work.code(), "refrigeration.non_positive_work");
        assert_eq!(work.category(), ErrorCategory::Input);

        let cyc = RefrigError::InconsistentCycle("h2 <= h1");
        assert_eq!(cyc.code(), "refrigeration.inconsistent_cycle");
        assert_eq!(cyc.category(), ErrorCategory::Thermodynamics);
    }

    #[test]
    fn display_strings_render() {
        let e = RefrigError::NonPositiveWork(-3.0);
        let s = format!("{e}");
        assert!(s.contains("-3"));
        assert!(s.contains("COP"));
    }
}
