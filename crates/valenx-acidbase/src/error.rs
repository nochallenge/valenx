//! Error taxonomy for `valenx-acidbase`.
//!
//! Every fallible public function returns [`Result<_>`]. The variants are
//! deliberately coarse: a caller usually only needs to know whether an
//! input concentration / constant was out of its physical domain. The
//! [`validate_concentration`], [`validate_ka`], and
//! [`validate_kw`] helpers are the single source of truth for the input
//! checks and are reused across every model module.

use thiserror::Error;

/// Errors produced by `valenx-acidbase`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum AcidBaseError {
    /// A concentration was not a strictly positive, finite number.
    /// Molar concentrations must satisfy `C > 0` (a zero or negative
    /// concentration has no defined pH), and `NaN` / infinity are
    /// rejected outright.
    #[error("bad concentration `{name}` = {value}: must be finite and strictly positive (mol/L)")]
    BadConcentration {
        /// Parameter name, for the caller to locate the offending input.
        name: &'static str,
        /// The rejected value, echoed back verbatim.
        value: f64,
    },

    /// An acid / base dissociation constant `Ka` (or `Kb`) was not a
    /// strictly positive, finite number. Equilibrium constants are
    /// positive by construction.
    #[error("bad equilibrium constant `{name}` = {value}: must be finite and strictly positive")]
    BadConstant {
        /// Parameter name (`"Ka"`, `"Kb"`, `"Kw"`, ...).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A buffer was specified with a non-positive acid or conjugate-base
    /// concentration, so the `[A-] / [HA]` ratio in
    /// Henderson-Hasselbalch is undefined (a log of zero or a negative
    /// argument).
    #[error(
        "degenerate buffer: [HA] = {acid} mol/L, [A-] = {base} mol/L; \
         both must be finite and strictly positive"
    )]
    DegenerateBuffer {
        /// The weak-acid (`HA`) concentration supplied.
        acid: f64,
        /// The conjugate-base (`A-`) concentration supplied.
        base: f64,
    },
}

/// Convenience alias for `Result<T, AcidBaseError>`.
pub type Result<T> = std::result::Result<T, AcidBaseError>;

/// Validate that a concentration is finite and strictly positive,
/// returning it unchanged on success.
///
/// `name` is echoed into [`AcidBaseError::BadConcentration`] so the
/// caller can tell which argument was rejected.
///
/// # Errors
///
/// Returns [`AcidBaseError::BadConcentration`] when `value` is `NaN`,
/// infinite, zero, or negative.
pub fn validate_concentration(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(AcidBaseError::BadConcentration { name, value })
    }
}

/// Validate that an equilibrium constant (`Ka`, `Kb`) is finite and
/// strictly positive, returning it unchanged on success.
///
/// # Errors
///
/// Returns [`AcidBaseError::BadConstant`] when `value` is `NaN`,
/// infinite, zero, or negative.
pub fn validate_ka(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(AcidBaseError::BadConstant { name, value })
    }
}

/// Validate the water autoionization constant `Kw`.
///
/// Identical domain to [`validate_ka`] (finite, strictly positive) but
/// reported against the `"Kw"` name for clearer diagnostics.
///
/// # Errors
///
/// Returns [`AcidBaseError::BadConstant`] when `value` is `NaN`,
/// infinite, zero, or negative.
pub fn validate_kw(value: f64) -> Result<f64> {
    validate_ka("Kw", value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_concentration_passes_through() {
        let c = validate_concentration("C", 0.1).expect("0.1 M is valid");
        assert!((c - 0.1).abs() < 1e-15);
    }

    #[test]
    fn zero_concentration_rejected() {
        let err = validate_concentration("C", 0.0).expect_err("0 M has no pH");
        assert_eq!(
            err,
            AcidBaseError::BadConcentration {
                name: "C",
                value: 0.0
            }
        );
    }

    #[test]
    fn negative_and_nonfinite_concentration_rejected() {
        assert!(validate_concentration("C", -1.0).is_err());
        assert!(validate_concentration("C", f64::NAN).is_err());
        assert!(validate_concentration("C", f64::INFINITY).is_err());
    }

    #[test]
    fn constants_validate_independently() {
        assert!(validate_ka("Ka", 1.8e-5).is_ok());
        assert!(validate_ka("Ka", 0.0).is_err());
        assert!(validate_kw(1.0e-14).is_ok());
        assert!(validate_kw(-1.0).is_err());
    }

    #[test]
    fn error_messages_name_the_field() {
        let msg = AcidBaseError::BadConstant {
            name: "Ka",
            value: 0.0,
        }
        .to_string();
        assert!(msg.contains("Ka"), "message should name the field: {msg}");
    }
}
