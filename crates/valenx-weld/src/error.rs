//! Error taxonomy for `valenx-weld`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, WeldError>`]. The variants are intentionally coarse — a
//! weld-strength caller usually only cares about two things:
//!
//! 1. Was a geometric or load input non-physical — a non-positive leg,
//!    throat, length, thickness, force or allowable stress, or a value
//!    that is not finite ([`WeldError::BadParameter`])?
//! 2. Is the weld geometry itself degenerate — for example a requested
//!    effective throat that exceeds the geometric throat the leg can
//!    provide ([`WeldError::Degenerate`])?
//!
//! Use [`WeldError::code`] for stable log / telemetry tagging and
//! [`WeldError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-gears`' `GearsError`.

use thiserror::Error;

/// Errors raised by weld-strength calculations.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum WeldError {
    /// A geometric or load parameter was non-physical: non-positive,
    /// `NaN` or infinite. `name` identifies the offending parameter and
    /// `reason` is a human-readable explanation surfaced verbatim in the
    /// UI.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name (e.g. `"leg"`, `"length"`, `"force"`).
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// The weld geometry is internally inconsistent — for instance an
    /// effective throat larger than the geometric throat a fillet of the
    /// given leg can supply. A property of the *geometry*, not of a
    /// single out-of-range scalar.
    #[error("degenerate weld geometry: {0}")]
    Degenerate(String),
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every [`WeldError`] variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A user-supplied scalar input was out of range or not finite.
    Input,
    /// The weld geometry as a whole is inconsistent.
    Geometry,
}

impl WeldError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"weld.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            WeldError::BadParameter { .. } => "weld.bad_parameter",
            WeldError::Degenerate(_) => "weld.degenerate",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            WeldError::BadParameter { .. } => ErrorCategory::Input,
            WeldError::Degenerate(_) => ErrorCategory::Geometry,
        }
    }

    /// Convenience constructor for [`WeldError::BadParameter`].
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        WeldError::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`WeldError::Degenerate`].
    pub fn degenerate(reason: impl Into<String>) -> Self {
        WeldError::Degenerate(reason.into())
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, WeldError>;

/// Validate that `value` is a strictly positive, finite quantity.
///
/// Returns the value unchanged on success, or a
/// [`WeldError::BadParameter`] naming `name` when `value` is `NaN`,
/// infinite, zero or negative. This is the single gate every public
/// constructor and formula in the crate runs its scalar inputs through,
/// so the error message wording stays consistent.
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`] when `value` is not a strictly
/// positive finite number.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if value.is_nan() {
        return Err(WeldError::bad_parameter(name, "value is NaN"));
    }
    if value.is_infinite() {
        return Err(WeldError::bad_parameter(name, "value is infinite"));
    }
    if value <= 0.0 {
        return Err(WeldError::bad_parameter(
            name,
            format!("must be strictly positive, got {value}"),
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = WeldError::bad_parameter("leg", "must be positive");
        assert_eq!(err.code(), "weld.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = WeldError::degenerate("throat exceeds geometric throat");
        assert_eq!(err.code(), "weld.degenerate");
        assert_eq!(err.category(), ErrorCategory::Geometry);
    }

    #[test]
    fn display_is_informative() {
        let msg = WeldError::bad_parameter("length", "must be positive, got -3").to_string();
        assert!(msg.contains("length"), "got: {msg}");
        assert!(msg.contains("must be positive"), "got: {msg}");

        let msg = WeldError::degenerate("throat 9 mm > geometric throat 5 mm").to_string();
        assert!(msg.contains("degenerate"), "got: {msg}");
        assert!(msg.contains("throat"), "got: {msg}");
    }

    #[test]
    fn require_positive_accepts_positive_finite() {
        let v = require_positive("leg", 6.0).expect("6.0 is positive");
        assert!((v - 6.0).abs() < 1e-12, "got: {v}");
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        let err = require_positive("leg", 0.0).expect_err("zero rejected");
        assert_eq!(err.code(), "weld.bad_parameter");

        let err = require_positive("leg", -1.0).expect_err("negative rejected");
        assert_eq!(err.code(), "weld.bad_parameter");
    }

    #[test]
    fn require_positive_rejects_non_finite() {
        let err = require_positive("force", f64::NAN).expect_err("NaN rejected");
        assert!(err.to_string().contains("NaN"), "got: {err}");

        let err = require_positive("force", f64::INFINITY).expect_err("inf rejected");
        assert!(err.to_string().contains("infinite"), "got: {err}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(WeldError::degenerate("bad"));
        assert!(err.to_string().contains("bad"), "got: {err}");
    }
}
