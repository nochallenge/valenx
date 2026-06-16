//! Reinforced-concrete beam error taxonomy.
//!
//! Every fallible entry point in this crate returns [`RcBeamError`].
//! The variants carry the offending parameter name and a human-readable
//! reason; [`RcBeamError::code`] gives a stable kebab-cased identifier
//! and [`RcBeamError::category`] a coarse bucket for telemetry.

use thiserror::Error;

/// Errors raised while validating beam inputs or computing capacity.
#[derive(Debug, Error)]
pub enum RcBeamError {
    /// A geometric or material parameter was not strictly positive.
    ///
    /// Section width `b`, effective depth `d`, concrete strength
    /// `fc'`, steel yield `fy` and steel area `As` must all be `> 0`
    /// for the flexure model to be physically meaningful.
    #[error("non-positive parameter `{name}`: got {value}, must be > 0")]
    NonPositive {
        /// Parameter name (e.g. `"b"`, `"d"`, `"fc"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A parameter was not a finite number (was `NaN` or `±∞`).
    #[error("non-finite parameter `{name}`: got {value}")]
    NonFinite {
        /// Parameter name.
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// The strength-reduction factor `phi` was outside `(0, 1]`.
    ///
    /// ACI-318 strength-design phi factors are positive and never
    /// exceed unity (tension-controlled flexure uses `phi = 0.90`).
    #[error("strength-reduction factor phi out of range: got {value}, must be in (0, 1]")]
    PhiOutOfRange {
        /// The offending value.
        value: f64,
    },

    /// The equivalent stress-block depth `a` met or exceeded the
    /// effective depth `d`, so the lever arm `d - a/2` would be
    /// non-positive and the section geometry is degenerate.
    ///
    /// Physically this means the supplied steel area cannot be
    /// developed within the section: `As*fy` exceeds the concrete
    /// compression block the section can supply over depth `d`.
    #[error("degenerate section: stress-block depth a = {a} >= effective depth d = {d}")]
    StressBlockExceedsDepth {
        /// Computed stress-block depth `a`.
        a: f64,
        /// Effective depth `d`.
        d: f64,
    },

    /// The target moment is larger than the section can carry as a
    /// singly-reinforced, tension-controlled member.
    ///
    /// Sizing the tension steel for this moment would need a stress block
    /// at or beyond the effective depth (a negative discriminant in the
    /// design quadratic). The remedy — a deeper or wider section, stronger
    /// concrete, or compression steel — is outside this singly-reinforced
    /// model, so the request is rejected rather than approximated.
    #[error("target moment {target_moment} exceeds the singly-reinforced capacity {max_moment}")]
    MomentExceedsCapacity {
        /// The requested nominal moment.
        target_moment: f64,
        /// The largest nominal moment reachable singly-reinforced (at the
        /// `a = d` limit), `0.85 * fc' * b * d^2 / 2`.
        max_moment: f64,
    },
}

/// Coarse error category, for routing / telemetry.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Bad user-supplied input (geometry, materials, factors).
    Input,
    /// The algorithm domain was violated (degenerate geometry).
    Algorithm,
}

impl RcBeamError {
    /// Construct a [`RcBeamError::NonPositive`] / [`RcBeamError::NonFinite`]
    /// guard for a named scalar that must be a finite, strictly
    /// positive number.
    ///
    /// Returns `Ok(value)` when the check passes, so call sites can
    /// write `let b = RcBeamError::require_positive("b", b)?;`.
    ///
    /// # Errors
    ///
    /// Returns [`RcBeamError::NonFinite`] when `value` is `NaN` or
    /// infinite, and [`RcBeamError::NonPositive`] when `value <= 0`.
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, Self> {
        if !value.is_finite() {
            return Err(RcBeamError::NonFinite { name, value });
        }
        if value <= 0.0 {
            return Err(RcBeamError::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Stable kebab-cased identifier for the variant.
    pub fn code(&self) -> &'static str {
        match self {
            RcBeamError::NonPositive { .. } => "rcbeam.non-positive",
            RcBeamError::NonFinite { .. } => "rcbeam.non-finite",
            RcBeamError::PhiOutOfRange { .. } => "rcbeam.phi-out-of-range",
            RcBeamError::StressBlockExceedsDepth { .. } => "rcbeam.stress-block-exceeds-depth",
            RcBeamError::MomentExceedsCapacity { .. } => "rcbeam.moment-exceeds-capacity",
        }
    }

    /// Coarse category bucket.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RcBeamError::NonPositive { .. }
            | RcBeamError::NonFinite { .. }
            | RcBeamError::PhiOutOfRange { .. } => ErrorCategory::Input,
            RcBeamError::StressBlockExceedsDepth { .. }
            | RcBeamError::MomentExceedsCapacity { .. } => ErrorCategory::Algorithm,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        let v = RcBeamError::require_positive("b", 250.0).unwrap();
        assert!((v - 250.0).abs() < 1e-12);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        let z = RcBeamError::require_positive("d", 0.0).unwrap_err();
        assert_eq!(z.code(), "rcbeam.non-positive");
        let n = RcBeamError::require_positive("d", -1.0).unwrap_err();
        assert_eq!(n.code(), "rcbeam.non-positive");
    }

    #[test]
    fn require_positive_rejects_non_finite() {
        let nan = RcBeamError::require_positive("fc", f64::NAN).unwrap_err();
        assert_eq!(nan.code(), "rcbeam.non-finite");
        let inf = RcBeamError::require_positive("fc", f64::INFINITY).unwrap_err();
        assert_eq!(inf.code(), "rcbeam.non-finite");
    }

    #[test]
    fn categories_route_correctly() {
        assert_eq!(
            RcBeamError::NonPositive {
                name: "b",
                value: 0.0
            }
            .category(),
            ErrorCategory::Input
        );
        assert_eq!(
            RcBeamError::StressBlockExceedsDepth { a: 600.0, d: 500.0 }.category(),
            ErrorCategory::Algorithm
        );
    }
}
