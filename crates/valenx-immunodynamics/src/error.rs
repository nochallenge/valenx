//! Error taxonomy for `valenx-immunodynamics`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, ImmunoError>`]. The variants are intentionally coarse —
//! a within-host-dynamics caller usually only cares about two things:
//!
//! 1. Did the caller pass a nonsense argument — a negative rate
//!    constant, a non-positive integration step, an initial population
//!    below zero ([`ImmunoError::Invalid`])?
//! 2. Did the integration fail to produce a finite, bounded state — did
//!    the step size let the explicit scheme blow up to a non-finite
//!    value ([`ImmunoError::NotFinite`])?
//!
//! Use [`ImmunoError::code`] for stable log / telemetry tagging and
//! [`ImmunoError::category`] to bucket failures into Input / Numerical
//! without matching every variant. The pattern mirrors
//! `valenx-popgen`'s `PopgenError` and `valenx-sysbio`'s `SysbioError`.

/// Errors produced by `valenx-immunodynamics`.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ImmunoError {
    /// Caller passed an argument the model cannot accept: a negative
    /// rate constant, a non-positive time step, a negative initial
    /// population, an end-time that does not exceed the start time. A
    /// property of the *call*, not of a numerical breakdown.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"beta"`, `"dt"`, `"t_end"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The integration produced a non-finite (`NaN` / `±∞`) state. With
    /// a far-too-large step the explicit RK4 scheme can diverge; this
    /// reports the step at which the state stopped being finite so the
    /// caller can shrink `dt` and retry.
    #[error("integration diverged to a non-finite state at step {step} (t = {t}); reduce dt")]
    NotFinite {
        /// Step index at which the state first became non-finite.
        step: usize,
        /// Simulation time at that step.
        t: f64,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (bad argument).
    Input,
    /// The numerical method broke down (divergence to a non-finite
    /// state).
    Numerical,
}

impl ImmunoError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"immuno.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            ImmunoError::Invalid { .. } => "immuno.invalid",
            ImmunoError::NotFinite { .. } => "immuno.not_finite",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            ImmunoError::Invalid { .. } => "input",
            ImmunoError::NotFinite { .. } => "numerical",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            ImmunoError::Invalid { .. } => ErrorCategory::Input,
            ImmunoError::NotFinite { .. } => ErrorCategory::Numerical,
        }
    }

    /// Convenience constructor for [`ImmunoError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        ImmunoError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`ImmunoError::NotFinite`].
    pub fn not_finite(step: usize, t: f64) -> Self {
        ImmunoError::NotFinite { step, t }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, ImmunoError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = ImmunoError::invalid("beta", "must be non-negative");
        assert_eq!(err.code(), "immuno.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = ImmunoError::not_finite(42, 4.2);
        assert_eq!(err.code(), "immuno.not_finite");
        assert_eq!(err.category(), "numerical");
        assert_eq!(err.category_enum(), ErrorCategory::Numerical);
    }

    #[test]
    fn display_is_informative() {
        let msg = ImmunoError::invalid("dt", "step must be positive").to_string();
        assert!(msg.contains("dt"), "got: {msg}");
        assert!(msg.contains("step must be positive"), "got: {msg}");

        let msg = ImmunoError::not_finite(7, 1.5).to_string();
        assert!(msg.contains('7'), "got: {msg}");
        // The reported time should appear in the message.
        assert!(msg.contains("1.5"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(ImmunoError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
