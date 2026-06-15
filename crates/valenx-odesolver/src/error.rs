//! Error taxonomy for the ODE-integration routines.
//!
//! Every fallible entry point in this crate returns [`OdeError`]. The
//! variants are constructed only through the validated constructors below,
//! so a constructed [`OdeError`] always describes a genuine, checked
//! violation (a non-finite parameter, a non-positive step, a zero step
//! count, or a ragged state vector).

use thiserror::Error;

/// Errors raised while configuring or running an integration.
///
/// Construct values through the associated constructors
/// ([`OdeError::bad_step`], [`OdeError::bad_step_count`],
/// [`OdeError::non_finite`], [`OdeError::dimension_mismatch`]) rather than
/// building the struct variants directly: the constructors perform the
/// validation that gives each variant its meaning.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum OdeError {
    /// The fixed step size `dt` was not a strictly positive, finite number.
    #[error("invalid step size `dt = {dt}`: {reason}")]
    BadStep {
        /// The offending step value.
        dt: f64,
        /// Human-readable reason the step was rejected.
        reason: &'static str,
    },

    /// The requested number of steps was zero.
    ///
    /// A fixed-step integration must advance at least one step; a zero-step
    /// request is almost always a caller bug rather than a no-op intent.
    #[error("step count must be >= 1, got 0")]
    BadStepCount,

    /// A supplied initial value (or a value produced mid-integration) was
    /// not finite (it was `NaN` or `±∞`).
    #[error("non-finite value for `{name}`: {value}")]
    NonFinite {
        /// Name of the quantity that was non-finite.
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// Two state vectors that must share a length did not.
    ///
    /// Raised when the derivative closure of a system returns a vector whose
    /// length differs from the supplied initial state, which would make the
    /// per-component update ill-defined.
    #[error("dimension mismatch: expected length {expected}, got {actual}")]
    DimensionMismatch {
        /// The length the state vector was expected to have.
        expected: usize,
        /// The length actually observed.
        actual: usize,
    },
}

impl OdeError {
    /// Validate a fixed step size, returning an [`OdeError::BadStep`] when it
    /// is not strictly positive and finite.
    ///
    /// A step is rejected when it is `NaN`, infinite, zero, or negative.
    /// Fixed-step explicit integrators require `dt > 0` to march forward in
    /// time with a well-defined truncation error.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_odesolver::OdeError;
    /// assert!(OdeError::bad_step(0.0).is_some());
    /// assert!(OdeError::bad_step(-1.0).is_some());
    /// assert!(OdeError::bad_step(f64::NAN).is_some());
    /// assert!(OdeError::bad_step(1e-3).is_none());
    /// ```
    #[must_use]
    pub fn bad_step(dt: f64) -> Option<OdeError> {
        if dt.is_nan() {
            Some(OdeError::BadStep {
                dt,
                reason: "step is NaN",
            })
        } else if dt.is_infinite() {
            Some(OdeError::BadStep {
                dt,
                reason: "step is infinite",
            })
        } else if dt <= 0.0 {
            Some(OdeError::BadStep {
                dt,
                reason: "step must be strictly positive",
            })
        } else {
            None
        }
    }

    /// Validate a step count, returning [`OdeError::BadStepCount`] when it is
    /// zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_odesolver::OdeError;
    /// assert!(OdeError::bad_step_count(0).is_some());
    /// assert!(OdeError::bad_step_count(1).is_none());
    /// ```
    #[must_use]
    pub fn bad_step_count(n: usize) -> Option<OdeError> {
        if n == 0 {
            Some(OdeError::BadStepCount)
        } else {
            None
        }
    }

    /// Validate that `value` is finite, returning [`OdeError::NonFinite`]
    /// (tagged with `name`) otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_odesolver::OdeError;
    /// assert!(OdeError::non_finite("y0", f64::INFINITY).is_some());
    /// assert!(OdeError::non_finite("y0", 3.5).is_none());
    /// ```
    #[must_use]
    pub fn non_finite(name: &'static str, value: f64) -> Option<OdeError> {
        if value.is_finite() {
            None
        } else {
            Some(OdeError::NonFinite { name, value })
        }
    }

    /// Validate that two lengths match, returning
    /// [`OdeError::DimensionMismatch`] otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_odesolver::OdeError;
    /// assert!(OdeError::dimension_mismatch(3, 2).is_some());
    /// assert!(OdeError::dimension_mismatch(3, 3).is_none());
    /// ```
    #[must_use]
    pub fn dimension_mismatch(expected: usize, actual: usize) -> Option<OdeError> {
        if expected == actual {
            None
        } else {
            Some(OdeError::DimensionMismatch { expected, actual })
        }
    }

    /// Stable kebab-cased identifier for logging / matching, independent of
    /// the human-readable message.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            OdeError::BadStep { .. } => "odesolver.bad_step",
            OdeError::BadStepCount => "odesolver.bad_step_count",
            OdeError::NonFinite { .. } => "odesolver.non_finite",
            OdeError::DimensionMismatch { .. } => "odesolver.dimension_mismatch",
        }
    }
}
