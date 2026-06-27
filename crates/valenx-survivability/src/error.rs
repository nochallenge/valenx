//! The crate error type.
//!
//! Every fallible entry point returns [`SurvivabilityError`]. The
//! survivability models are full of divides and square roots that are only
//! defined on a physical domain (positive charge mass, positive stand-off,
//! positive plate area, …), so the policy is **fail loud, never panic**: a
//! degenerate or non-physical input becomes an [`Err`], not a `NaN` that
//! silently corrupts a downstream trade study.

use thiserror::Error;

/// An error from a survivability / protection calculation.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum SurvivabilityError {
    /// A scalar input was non-physical — outside the domain on which the
    /// model is defined (e.g. a non-positive charge mass, stand-off, mass,
    /// area, velocity, or a `NaN`/`±∞`). The message names the offending
    /// quantity and the value seen.
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),

    /// The requested scaled distance `Z = R / W^(1/3)` falls outside the
    /// validated fit range of the chosen empirical blast correlation. The
    /// correlations are curve fits to test data over a finite band of `Z`;
    /// extrapolating far past it is not trustworthy, so we refuse rather than
    /// return a confidently-wrong number.
    #[error(
        "scaled distance Z = {z:.4} m/kg^(1/3) is outside the validated range \
         [{min:.3}, {max:.3}] for the {model} fit"
    )]
    ScaledDistanceOutOfRange {
        /// The scaled distance that was requested.
        z: f64,
        /// Lower bound of the fit's validated range.
        min: f64,
        /// Upper bound of the fit's validated range.
        max: f64,
        /// Name of the empirical model whose range was exceeded.
        model: &'static str,
    },

    /// The underlying [`valenx_fem`] transient survivability solver failed.
    /// The structural-response models defer to that crate's Friedlander /
    /// Newmark integrator; its error is wrapped here verbatim.
    #[error("structural transient solve failed: {0}")]
    Transient(String),
}

impl From<valenx_fem::SurvivabilityError> for SurvivabilityError {
    fn from(e: valenx_fem::SurvivabilityError) -> Self {
        SurvivabilityError::Transient(e.to_string())
    }
}
