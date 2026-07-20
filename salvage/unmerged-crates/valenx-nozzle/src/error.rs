//! Error type for the isentropic-nozzle calculations.
//!
//! Every public entry point validates its arguments before touching a
//! `sqrt`, `powf`, or division, so a bad input surfaces as a typed
//! [`NozzleError`] rather than a silent `NaN`/`Inf`. Each variant carries
//! a stable [`NozzleError::code`] string suitable for logging or matching
//! in a UI layer.

use thiserror::Error;

/// Shorthand for `Result<T, NozzleError>`.
pub type Result<T> = core::result::Result<T, NozzleError>;

/// Anything that can go wrong validating or evaluating an isentropic-flow
/// quantity.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in a future
/// release without it being a breaking change, so downstream `match` arms
/// must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum NozzleError {
    /// An argument was not a finite number (it was `NaN` or `±∞`).
    #[error("{name} must be finite, got {value}")]
    NotFinite {
        /// Which argument was bad (e.g. `"mach"`, `"gamma"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// The Mach number was negative. Mach is a non-negative ratio of flow
    /// speed to the local speed of sound; subsonic, sonic, and supersonic
    /// flow all have `M >= 0`.
    #[error("mach must be >= 0, got {value}")]
    NegativeMach {
        /// The offending Mach value.
        value: f64,
    },

    /// The ratio of specific heats `gamma` was outside the physical open
    /// interval `(1, 5/3]`. Real gases sit between the diatomic value
    /// `7/5 = 1.4` (air) and the monatomic limit `5/3 ≈ 1.667`; `gamma`
    /// must exceed 1 for the relations to be defined (the exponents carry a
    /// `gamma - 1` in the denominator).
    #[error("gamma must be in (1, 5/3], got {value}")]
    GammaOutOfRange {
        /// The offending `gamma` value.
        value: f64,
    },

    /// An area ratio `A/A*` was supplied that is below the choked minimum of
    /// `1.0`. The area-Mach relation has a global minimum of unity at the
    /// sonic throat (`M = 1`); no real Mach number maps to `A/A* < 1`.
    #[error("area ratio A/A* must be >= 1, got {value}")]
    AreaRatioBelowChoke {
        /// The offending area ratio.
        value: f64,
    },

    /// A pressure ratio `p/p0` handed to an inverse routine was outside the
    /// physically attainable open interval `(0, 1]`: the static pressure of
    /// an expanding flow is always positive and never exceeds the
    /// stagnation pressure.
    #[error("pressure ratio p/p0 must be in (0, 1], got {value}")]
    PressureRatioOutOfRange {
        /// The offending pressure ratio.
        value: f64,
    },

    /// The Newton iteration in an inverse routine failed to converge within
    /// its iteration budget. Carries the residual it reached.
    #[error("inversion did not converge (residual {residual:e} after {iters} iterations)")]
    NotConverged {
        /// The absolute residual reached at the last step.
        residual: f64,
        /// How many iterations were spent.
        iters: u32,
    },
}

impl NozzleError {
    /// A short, stable machine-readable code for this error.
    ///
    /// Codes are part of the public contract and will not change for an
    /// existing variant; they are handy for log filters or for surfacing a
    /// category to a UI without formatting the full message.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            NozzleError::NotFinite { .. } => "not_finite",
            NozzleError::NegativeMach { .. } => "negative_mach",
            NozzleError::GammaOutOfRange { .. } => "gamma_out_of_range",
            NozzleError::AreaRatioBelowChoke { .. } => "area_ratio_below_choke",
            NozzleError::PressureRatioOutOfRange { .. } => "pressure_ratio_out_of_range",
            NozzleError::NotConverged { .. } => "not_converged",
        }
    }
}

/// Validate a `gamma` (ratio of specific heats): finite and in `(1, 5/3]`.
///
/// # Errors
///
/// Returns [`NozzleError::NotFinite`] if `gamma` is `NaN`/`±∞`, or
/// [`NozzleError::GammaOutOfRange`] if it is outside `(1, 5/3]`.
pub fn check_gamma(gamma: f64) -> Result<()> {
    if !gamma.is_finite() {
        return Err(NozzleError::NotFinite {
            name: "gamma",
            value: gamma,
        });
    }
    // 5/3 is the monatomic-ideal-gas upper limit; allow a hair of slack so a
    // literal `5.0/3.0` (which rounds to 1.6666666666666667) is accepted.
    if gamma <= 1.0 || gamma > 5.0 / 3.0 + 1e-12 {
        return Err(NozzleError::GammaOutOfRange { value: gamma });
    }
    Ok(())
}

/// Validate a Mach number: finite and non-negative.
///
/// # Errors
///
/// Returns [`NozzleError::NotFinite`] if `mach` is `NaN`/`±∞`, or
/// [`NozzleError::NegativeMach`] if it is negative.
pub fn check_mach(mach: f64) -> Result<()> {
    if !mach.is_finite() {
        return Err(NozzleError::NotFinite {
            name: "mach",
            value: mach,
        });
    }
    if mach < 0.0 {
        return Err(NozzleError::NegativeMach { value: mach });
    }
    Ok(())
}
