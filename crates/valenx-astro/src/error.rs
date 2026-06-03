//! Error type for the ascent simulator.

use thiserror::Error;

/// Shorthand for `Result<T, AstroError>`.
pub type Result<T> = core::result::Result<T, AstroError>;

/// Anything that can go wrong building or running an ascent case.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum AstroError {
    /// A vehicle was supplied with no stages.
    #[error("vehicle has no stages")]
    NoStages,

    /// A stage carried a non-physical mass (negative, zero, or NaN).
    #[error("stage {index}: invalid mass ({field} = {value})")]
    InvalidMass {
        /// Zero-based stage index.
        index: usize,
        /// Which mass field was bad (`dry_mass` / `propellant_mass`).
        field: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A stage carried a non-physical propulsion parameter.
    #[error("stage {index}: invalid propulsion ({field} = {value})")]
    InvalidPropulsion {
        /// Zero-based stage index.
        index: usize,
        /// Which field was bad (`thrust_vac`, `isp_vac`, …).
        field: &'static str,
        /// The offending value.
        value: f64,
    },

    /// The aerodynamic reference area or drag table was invalid.
    #[error("invalid aerodynamics: {0}")]
    InvalidAero(&'static str),

    /// A guidance / launch parameter was out of range.
    #[error("invalid guidance: {0}")]
    InvalidGuidance(&'static str),

    /// A state vector or set of orbital elements was non-physical in a
    /// way that would otherwise produce a silent `NaN`/`Inf` result —
    /// e.g. a zero position/velocity vector, zero angular momentum
    /// (rectilinear motion), the parabolic energy singularity, or a
    /// non-positive semi-latus rectum. Carries a short reason.
    #[error("non-physical state: {0}")]
    NonPhysicalState(&'static str),

    /// The integration step or duration was non-positive.
    #[error("invalid integration setting: {0}")]
    InvalidIntegration(&'static str),

    /// The simulation hit its wall-clock step budget before any
    /// termination condition (orbit / impact / burnout-coast end) fired.
    #[error("simulation exceeded the maximum step count ({0})")]
    StepBudgetExceeded(u64),

    /// A count or size argument exceeded its hard upper bound. Used to
    /// reject inputs that would otherwise drive the simulator into an
    /// effectively-unbounded loop or allocation (e.g. an absurd step or
    /// sample count).
    #[error("{what} = {value} exceeds the maximum allowed ({max})")]
    OutOfRange {
        /// Which quantity was out of range (e.g. `"steps"`, `"samples"`).
        what: &'static str,
        /// The offending value.
        value: u64,
        /// The hard upper bound.
        max: u64,
    },

    /// A scalar design / mission parameter was non-physical in a way that
    /// would otherwise feed a silent `NaN`/`Inf` into a `√`/`ln`/`exp` or
    /// a division (e.g. a non-finite or out-of-`(0,1)` structural mass
    /// fraction, a non-positive specific impulse, density, ballistic
    /// coefficient, nose radius or flight-path angle, or an infeasible
    /// `Δv` budget that has no positive payload solution). Carries a short
    /// reason.
    #[error("invalid parameter: {0}")]
    InvalidParameter(&'static str),

    /// An analytic / iterative solution did not converge, or its
    /// closed-form is singular at the requested point — e.g. the
    /// Clohessy–Wiltshire two-impulse transfer when the transfer time is a
    /// whole number of orbital periods (`sin(nT) = 0`), where the
    /// position-to-velocity block of the state-transition matrix is not
    /// invertible. Carries a short reason.
    #[error("non-convergent / singular solution: {0}")]
    NonConvergent(&'static str),
}
