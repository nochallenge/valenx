//! Wind-turbine error taxonomy.

use thiserror::Error;

/// Errors raised by wind-turbine aerodynamic calculations.
///
/// Every fallible constructor and free function in this crate returns
/// [`WindTurbineError`]. The variants distinguish a plainly invalid
/// physical input (a non-positive density, a negative wind speed) from a
/// logically inconsistent power-curve configuration (cut-out below
/// cut-in, a power coefficient exceeding the Betz limit).
#[derive(Debug, Error)]
pub enum WindTurbineError {
    /// A scalar parameter fell outside its admissible physical range.
    ///
    /// `name` is the parameter; `reason` explains the violated bound
    /// (e.g. `"must be > 0"`).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name (e.g. `"air_density"`, `"wind_speed"`).
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A power coefficient `Cp` exceeded the Betz limit `16/27`.
    ///
    /// An ideal actuator disc cannot extract more than `16/27 ~ 0.593`
    /// of the wind's kinetic-energy flux; a `Cp` above that is
    /// physically impossible (momentum theory), so it is rejected
    /// rather than silently producing super-Betz power.
    #[error("power coefficient Cp = {cp} exceeds the Betz limit {betz}")]
    AboveBetz {
        /// The offending power coefficient.
        cp: f64,
        /// The Betz limit `16/27`.
        betz: f64,
    },

    /// The characteristic wind speeds of a power curve are out of order.
    ///
    /// A well-formed curve requires
    /// `0 < cut_in < rated < cut_out`. This variant reports the first
    /// ordering constraint that was violated.
    #[error("inconsistent power curve: {0}")]
    InconsistentCurve(String),
}

/// Coarse category for a [`WindTurbineError`], for UI grouping / logging.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid physical input value.
    Input,
    /// A tunable configuration knob is inconsistent.
    Config,
    /// A modelling / domain constraint was violated.
    Algorithm,
}

impl WindTurbineError {
    /// Stable kebab-cased identifier, suitable for matching in logs or
    /// machine-readable error reporting.
    ///
    /// ```
    /// use valenx_windturbine::WindTurbineError;
    /// let e = WindTurbineError::BadParameter {
    ///     name: "wind_speed",
    ///     reason: "must be >= 0".to_string(),
    /// };
    /// assert_eq!(e.code(), "windturbine.bad_parameter");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            WindTurbineError::BadParameter { .. } => "windturbine.bad_parameter",
            WindTurbineError::AboveBetz { .. } => "windturbine.above_betz",
            WindTurbineError::InconsistentCurve(_) => "windturbine.inconsistent_curve",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            WindTurbineError::BadParameter { .. } => ErrorCategory::Input,
            WindTurbineError::AboveBetz { .. } => ErrorCategory::Algorithm,
            WindTurbineError::InconsistentCurve(_) => ErrorCategory::Config,
        }
    }
}
