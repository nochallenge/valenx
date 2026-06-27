//! Error type for `valenx-thermo`.

use std::fmt;

/// Universal gas constant, J/(mol·K) (CODATA 2018, exact since the 2019 SI redefinition).
pub const GAS_CONSTANT: f64 = 8.314_462_618;

/// Errors produced by the thermodynamic models.
///
/// Every fallible public entry point validates its inputs and returns one of
/// these variants rather than panicking, so callers can recover from bad data
/// (negative temperatures, non-physical critical constants, non-converging
/// iterations) loudly and explicitly.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermoError {
    /// A quantity that must be strictly positive (temperature, pressure,
    /// critical constants, molar volume) was zero or negative.
    ///
    /// `name` identifies the offending quantity and `value` is what was passed.
    NonPositive {
        /// Human-readable name of the quantity (e.g. `"temperature"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },
    /// A parameter was outside its physically meaningful range.
    OutOfRange {
        /// Human-readable name of the quantity.
        name: &'static str,
        /// The offending value.
        value: f64,
        /// Description of the expected range (e.g. `"acentric factor in [-1, 2]"`).
        expected: &'static str,
    },
    /// An iterative solver (e.g. the saturation-pressure Newton loop) failed to
    /// converge within the allotted iterations.
    NotConverged {
        /// Name of the solver that failed.
        solver: &'static str,
        /// Number of iterations attempted.
        iterations: usize,
        /// Residual at the last iterate.
        residual: f64,
    },
    /// The requested phase root does not exist for the given state (e.g. asking
    /// for a vapor root well above the critical point where only one root is
    /// real).
    NoSuchRoot {
        /// Description of what was requested.
        what: &'static str,
    },
}

impl fmt::Display for ThermoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ThermoError::NonPositive { name, value } => {
                write!(f, "{name} must be strictly positive, got {value}")
            }
            ThermoError::OutOfRange {
                name,
                value,
                expected,
            } => write!(f, "{name} = {value} out of range ({expected})"),
            ThermoError::NotConverged {
                solver,
                iterations,
                residual,
            } => write!(
                f,
                "{solver} did not converge in {iterations} iterations (residual {residual:e})"
            ),
            ThermoError::NoSuchRoot { what } => {
                write!(f, "no real root for {what} at the requested state")
            }
        }
    }
}

impl std::error::Error for ThermoError {}

/// Convenience result alias used throughout the crate.
pub type Result<T> = std::result::Result<T, ThermoError>;

/// Validate that `value` is finite and strictly positive.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ThermoError::NonPositive { name, value })
    }
}
