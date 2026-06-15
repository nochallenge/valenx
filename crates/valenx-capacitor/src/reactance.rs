//! Capacitive reactance in a sinusoidal steady state.
//!
//! ## Model
//!
//! An ideal capacitor driven by a sinusoid of frequency `f` (Hz) presents
//! a frequency-dependent opposition to current called the capacitive
//! reactance:
//!
//! ```text
//! X_C = 1 / (2 * pi * f * C)
//! ```
//!
//! measured in ohms. Reactance falls as either frequency or capacitance
//! rises: at DC (`f -> 0`) it diverges (an ideal capacitor blocks DC), and
//! at high frequency it tends to zero (the capacitor approaches a short).
//!
//! ## Honest scope
//!
//! This treats the capacitor as a pure, lossless reactance. It ignores
//! equivalent series resistance (ESR), equivalent series inductance (ESL)
//! and the resulting self-resonant frequency, dielectric loss
//! (`tan delta`) and leakage. Real components only behave like this well
//! below their self-resonant frequency.

use crate::error::{CapacitorError, Result};

/// Mathematical constant `tau = 2 * pi`, the radians in one full turn.
const TAU: f64 = core::f64::consts::TAU;

/// Capacitive reactance `X_C` of an ideal capacitor, in ohms.
///
/// Computes `X_C = 1 / (2 * pi * f * C)`.
///
/// # Parameters
///
/// - `frequency_hz` — drive frequency `f`, in hertz. Must be `> 0` (the
///   reactance is undefined at DC, where it diverges).
/// - `capacitance_f` — capacitance `C`, in farads. Must be `> 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if either argument is not
/// strictly positive and finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::reactance::reactance;
///
/// // 1 uF at 1 kHz: X_C = 1 / (2 pi * 1000 * 1e-6) ~= 159.15 ohm.
/// let xc = reactance(1.0e3, 1.0e-6).unwrap();
/// assert!((xc - 159.154_943_091_9).abs() < 1e-6);
/// ```
pub fn reactance(frequency_hz: f64, capacitance_f: f64) -> Result<f64> {
    let f = CapacitorError::require_positive(
        "frequency_hz",
        frequency_hz,
        "frequency must be positive (reactance diverges at DC)",
    )?;
    let c = CapacitorError::require_positive(
        "capacitance_f",
        capacitance_f,
        "capacitance must be positive",
    )?;
    Ok(1.0 / (TAU * f * c))
}

/// Angular-frequency form of the capacitive reactance, in ohms.
///
/// Computes `X_C = 1 / (omega * C)` where `omega = 2 * pi * f` is the
/// angular frequency in radians per second. Equivalent to
/// [`reactance`] when `omega = 2 pi f`, but convenient when the caller
/// already works in `omega`.
///
/// # Parameters
///
/// - `omega_rad_s` — angular frequency `omega`, in rad/s. Must be `> 0`.
/// - `capacitance_f` — capacitance `C`, in farads. Must be `> 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if either argument is not
/// strictly positive and finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::reactance::reactance_omega;
///
/// let xc = reactance_omega(1000.0, 1.0e-6).unwrap();
/// assert!((xc - 1000.0).abs() < 1e-9);
/// ```
pub fn reactance_omega(omega_rad_s: f64, capacitance_f: f64) -> Result<f64> {
    let w = CapacitorError::require_positive(
        "omega_rad_s",
        omega_rad_s,
        "angular frequency must be positive",
    )?;
    let c = CapacitorError::require_positive(
        "capacitance_f",
        capacitance_f,
        "capacitance must be positive",
    )?;
    Ok(1.0 / (w * c))
}
