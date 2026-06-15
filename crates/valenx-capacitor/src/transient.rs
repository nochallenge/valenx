//! Single-pole RC charge and discharge transients.
//!
//! ## Model
//!
//! A capacitor `C` charging through a series resistance `R` from a source
//! `V0`, starting from zero charge, follows the first-order response
//!
//! ```text
//! V(t) = V0 * (1 - exp(-t / (R * C)))
//! ```
//!
//! and a charged capacitor discharging through `R` decays as
//!
//! ```text
//! V(t) = V0 * exp(-t / (R * C))
//! ```
//!
//! The product [`tau = R * C`](time_constant) is the time constant, in
//! seconds. After one time constant the charging voltage has reached
//! `1 - 1/e ~= 63.2 %` of `V0`; after five it is within ~0.7 % of the
//! final value.
//!
//! ## Honest scope
//!
//! This is the ideal lumped first-order solution: a constant source, a
//! single linear resistor and a single linear capacitor with no leakage,
//! no ESR/ESL and no dielectric absorption ("soakage"). Real circuits show
//! additional slow tails and source/wiring impedance not captured here.

use crate::error::{CapacitorError, Result};

/// Fraction `1 - 1/e` of the final value reached after one time constant
/// while charging (the familiar "63 %" figure), to full `f64` precision.
pub const CHARGE_FRACTION_ONE_TAU: f64 = 0.632_120_558_828_557_7;

/// RC time constant `tau = R * C`, in seconds.
///
/// # Parameters
///
/// - `resistance_ohm` — series resistance `R`, in ohms. Must be `> 0`.
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
/// use valenx_capacitor::transient::time_constant;
///
/// // 1 kohm with 1 uF gives tau = 1 ms.
/// let tau = time_constant(1.0e3, 1.0e-6).unwrap();
/// assert!((tau - 1.0e-3).abs() < 1e-15);
/// ```
pub fn time_constant(resistance_ohm: f64, capacitance_f: f64) -> Result<f64> {
    let r = CapacitorError::require_positive(
        "resistance_ohm",
        resistance_ohm,
        "resistance must be positive",
    )?;
    let c = CapacitorError::require_positive(
        "capacitance_f",
        capacitance_f,
        "capacitance must be positive",
    )?;
    Ok(r * c)
}

/// Capacitor voltage while *charging* through a series resistance, in
/// volts.
///
/// Computes `V(t) = V0 * (1 - exp(-t / (R * C)))`, the rising response of
/// an initially-uncharged capacitor connected to a constant source `V0`.
///
/// # Parameters
///
/// - `v0_v` — source / final voltage `V0`, in volts. Any finite value is
///   accepted (it sets the asymptote and may be negative).
/// - `resistance_ohm` — series resistance `R`, in ohms. Must be `> 0`.
/// - `capacitance_f` — capacitance `C`, in farads. Must be `> 0`.
/// - `time_s` — elapsed time `t`, in seconds. Must be `>= 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `v0_v` is non-finite,
/// if `resistance_ohm` or `capacitance_f` is not strictly positive, or if
/// `time_s` is negative.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::transient::charging_voltage;
///
/// // At one time constant (tau = 1 ms here) the voltage is ~63.2 % of V0.
/// let v = charging_voltage(5.0, 1.0e3, 1.0e-6, 1.0e-3).unwrap();
/// assert!((v - 5.0 * 0.632_120_558_828_557_7).abs() < 1e-9);
/// ```
pub fn charging_voltage(
    v0_v: f64,
    resistance_ohm: f64,
    capacitance_f: f64,
    time_s: f64,
) -> Result<f64> {
    let (v0, t, tau) = validate_transient(v0_v, resistance_ohm, capacitance_f, time_s)?;
    Ok(v0 * (1.0 - (-t / tau).exp()))
}

/// Capacitor voltage while *discharging* through a series resistance, in
/// volts.
///
/// Computes `V(t) = V0 * exp(-t / (R * C))`, the decaying response of a
/// capacitor initially charged to `V0` and released into the resistor.
///
/// # Parameters
///
/// - `v0_v` — initial voltage `V0`, in volts. Any finite value is
///   accepted.
/// - `resistance_ohm` — series resistance `R`, in ohms. Must be `> 0`.
/// - `capacitance_f` — capacitance `C`, in farads. Must be `> 0`.
/// - `time_s` — elapsed time `t`, in seconds. Must be `>= 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `v0_v` is non-finite,
/// if `resistance_ohm` or `capacitance_f` is not strictly positive, or if
/// `time_s` is negative.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::transient::discharging_voltage;
///
/// // After one time constant the voltage has fallen to V0 / e (~36.8 %).
/// let v = discharging_voltage(5.0, 1.0e3, 1.0e-6, 1.0e-3).unwrap();
/// assert!((v - 5.0 * (-1.0f64).exp()).abs() < 1e-9);
/// ```
pub fn discharging_voltage(
    v0_v: f64,
    resistance_ohm: f64,
    capacitance_f: f64,
    time_s: f64,
) -> Result<f64> {
    let (v0, t, tau) = validate_transient(v0_v, resistance_ohm, capacitance_f, time_s)?;
    Ok(v0 * (-t / tau).exp())
}

/// Shared validation for the transient solvers. Returns the validated
/// `(v0, t, tau)` triple on success.
fn validate_transient(
    v0_v: f64,
    resistance_ohm: f64,
    capacitance_f: f64,
    time_s: f64,
) -> Result<(f64, f64, f64)> {
    if !v0_v.is_finite() {
        return Err(CapacitorError::InvalidParameter {
            name: "v0_v",
            value: v0_v,
            reason: "source voltage must be finite",
        });
    }
    let t =
        CapacitorError::require_non_negative("time_s", time_s, "elapsed time cannot be negative")?;
    let tau = time_constant(resistance_ohm, capacitance_f)?;
    Ok((v0_v, t, tau))
}
