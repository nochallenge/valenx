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

/// Time (seconds) for a *charging* capacitor to first reach a target
/// voltage `target_voltage_v`, the inverse of [`charging_voltage`].
///
/// Inverting `V = V0 (1 - exp(-t/RC))` gives
/// `t = R C * ln( V0 / (V0 - V) ) = -R C * ln(1 - V/V0)`. The target must
/// lie between `0` and `V0` (same sign as `V0`, strictly below it, since
/// the asymptote is only approached as `t -> ∞`). At `V = 0` the time is
/// `0`; the design rule "5 time constants to settle" comes straight out of
/// this (reaching `99.33 %` needs `t = 5 R C`).
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `v0_v`/`target_voltage_v`
/// is non-finite, if `resistance_ohm` or `capacitance_f` is not strictly
/// positive, or if the target fraction `V/V0` is not in `[0, 1)`.
pub fn time_to_charge(
    v0_v: f64,
    resistance_ohm: f64,
    capacitance_f: f64,
    target_voltage_v: f64,
) -> Result<f64> {
    let tau = validate_inverse(v0_v, resistance_ohm, capacitance_f, target_voltage_v)?;
    let f = target_voltage_v / v0_v;
    if !(0.0..1.0).contains(&f) {
        return Err(CapacitorError::InvalidParameter {
            name: "target_voltage_v",
            value: target_voltage_v,
            reason: "charging target must be between 0 and V0 (V/V0 in [0, 1))",
        });
    }
    Ok(-tau * (1.0 - f).ln())
}

/// Time (seconds) for a *discharging* capacitor to fall from `V0` to a
/// target voltage `target_voltage_v`, the inverse of [`discharging_voltage`].
///
/// Inverting `V = V0 exp(-t/RC)` gives `t = R C * ln(V0 / V)`. The target
/// must lie between `0` (exclusive — fully discharged is only reached as
/// `t -> ∞`) and `V0` (inclusive, where `t = 0`).
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `v0_v`/`target_voltage_v`
/// is non-finite, if `resistance_ohm` or `capacitance_f` is not strictly
/// positive, or if the target fraction `V/V0` is not in `(0, 1]`.
pub fn time_to_discharge(
    v0_v: f64,
    resistance_ohm: f64,
    capacitance_f: f64,
    target_voltage_v: f64,
) -> Result<f64> {
    let tau = validate_inverse(v0_v, resistance_ohm, capacitance_f, target_voltage_v)?;
    let f = target_voltage_v / v0_v;
    if !f.is_finite() || f <= 0.0 || f > 1.0 {
        return Err(CapacitorError::InvalidParameter {
            name: "target_voltage_v",
            value: target_voltage_v,
            reason: "discharge target must be between 0 and V0 (V/V0 in (0, 1])",
        });
    }
    Ok(-tau * f.ln())
}

/// Shared validation for the inverse (time-to-voltage) solvers: checks the
/// voltages are finite and returns the validated time constant `tau`.
fn validate_inverse(
    v0_v: f64,
    resistance_ohm: f64,
    capacitance_f: f64,
    target_voltage_v: f64,
) -> Result<f64> {
    if !v0_v.is_finite() {
        return Err(CapacitorError::InvalidParameter {
            name: "v0_v",
            value: v0_v,
            reason: "source voltage must be finite",
        });
    }
    if !target_voltage_v.is_finite() {
        return Err(CapacitorError::InvalidParameter {
            name: "target_voltage_v",
            value: target_voltage_v,
            reason: "target voltage must be finite",
        });
    }
    time_constant(resistance_ohm, capacitance_f)
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

#[cfg(test)]
mod inverse_tests {
    use super::*;

    // 1 kohm, 1 uF -> tau = 1 ms.
    const R: f64 = 1.0e3;
    const C: f64 = 1.0e-6;
    const TAU: f64 = 1.0e-3;

    /// Reaching 63.2 % while charging takes exactly one time constant.
    #[test]
    fn charge_to_one_tau_fraction_takes_one_tau() {
        let v0 = 5.0;
        let t = time_to_charge(v0, R, C, v0 * CHARGE_FRACTION_ONE_TAU).unwrap();
        assert!((t - TAU).abs() < 1e-15, "t = {t}");
    }

    /// Falling to V0/e (~36.8 %) while discharging takes one time constant.
    #[test]
    fn discharge_to_one_over_e_takes_one_tau() {
        let v0 = 5.0;
        let t = time_to_discharge(v0, R, C, v0 * (-1.0f64).exp()).unwrap();
        assert!((t - TAU).abs() < 1e-15, "t = {t}");
    }

    /// A known value: charging to 50 % takes tau*ln(2).
    #[test]
    fn charge_to_half_takes_tau_ln2() {
        let t = time_to_charge(10.0, R, C, 5.0).unwrap();
        assert!((t - TAU * 2.0_f64.ln()).abs() < 1e-15, "t = {t}");
    }

    /// time_to_charge inverts charging_voltage over a sweep of times.
    #[test]
    fn charge_inverse_round_trips() {
        let v0 = 12.0;
        for &t in &[1.0e-4, 5.0e-4, 1.0e-3, 3.0e-3] {
            let v = charging_voltage(v0, R, C, t).unwrap();
            let back = time_to_charge(v0, R, C, v).unwrap();
            assert!((back - t).abs() < 1e-12, "t {t} -> v {v} -> {back}");
        }
    }

    /// time_to_discharge inverts discharging_voltage over a sweep of times.
    #[test]
    fn discharge_inverse_round_trips() {
        let v0 = 12.0;
        for &t in &[1.0e-4, 5.0e-4, 1.0e-3, 3.0e-3] {
            let v = discharging_voltage(v0, R, C, t).unwrap();
            let back = time_to_discharge(v0, R, C, v).unwrap();
            assert!((back - t).abs() < 1e-12, "t {t} -> v {v} -> {back}");
        }
    }

    /// The relations are sign-agnostic: a negative source charges by the
    /// same fraction in the same time.
    #[test]
    fn negative_source_uses_the_fraction() {
        let v0 = -5.0;
        let t = time_to_charge(v0, R, C, v0 * CHARGE_FRACTION_ONE_TAU).unwrap();
        assert!((t - TAU).abs() < 1e-15, "t = {t}");
    }

    #[test]
    fn rejects_out_of_range_targets_and_bad_rc() {
        // Charging cannot reach or exceed V0.
        assert!(time_to_charge(5.0, R, C, 5.0).is_err());
        assert!(time_to_charge(5.0, R, C, 6.0).is_err());
        assert!(time_to_charge(5.0, R, C, -1.0).is_err()); // wrong sign -> f < 0
                                                           // Discharge cannot reach 0 or exceed V0.
        assert!(time_to_discharge(5.0, R, C, 0.0).is_err());
        assert!(time_to_discharge(5.0, R, C, 6.0).is_err());
        // Bad R / C propagate from time_constant.
        assert!(time_to_charge(5.0, 0.0, C, 2.0).is_err());
        assert!(time_to_discharge(5.0, R, -1.0, 2.0).is_err());
    }
}
