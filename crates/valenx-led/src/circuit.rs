//! Series-resistor LED circuit models and their closed-form DC solution.
//!
//! Two configurations are provided:
//!
//! - [`LedCircuit`] — one LED in series with a current-limiting resistor.
//! - [`LedString`] — `n` identical LEDs in series (their forward voltages
//!   add) in series with one resistor.
//!
//! Both reduce to the same single-loop Kirchhoff voltage equation. With a
//! supply voltage `Vs`, a total LED forward voltage `Vf` (a single drop, or
//! the sum over a string), and a desired LED current `I`, the resistor that
//! sets that current is
//!
//! ```text
//! R = (Vs - Vf) / I
//! ```
//!
//! and the steady-state power split is
//!
//! ```text
//! P_led      = Vf * I
//! P_resistor = (Vs - Vf) * I
//! P_total    = Vs * I = P_led + P_resistor
//! ```
//!
//! The LED forward voltage is treated as a fixed parameter (the standard
//! constant-drop diode model). A real diode's `Vf` varies weakly with
//! current and temperature; see the crate-level "Honest scope" note.

use crate::error::{require_non_negative, require_positive, LedError};
use serde::{Deserialize, Serialize};

/// A single LED driven through a series current-limiting resistor.
///
/// The fields are the three independent inputs of the design problem: the
/// supply voltage, the LED forward voltage, and the target LED current. The
/// resistor value and all power figures are *derived* on demand from these,
/// so the stored state is always internally consistent.
///
/// Construct with [`LedCircuit::new`], which validates the inputs and
/// guarantees `supply_v > forward_v` (positive resistor headroom).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LedCircuit {
    /// Supply (source) voltage `Vs`, in volts. Strictly positive.
    pub supply_v: f64,
    /// LED forward voltage `Vf`, in volts. Non-negative and `< supply_v`.
    pub forward_v: f64,
    /// Target LED current `I`, in amperes. Strictly positive.
    pub current_a: f64,
}

impl LedCircuit {
    /// Build a validated single-LED circuit.
    ///
    /// # Parameters
    ///
    /// - `supply_v`: supply voltage `Vs` (volts), must be finite and `> 0`.
    /// - `forward_v`: LED forward voltage `Vf` (volts), must be finite and
    ///   `>= 0`.
    /// - `current_a`: target LED current `I` (amperes), must be finite and
    ///   `> 0`.
    ///
    /// # Errors
    ///
    /// Returns [`LedError::NotFinite`], [`LedError::NonPositive`], or
    /// [`LedError::Negative`] for an out-of-domain individual parameter, and
    /// [`LedError::InsufficientHeadroom`] when `supply_v <= forward_v` (no
    /// positive voltage is left for the resistor, so no operating point
    /// exists).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_led::circuit::LedCircuit;
    ///
    /// // 5 V supply, 2 V red LED, 20 mA target.
    /// let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
    /// // R = (5 - 2) / 0.020 = 150 ohm.
    /// assert!((c.resistor_ohm() - 150.0).abs() < 1e-9);
    /// ```
    pub fn new(supply_v: f64, forward_v: f64, current_a: f64) -> Result<Self, LedError> {
        let supply_v = require_positive("supply_v", supply_v)?;
        let forward_v = require_non_negative("forward_v", forward_v)?;
        let current_a = require_positive("current_a", current_a)?;

        let headroom = supply_v - forward_v;
        if headroom <= 0.0 {
            return Err(LedError::InsufficientHeadroom {
                supply_v,
                forward_v,
                headroom,
            });
        }

        Ok(Self {
            supply_v,
            forward_v,
            current_a,
        })
    }

    /// Build a circuit from a **chosen series resistor** instead of a
    /// target current — the inverse of the design equation,
    /// `I = (Vs - Vf) / R`.
    ///
    /// Use this when you have already picked a resistor (e.g. a standard
    /// E-series value) and want the current it actually sets. It
    /// round-trips with [`LedCircuit::resistor_ohm`]: constructing from a
    /// resistor and reading the resistor back returns the same value (and
    /// likewise current ↔ resistor).
    ///
    /// # Errors
    ///
    /// Returns the same per-parameter domain errors as
    /// [`LedCircuit::new`] for `supply_v` / `forward_v`,
    /// [`LedError::NonPositive`] / [`LedError::NotFinite`] for
    /// `resistor_ohm`, and [`LedError::InsufficientHeadroom`] when
    /// `supply_v <= forward_v`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_led::circuit::LedCircuit;
    ///
    /// // 5 V supply, 2 V LED, a 150 ohm resistor -> 20 mA.
    /// let c = LedCircuit::from_resistor(5.0, 2.0, 150.0).unwrap();
    /// assert!((c.current_a - 0.020).abs() < 1e-9);
    /// ```
    pub fn from_resistor(
        supply_v: f64,
        forward_v: f64,
        resistor_ohm: f64,
    ) -> Result<Self, LedError> {
        let supply_v = require_positive("supply_v", supply_v)?;
        let forward_v = require_non_negative("forward_v", forward_v)?;
        let resistor_ohm = require_positive("resistor_ohm", resistor_ohm)?;
        let headroom = supply_v - forward_v;
        if headroom <= 0.0 {
            return Err(LedError::InsufficientHeadroom {
                supply_v,
                forward_v,
                headroom,
            });
        }
        Self::new(supply_v, forward_v, headroom / resistor_ohm)
    }

    /// Total forward-voltage drop across the LED portion of the loop, in
    /// volts.
    ///
    /// For a single LED this is simply its forward voltage; the method
    /// exists so [`LedCircuit`] and [`LedString`] share a common accessor.
    pub fn total_forward_v(&self) -> f64 {
        self.forward_v
    }

    /// Voltage developed across the series resistor, in volts.
    ///
    /// This is the Kirchhoff headroom `Vs - Vf`, guaranteed strictly
    /// positive by the constructor.
    pub fn resistor_voltage(&self) -> f64 {
        self.supply_v - self.total_forward_v()
    }

    /// Required series resistance, in ohms: `R = (Vs - Vf) / I`.
    ///
    /// This is the resistor that makes the loop carry exactly `current_a`.
    pub fn resistor_ohm(&self) -> f64 {
        self.resistor_voltage() / self.current_a
    }

    /// Power dissipated in the LED, in watts: `P_led = Vf * I`.
    pub fn led_power_w(&self) -> f64 {
        self.total_forward_v() * self.current_a
    }

    /// Power dissipated in the series resistor, in watts:
    /// `P_resistor = (Vs - Vf) * I`.
    pub fn resistor_power_w(&self) -> f64 {
        self.resistor_voltage() * self.current_a
    }

    /// Total power drawn from the supply, in watts: `P_total = Vs * I`.
    ///
    /// By conservation this equals [`Self::led_power_w`] plus
    /// [`Self::resistor_power_w`] (to within floating-point rounding).
    pub fn total_power_w(&self) -> f64 {
        self.supply_v * self.current_a
    }

    /// Fraction of supply power delivered to the LED (dimensionless, `0..1`).
    ///
    /// Equal to `Vf / Vs`. The remainder is burned as heat in the resistor;
    /// this ratio is the headline efficiency limitation of the simple
    /// series-resistor drive.
    pub fn led_efficiency(&self) -> f64 {
        self.led_power_w() / self.total_power_w()
    }
}

/// `n` identical LEDs in series, driven through one current-limiting
/// resistor.
///
/// In a series string the same current flows through every LED, and their
/// forward voltages add: the total drop is `count * forward_v_each`. The
/// resistor and power expressions are then identical to the single-LED case
/// with `Vf` replaced by that sum.
///
/// Construct with [`LedString::new`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LedString {
    /// Number of LEDs in the series string. At least one.
    pub count: usize,
    /// Forward voltage of one LED `Vf`, in volts. Non-negative.
    pub forward_v_each: f64,
    /// Supply (source) voltage `Vs`, in volts. Strictly positive.
    pub supply_v: f64,
    /// Target string current `I`, in amperes. Strictly positive.
    pub current_a: f64,
}

impl LedString {
    /// Build a validated series LED string.
    ///
    /// # Parameters
    ///
    /// - `count`: number of LEDs in series, must be `>= 1`.
    /// - `forward_v_each`: forward voltage of a single LED `Vf` (volts),
    ///   finite and `>= 0`.
    /// - `supply_v`: supply voltage `Vs` (volts), finite and `> 0`.
    /// - `current_a`: target current `I` (amperes), finite and `> 0`.
    ///
    /// # Errors
    ///
    /// Returns [`LedError::EmptyString`] for `count == 0`; the per-parameter
    /// domain errors for non-finite or out-of-range voltages and current;
    /// and [`LedError::InsufficientHeadroom`] when the supply does not exceed
    /// the summed forward voltage `count * forward_v_each`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_led::circuit::LedString;
    ///
    /// // Three 3.2 V white LEDs (9.6 V total) on a 12 V supply at 20 mA.
    /// let s = LedString::new(3, 3.2, 12.0, 0.020).unwrap();
    /// assert!((s.total_forward_v() - 9.6).abs() < 1e-9);
    /// // R = (12 - 9.6) / 0.020 = 120 ohm.
    /// assert!((s.resistor_ohm() - 120.0).abs() < 1e-9);
    /// ```
    pub fn new(
        count: usize,
        forward_v_each: f64,
        supply_v: f64,
        current_a: f64,
    ) -> Result<Self, LedError> {
        if count == 0 {
            return Err(LedError::EmptyString { count });
        }
        let forward_v_each = require_non_negative("forward_v_each", forward_v_each)?;
        let supply_v = require_positive("supply_v", supply_v)?;
        let current_a = require_positive("current_a", current_a)?;

        let total_forward = forward_v_each * count as f64;
        let headroom = supply_v - total_forward;
        if headroom <= 0.0 {
            return Err(LedError::InsufficientHeadroom {
                supply_v,
                forward_v: total_forward,
                headroom,
            });
        }

        Ok(Self {
            count,
            forward_v_each,
            supply_v,
            current_a,
        })
    }

    /// Build a string from a **chosen series resistor** instead of a
    /// target current — the inverse of the design equation,
    /// `I = (Vs - n * Vf) / R`.
    ///
    /// The string analogue of [`LedCircuit::from_resistor`]: given the
    /// resistor you have, it returns the current the string will carry.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`LedString::new`] for `count` /
    /// `forward_v_each` / `supply_v`, a domain error for `resistor_ohm`,
    /// and [`LedError::InsufficientHeadroom`] when the supply does not
    /// exceed the summed forward voltage `count * forward_v_each`.
    pub fn from_resistor(
        count: usize,
        forward_v_each: f64,
        supply_v: f64,
        resistor_ohm: f64,
    ) -> Result<Self, LedError> {
        if count == 0 {
            return Err(LedError::EmptyString { count });
        }
        let forward_v_each = require_non_negative("forward_v_each", forward_v_each)?;
        let supply_v = require_positive("supply_v", supply_v)?;
        let resistor_ohm = require_positive("resistor_ohm", resistor_ohm)?;
        let total_forward = forward_v_each * count as f64;
        let headroom = supply_v - total_forward;
        if headroom <= 0.0 {
            return Err(LedError::InsufficientHeadroom {
                supply_v,
                forward_v: total_forward,
                headroom,
            });
        }
        Self::new(count, forward_v_each, supply_v, headroom / resistor_ohm)
    }

    /// Summed forward voltage of the whole string, in volts:
    /// `count * forward_v_each`.
    pub fn total_forward_v(&self) -> f64 {
        self.forward_v_each * self.count as f64
    }

    /// Collapse the string into the equivalent single-loop [`LedCircuit`].
    ///
    /// The returned circuit has the same supply, current, and *total*
    /// forward voltage, so every derived quantity (resistor, powers,
    /// efficiency) matches the string. This is exact because a series string
    /// and a single LED with the summed `Vf` obey the identical loop
    /// equation.
    ///
    /// Never fails in practice: the string invariants already imply the
    /// circuit invariants, but the [`Result`] is propagated for totality.
    pub fn as_circuit(&self) -> Result<LedCircuit, LedError> {
        LedCircuit::new(self.supply_v, self.total_forward_v(), self.current_a)
    }

    /// Voltage developed across the series resistor, in volts:
    /// `Vs - count * forward_v_each`.
    pub fn resistor_voltage(&self) -> f64 {
        self.supply_v - self.total_forward_v()
    }

    /// Required series resistance, in ohms: `R = (Vs - n * Vf) / I`.
    pub fn resistor_ohm(&self) -> f64 {
        self.resistor_voltage() / self.current_a
    }

    /// Total LED power for the whole string, in watts:
    /// `P_leds = (n * Vf) * I`.
    pub fn led_power_w(&self) -> f64 {
        self.total_forward_v() * self.current_a
    }

    /// Power dissipated in the series resistor, in watts.
    pub fn resistor_power_w(&self) -> f64 {
        self.resistor_voltage() * self.current_a
    }

    /// Total power drawn from the supply, in watts: `P_total = Vs * I`.
    pub fn total_power_w(&self) -> f64 {
        self.supply_v * self.current_a
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::LedError;

    /// Absolute tolerance for floating-point comparisons throughout.
    const EPS: f64 = 1e-9;

    // -- ground-truth: R = (Vs - Vf) / I --------------------------------

    #[test]
    fn resistor_matches_closed_form_red_led() {
        // Canonical textbook example: 5 V rail, 2.0 V red LED, 20 mA.
        // R = (5 - 2) / 0.020 = 150 ohm exactly.
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        assert!((c.resistor_ohm() - 150.0).abs() < EPS);
        assert!((c.resistor_voltage() - 3.0).abs() < EPS);
    }

    #[test]
    fn resistor_matches_closed_form_blue_led() {
        // 12 V supply, 3.4 V blue LED, 10 mA.
        // R = (12 - 3.4) / 0.010 = 860 ohm.
        let c = LedCircuit::new(12.0, 3.4, 0.010).unwrap();
        assert!((c.resistor_ohm() - 860.0).abs() < EPS);
    }

    #[test]
    fn resistor_formula_holds_for_swept_inputs() {
        // Independently recompute R = (Vs - Vf)/I and compare.
        for &vs in &[3.3_f64, 5.0, 9.0, 12.0, 24.0] {
            for &vf in &[1.8_f64, 2.0, 3.2] {
                for &i in &[0.005_f64, 0.010, 0.020, 0.035] {
                    if vs <= vf {
                        continue;
                    }
                    let c = LedCircuit::new(vs, vf, i).unwrap();
                    let expected = (vs - vf) / i;
                    assert!(
                        (c.resistor_ohm() - expected).abs() < EPS,
                        "R mismatch at vs={vs} vf={vf} i={i}: \
                         got {} expected {expected}",
                        c.resistor_ohm()
                    );
                }
            }
        }
    }

    // -- ground-truth: LED power = Vf * I -------------------------------

    #[test]
    fn led_power_is_vf_times_i() {
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        // 2.0 V * 0.020 A = 0.040 W = 40 mW.
        assert!((c.led_power_w() - 0.040).abs() < EPS);
    }

    // -- ground-truth: resistor power = (Vs - Vf) * I -------------------

    #[test]
    fn resistor_power_is_headroom_times_i() {
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        // (5 - 2) V * 0.020 A = 0.060 W = 60 mW.
        assert!((c.resistor_power_w() - 0.060).abs() < EPS);
    }

    #[test]
    fn resistor_power_equals_i_squared_r() {
        // Cross-check the (Vs - Vf)*I form against the equivalent I^2 * R.
        let c = LedCircuit::new(9.0, 2.1, 0.015).unwrap();
        let p_vi = c.resistor_power_w();
        let p_i2r = c.current_a * c.current_a * c.resistor_ohm();
        assert!((p_vi - p_i2r).abs() < EPS);
    }

    // -- conservation: P_led + P_resistor = P_total = Vs * I ------------

    #[test]
    fn powers_sum_to_supply_power() {
        let c = LedCircuit::new(12.0, 3.4, 0.025).unwrap();
        let sum = c.led_power_w() + c.resistor_power_w();
        assert!((sum - c.total_power_w()).abs() < EPS);
        assert!((c.total_power_w() - 12.0 * 0.025).abs() < EPS);
    }

    #[test]
    fn efficiency_is_vf_over_vs() {
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        // Vf / Vs = 2 / 5 = 0.4.
        assert!((c.led_efficiency() - 0.4).abs() < EPS);
    }

    // -- ground-truth: series LEDs sum their forward voltages -----------

    #[test]
    fn series_string_sums_forward_voltage() {
        let s = LedString::new(3, 3.2, 12.0, 0.020).unwrap();
        // 3 * 3.2 = 9.6 V total drop.
        assert!((s.total_forward_v() - 9.6).abs() < EPS);
        // R = (12 - 9.6) / 0.020 = 120 ohm.
        assert!((s.resistor_ohm() - 120.0).abs() < EPS);
    }

    #[test]
    fn single_led_string_matches_single_circuit() {
        // A 1-LED string must be identical to the single-LED circuit.
        let s = LedString::new(1, 2.0, 5.0, 0.020).unwrap();
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        assert!((s.total_forward_v() - c.total_forward_v()).abs() < EPS);
        assert!((s.resistor_ohm() - c.resistor_ohm()).abs() < EPS);
        assert!((s.led_power_w() - c.led_power_w()).abs() < EPS);
        assert!((s.resistor_power_w() - c.resistor_power_w()).abs() < EPS);
    }

    #[test]
    fn string_as_circuit_is_consistent() {
        let s = LedString::new(4, 2.0, 12.0, 0.030).unwrap();
        let c = s.as_circuit().unwrap();
        // 4 * 2.0 = 8.0 V; R = (12 - 8)/0.030.
        assert!((c.total_forward_v() - 8.0).abs() < EPS);
        assert!((c.resistor_ohm() - s.resistor_ohm()).abs() < EPS);
        assert!((c.led_power_w() - s.led_power_w()).abs() < EPS);
        assert!((c.resistor_power_w() - s.resistor_power_w()).abs() < EPS);
        assert!((c.total_power_w() - s.total_power_w()).abs() < EPS);
    }

    #[test]
    fn string_led_power_is_total_vf_times_i() {
        let s = LedString::new(2, 3.0, 9.0, 0.020).unwrap();
        // (2 * 3.0) V * 0.020 A = 0.120 W.
        assert!((s.led_power_w() - 0.120).abs() < EPS);
    }

    #[test]
    fn string_powers_sum_to_supply_power() {
        let s = LedString::new(5, 1.8, 12.0, 0.018).unwrap();
        let sum = s.led_power_w() + s.resistor_power_w();
        assert!((sum - s.total_power_w()).abs() < EPS);
    }

    // -- ground-truth: current is set by the resistor -------------------
    //
    // Round-trip: design R for a target current, then recover that current
    // from the loop equation I = (Vs - Vf) / R. They must agree.

    #[test]
    fn current_is_set_by_resistor() {
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        let r = c.resistor_ohm();
        let recovered_i = (c.supply_v - c.total_forward_v()) / r;
        assert!((recovered_i - 0.020).abs() < EPS);
    }

    #[test]
    fn current_recovered_across_sweep() {
        for &vs in &[3.3_f64, 5.0, 12.0] {
            for &vf in &[1.8_f64, 2.0] {
                for &i in &[0.005_f64, 0.020, 0.040] {
                    let c = LedCircuit::new(vs, vf, i).unwrap();
                    let recovered = (vs - vf) / c.resistor_ohm();
                    assert!((recovered - i).abs() < EPS);
                }
            }
        }
    }

    #[test]
    fn larger_resistor_implies_smaller_current() {
        // Two designs on the same rail/LED: the higher target current needs
        // the smaller resistor (R is inversely proportional to I).
        let lo_i = LedCircuit::new(5.0, 2.0, 0.010).unwrap();
        let hi_i = LedCircuit::new(5.0, 2.0, 0.030).unwrap();
        assert!(hi_i.resistor_ohm() < lo_i.resistor_ohm());
        // Specifically R scales as 1/I: R(0.010)/R(0.030) == 3.
        assert!((lo_i.resistor_ohm() / hi_i.resistor_ohm() - 3.0).abs() < EPS);
    }

    // -- from_resistor: the resistor -> current inverse constructor -----

    #[test]
    fn from_resistor_matches_closed_form() {
        // 5 V, 2 V LED, 150 ohm -> I = (5 - 2)/150 = 20 mA.
        let c = LedCircuit::from_resistor(5.0, 2.0, 150.0).unwrap();
        assert!((c.current_a - 0.020).abs() < EPS, "I = {}", c.current_a);
        // ...and the resistor reads straight back.
        assert!((c.resistor_ohm() - 150.0).abs() < EPS);
    }

    #[test]
    fn from_resistor_inverts_new_across_sweep() {
        // Design R from a target current, then rebuild from that R: the
        // recovered circuit must match the original current and resistor.
        for &vs in &[3.3_f64, 5.0, 12.0, 24.0] {
            for &vf in &[1.8_f64, 2.0, 3.2] {
                for &i in &[0.005_f64, 0.020, 0.040] {
                    if vs <= vf {
                        continue;
                    }
                    let designed = LedCircuit::new(vs, vf, i).unwrap();
                    let rebuilt =
                        LedCircuit::from_resistor(vs, vf, designed.resistor_ohm()).unwrap();
                    assert!((rebuilt.current_a - i).abs() < EPS, "I at vs={vs} vf={vf}");
                    assert!((rebuilt.resistor_ohm() - designed.resistor_ohm()).abs() < EPS);
                }
            }
        }
    }

    #[test]
    fn string_from_resistor_matches_closed_form_and_round_trips() {
        // 3 x 3.2 V on 12 V, 120 ohm -> I = (12 - 9.6)/120 = 20 mA.
        let s = LedString::from_resistor(3, 3.2, 12.0, 120.0).unwrap();
        assert!((s.current_a - 0.020).abs() < EPS, "I = {}", s.current_a);
        assert!((s.resistor_ohm() - 120.0).abs() < EPS);
        // Round-trip against the design constructor.
        let designed = LedString::new(4, 2.0, 12.0, 0.030).unwrap();
        let rebuilt = LedString::from_resistor(4, 2.0, 12.0, designed.resistor_ohm()).unwrap();
        assert!((rebuilt.current_a - 0.030).abs() < EPS);
    }

    #[test]
    fn from_resistor_rejects_bad_inputs() {
        // No headroom (Vs == Vf).
        assert!(matches!(
            LedCircuit::from_resistor(2.0, 2.0, 150.0).unwrap_err(),
            LedError::InsufficientHeadroom { .. }
        ));
        // Non-positive / non-finite resistor.
        assert!(matches!(
            LedCircuit::from_resistor(5.0, 2.0, 0.0).unwrap_err(),
            LedError::NonPositive {
                name: "resistor_ohm",
                ..
            }
        ));
        assert!(matches!(
            LedCircuit::from_resistor(5.0, 2.0, f64::NAN).unwrap_err(),
            LedError::NotFinite { .. }
        ));
        // String: empty count and no headroom.
        assert!(matches!(
            LedString::from_resistor(0, 2.0, 5.0, 150.0).unwrap_err(),
            LedError::EmptyString { count: 0 }
        ));
        assert!(matches!(
            LedString::from_resistor(4, 3.2, 12.0, 100.0).unwrap_err(),
            LedError::InsufficientHeadroom { .. }
        ));
    }

    // -- boundary / validation ------------------------------------------

    #[test]
    fn zero_headroom_single_is_rejected() {
        // Vs == Vf: no voltage for the resistor.
        let err = LedCircuit::new(2.0, 2.0, 0.020).unwrap_err();
        assert!(matches!(err, LedError::InsufficientHeadroom { .. }));
        assert_eq!(err.code(), "led.insufficient-headroom");
    }

    #[test]
    fn supply_below_forward_is_rejected() {
        let err = LedCircuit::new(1.5, 2.0, 0.020).unwrap_err();
        assert!(matches!(
            err,
            LedError::InsufficientHeadroom { headroom, .. } if headroom < 0.0
        ));
    }

    #[test]
    fn string_zero_headroom_is_rejected() {
        // 4 * 3.2 = 12.8 V exceeds a 12 V supply.
        let err = LedString::new(4, 3.2, 12.0, 0.020).unwrap_err();
        assert!(matches!(err, LedError::InsufficientHeadroom { .. }));
    }

    #[test]
    fn empty_string_is_rejected() {
        let err = LedString::new(0, 2.0, 5.0, 0.020).unwrap_err();
        assert!(matches!(err, LedError::EmptyString { count: 0 }));
        assert_eq!(err.code(), "led.empty-string");
    }

    #[test]
    fn non_positive_current_is_rejected() {
        let err = LedCircuit::new(5.0, 2.0, 0.0).unwrap_err();
        assert!(matches!(
            err,
            LedError::NonPositive {
                name: "current_a",
                ..
            }
        ));

        let neg = LedCircuit::new(5.0, 2.0, -0.020).unwrap_err();
        assert!(matches!(neg, LedError::NonPositive { .. }));
    }

    #[test]
    fn non_positive_supply_is_rejected() {
        let err = LedCircuit::new(0.0, 2.0, 0.020).unwrap_err();
        assert!(matches!(
            err,
            LedError::NonPositive {
                name: "supply_v",
                ..
            }
        ));
    }

    #[test]
    fn negative_forward_voltage_is_rejected() {
        let err = LedCircuit::new(5.0, -2.0, 0.020).unwrap_err();
        assert!(matches!(
            err,
            LedError::Negative {
                name: "forward_v",
                ..
            }
        ));
    }

    #[test]
    fn zero_forward_voltage_is_allowed() {
        // Vf = 0 is degenerate but physically admissible (a wire/short);
        // the full supply appears across the resistor.
        let c = LedCircuit::new(5.0, 0.0, 0.020).unwrap();
        assert!((c.resistor_voltage() - 5.0).abs() < EPS);
        assert!(c.led_power_w().abs() < EPS);
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        assert!(matches!(
            LedCircuit::new(f64::NAN, 2.0, 0.020).unwrap_err(),
            LedError::NotFinite { .. }
        ));
        assert!(matches!(
            LedCircuit::new(5.0, f64::INFINITY, 0.020).unwrap_err(),
            LedError::NotFinite { .. }
        ));
        assert!(matches!(
            LedCircuit::new(5.0, 2.0, f64::NAN).unwrap_err(),
            LedError::NotFinite { .. }
        ));
    }

    #[test]
    fn serde_round_trips_circuit() {
        let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
        let json = serde_json::to_string(&c).unwrap();
        let back: LedCircuit = serde_json::from_str(&json).unwrap();
        assert_eq!(c, back);
    }
}
