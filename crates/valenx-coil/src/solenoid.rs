//! Long-solenoid inductor: geometry, inductance, energy, reactance and
//! the series time constant.
//!
//! ## Model
//!
//! For an ideal long air-cored solenoid of `N` turns wound uniformly
//! over an axial length `l` (metres) around a cross-sectional area `A`
//! (square metres), the self-inductance is
//!
//! ```text
//! L = mu0 * mu_r * N^2 * A / l
//! ```
//!
//! where `mu0` is the [vacuum permeability](VACUUM_PERMEABILITY) and
//! `mu_r` is the relative permeability of the core (1 for air / vacuum).
//! This is the classic textbook result obtained from Ampère's law under
//! the assumption that the field is uniform inside the coil and zero
//! outside, i.e. that `l` is large compared with the coil diameter so
//! end (fringing) effects are negligible.
//!
//! The downstream relations are exact consequences of `L`:
//!
//! ```text
//! E   = 0.5 * L * I^2          (energy stored in the field, joules)
//! X_L = 2 * pi * f * L         (inductive reactance, ohms)
//! tau = L / R                  (series-RL time constant, seconds)
//! i(t) energising = (V/R)(1 - exp(-t/tau))   (rising current, amperes)
//! i(t) decaying   = I0 * exp(-t/tau)         (source-free decay)
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. The long-solenoid formula
//! systematically *over*-estimates `L` for short, fat coils because it
//! ignores the Nagaoka end-correction factor; it also ignores the
//! proximity / skin effect, core saturation, hysteresis and
//! frequency-dependent core loss. Treat the numbers as first-order
//! textbook estimates, not production magnetics design values.

use crate::error::{require_positive, CoilError};
use serde::{Deserialize, Serialize};

/// Vacuum magnetic permeability `mu0`, in henries per metre (H/m).
///
/// Value: `4 * pi * 1e-7 ≈ 1.256_637_062e-6 H/m`. This is the classical
/// (pre-2019-SI) exact definition `mu0 = 4 * pi * 1e-7`, which remains
/// the conventional value used in textbook closed-form inductance
/// formulas and is accurate to well within the modelling error of the
/// long-solenoid approximation itself.
pub const VACUUM_PERMEABILITY: f64 = 4.0 * std::f64::consts::PI * 1.0e-7;

/// A uniform-permeability long solenoid, validated at construction.
///
/// All geometry is SI: `turns` is dimensionless, `area` is in square
/// metres, `length` is in metres, and `relative_permeability` is the
/// dimensionless `mu_r` of the core material (1 for air / vacuum). Build
/// one with [`Solenoid::new`]; the constructor rejects non-physical
/// inputs with a [`CoilError`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Solenoid {
    /// Number of turns `N` (dimensionless, `>= 0`).
    pub turns: f64,
    /// Cross-sectional area `A` enclosed by each turn, in square metres
    /// (`> 0`).
    pub area: f64,
    /// Axial length `l` of the winding, in metres (`> 0`).
    pub length: f64,
    /// Relative permeability `mu_r` of the core (dimensionless, `>= 1`).
    pub relative_permeability: f64,
}

impl Solenoid {
    /// Construct and validate a solenoid.
    ///
    /// # Errors
    ///
    /// Returns a [`CoilError`] if `area` or `length` is not strictly
    /// positive and finite, if `turns` is negative or non-finite, or if
    /// `relative_permeability` is below 1 or non-finite.
    pub fn new(
        turns: f64,
        area: f64,
        length: f64,
        relative_permeability: f64,
    ) -> Result<Self, CoilError> {
        if !turns.is_finite() {
            return Err(CoilError::NotFinite {
                name: "N",
                value: turns,
            });
        }
        if turns < 0.0 {
            return Err(CoilError::NegativeTurns { value: turns });
        }
        let area = require_positive("A", area)?;
        let length = require_positive("l", length)?;
        if !relative_permeability.is_finite() {
            return Err(CoilError::NotFinite {
                name: "mu_r",
                value: relative_permeability,
            });
        }
        if relative_permeability < 1.0 {
            return Err(CoilError::BadPermeability {
                value: relative_permeability,
            });
        }
        Ok(Self {
            turns,
            area,
            length,
            relative_permeability,
        })
    }

    /// Convenience constructor for an air-cored solenoid (`mu_r = 1`).
    ///
    /// # Errors
    ///
    /// Same validation as [`Solenoid::new`], with
    /// `relative_permeability` fixed to 1.
    pub fn air_core(turns: f64, area: f64, length: f64) -> Result<Self, CoilError> {
        Self::new(turns, area, length, 1.0)
    }

    /// The self-inductance `L`, in henries (H).
    ///
    /// Evaluates `L = mu0 * mu_r * N^2 * A / l`. Because every field is
    /// validated at construction, this getter is total and never fails.
    ///
    /// The dependence is quadratic in `N` and linear in `A`, `mu_r` and
    /// `1 / l`: doubling the turn count quadruples `L`, doubling the
    /// area or permeability doubles it, and doubling the length halves
    /// it.
    pub fn inductance_henries(&self) -> f64 {
        VACUUM_PERMEABILITY * self.relative_permeability * self.turns * self.turns * self.area
            / self.length
    }

    /// The magnetic energy `E = 0.5 * L * I^2` stored in the field at a
    /// steady current `current` (amperes), in joules (J).
    ///
    /// Delegates to [`energy_joules`] with this solenoid's inductance.
    pub fn energy_joules(&self, current: f64) -> f64 {
        energy_joules(self.inductance_henries(), current)
    }

    /// The inductive reactance `X_L = 2 * pi * f * L` at frequency
    /// `frequency` (hertz), in ohms.
    ///
    /// # Errors
    ///
    /// Returns a [`CoilError`] if `frequency` is not strictly positive
    /// and finite (see [`reactance_ohms`]).
    pub fn reactance_ohms(&self, frequency: f64) -> Result<f64, CoilError> {
        reactance_ohms(self.inductance_henries(), frequency)
    }

    /// The series-RL time constant `tau = L / R` (seconds) for this
    /// solenoid in series with a resistance `resistance` (ohms).
    ///
    /// # Errors
    ///
    /// Returns a [`CoilError`] if `resistance` is not strictly positive
    /// and finite (see [`time_constant_seconds`]).
    pub fn time_constant_seconds(&self, resistance: f64) -> Result<f64, CoilError> {
        time_constant_seconds(self.inductance_henries(), resistance)
    }

    /// The (fractional) number of turns `N` needed to reach a target
    /// self-inductance for a given geometry and core — the coil-winding
    /// design inverse of [`Solenoid::inductance_henries`].
    ///
    /// Inverting `L = mu0 mu_r N^2 A / l` gives
    ///
    /// ```text
    /// N = sqrt( L * l / (mu0 * mu_r * A) ).
    /// ```
    ///
    /// The result is the exact real turn count; round it up to a whole turn
    /// in practice. Building a [`Solenoid::new`] at this `N` (same `A`, `l`,
    /// `mu_r`) reproduces the target inductance. A higher-permeability core
    /// or shorter winding needs fewer turns.
    ///
    /// # Errors
    ///
    /// Returns a [`CoilError`] if `target_henries`, `area` or `length` is not
    /// strictly positive and finite, or if `relative_permeability` is below 1
    /// or non-finite.
    pub fn turns_for_inductance(
        target_henries: f64,
        area: f64,
        length: f64,
        relative_permeability: f64,
    ) -> Result<f64, CoilError> {
        let l_target = require_positive("L", target_henries)?;
        let area = require_positive("A", area)?;
        let length = require_positive("l", length)?;
        if !relative_permeability.is_finite() {
            return Err(CoilError::NotFinite {
                name: "mu_r",
                value: relative_permeability,
            });
        }
        if relative_permeability < 1.0 {
            return Err(CoilError::BadPermeability {
                value: relative_permeability,
            });
        }
        Ok((l_target * length / (VACUUM_PERMEABILITY * relative_permeability * area)).sqrt())
    }
}

/// The magnetic energy `E = 0.5 * L * I^2` stored in an inductor, in
/// joules (J).
///
/// `inductance` is in henries and `current` in amperes. The result is
/// quadratic in the current and is symmetric in its sign (the energy of
/// `+I` equals the energy of `-I`).
pub fn energy_joules(inductance: f64, current: f64) -> f64 {
    0.5 * inductance * current * current
}

/// The inductive reactance `X_L = 2 * pi * f * L`, in ohms.
///
/// `inductance` is in henries and `frequency` in hertz. Reactance is the
/// imaginary part of the inductor's impedance and grows linearly with
/// both `frequency` and `inductance`.
///
/// # Errors
///
/// Returns [`CoilError::NonPositive`] / [`CoilError::NotFinite`] if
/// `frequency` is not a strictly positive, finite value.
pub fn reactance_ohms(inductance: f64, frequency: f64) -> Result<f64, CoilError> {
    let frequency = require_positive("f", frequency)?;
    Ok(2.0 * std::f64::consts::PI * frequency * inductance)
}

/// The series-RL time constant `tau = L / R`, in seconds.
///
/// `inductance` is in henries and `resistance` in ohms. After one time
/// constant the current in a switched series-RL circuit reaches
/// `1 - 1/e ≈ 63.2 %` of its final value.
///
/// # Errors
///
/// Returns [`CoilError::NonPositive`] / [`CoilError::NotFinite`] if
/// `resistance` is not a strictly positive, finite value.
pub fn time_constant_seconds(inductance: f64, resistance: f64) -> Result<f64, CoilError> {
    let resistance = require_positive("R", resistance)?;
    Ok(inductance / resistance)
}

/// Validate a non-negative, finite time. `t = 0` is admissible.
fn require_time(time: f64) -> Result<f64, CoilError> {
    if !time.is_finite() {
        return Err(CoilError::NotFinite {
            name: "t",
            value: time,
        });
    }
    if time < 0.0 {
        return Err(CoilError::NonPositive {
            name: "t",
            value: time,
        });
    }
    Ok(time)
}

/// The current in a switched series-RL circuit **energising** from a DC
/// source `v_source` (volts) through resistance `resistance` (ohms) with
/// the coil's `inductance` (henries), at time `time` (seconds):
///
/// ```text
/// i(t) = (V / R) * (1 - exp(-t * R / L)).
/// ```
///
/// The current rises from `0` toward the steady value `V/R` with the
/// time constant `tau = L/R`; after one time constant it has reached
/// `1 - 1/e ≈ 63.2 %` of the final value. This is the inductive dual of
/// the capacitor's RC charging voltage.
///
/// # Errors
///
/// Returns [`CoilError::NotFinite`] if `v_source` or `time` is non-finite,
/// [`CoilError::NonPositive`] if `time` is negative or if `resistance` /
/// `inductance` is not strictly positive.
pub fn rising_current(
    v_source: f64,
    resistance: f64,
    inductance: f64,
    time: f64,
) -> Result<f64, CoilError> {
    if !v_source.is_finite() {
        return Err(CoilError::NotFinite {
            name: "V",
            value: v_source,
        });
    }
    let r = require_positive("R", resistance)?;
    let l = require_positive("L", inductance)?;
    let t = require_time(time)?;
    Ok((v_source / r) * (1.0 - (-t * r / l).exp()))
}

/// The current in a source-free series-RL circuit **decaying** from an
/// initial current `initial_current` (amperes) through resistance
/// `resistance` (ohms) with `inductance` (henries), at time `time`
/// (seconds):
///
/// ```text
/// i(t) = I0 * exp(-t * R / L).
/// ```
///
/// The current falls from `I0` toward zero with the time constant
/// `tau = L/R`; after one time constant it has fallen to `1/e ≈ 36.8 %`
/// of `I0`. The dual of the capacitor's RC discharge.
///
/// # Errors
///
/// Returns [`CoilError::NotFinite`] if `initial_current` or `time` is
/// non-finite, [`CoilError::NonPositive`] if `time` is negative or if
/// `resistance` / `inductance` is not strictly positive.
pub fn decaying_current(
    initial_current: f64,
    resistance: f64,
    inductance: f64,
    time: f64,
) -> Result<f64, CoilError> {
    if !initial_current.is_finite() {
        return Err(CoilError::NotFinite {
            name: "I0",
            value: initial_current,
        });
    }
    let r = require_positive("R", resistance)?;
    let l = require_positive("L", inductance)?;
    let t = require_time(time)?;
    Ok(initial_current * (-t * r / l).exp())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    /// Absolute tolerance for floating-point comparisons.
    const EPS: f64 = 1.0e-12;

    /// Vacuum permeability matches `4 * pi * 1e-7` to machine precision.
    #[test]
    fn mu0_value() {
        let expected = 1.256_637_061_435_917_3e-6;
        assert!((VACUUM_PERMEABILITY - expected).abs() < 1.0e-18);
    }

    /// Ground-truth inductance: N=100, A=1e-4 m^2, l=0.1 m, air core.
    ///
    /// L = mu0 * 100^2 * 1e-4 / 0.1
    ///   = (4 pi e-7) * 10000 * 1e-4 / 0.1
    ///   = (4 pi e-7) * (1.0 / 0.1)
    ///   = (4 pi e-7) * 10
    ///   = 4 pi e-6 ≈ 1.256_637e-5 H.
    #[test]
    fn inductance_ground_truth() {
        let coil = Solenoid::air_core(100.0, 1.0e-4, 0.1).unwrap();
        let expected = 4.0 * std::f64::consts::PI * 1.0e-6;
        assert!((coil.inductance_henries() - expected).abs() < EPS * expected.max(1.0));
    }

    /// L scales as N^2: doubling the turns quadruples the inductance.
    #[test]
    fn inductance_scales_n_squared() {
        let base = Solenoid::air_core(50.0, 2.0e-4, 0.05).unwrap();
        let doubled = Solenoid::air_core(100.0, 2.0e-4, 0.05).unwrap();
        let ratio = doubled.inductance_henries() / base.inductance_henries();
        assert!((ratio - 4.0).abs() < 1.0e-9);
    }

    /// Doubling the turns explicitly gives 4x L (the headline check).
    #[test]
    fn doubling_turns_quadruples_inductance() {
        let n = 37.0;
        let a = 3.3e-4;
        let l = 0.12;
        let single = Solenoid::air_core(n, a, l).unwrap().inductance_henries();
        let double = Solenoid::air_core(2.0 * n, a, l)
            .unwrap()
            .inductance_henries();
        assert!((double - 4.0 * single).abs() < 1.0e-12);
    }

    /// L scales linearly with the cross-sectional area A.
    #[test]
    fn inductance_scales_linearly_with_area() {
        let base = Solenoid::air_core(80.0, 1.0e-4, 0.2).unwrap();
        let triple_area = Solenoid::air_core(80.0, 3.0e-4, 0.2).unwrap();
        let ratio = triple_area.inductance_henries() / base.inductance_henries();
        assert!((ratio - 3.0).abs() < 1.0e-9);
    }

    /// L scales inversely with the length l: 2x length -> 0.5x L.
    #[test]
    fn inductance_scales_inverse_length() {
        let short = Solenoid::air_core(80.0, 1.0e-4, 0.1).unwrap();
        let long = Solenoid::air_core(80.0, 1.0e-4, 0.2).unwrap();
        let ratio = long.inductance_henries() / short.inductance_henries();
        assert!((ratio - 0.5).abs() < 1.0e-9);
    }

    /// A high-permeability core multiplies L by mu_r relative to air.
    #[test]
    fn inductance_scales_with_permeability() {
        let air = Solenoid::air_core(120.0, 5.0e-5, 0.08).unwrap();
        let cored = Solenoid::new(120.0, 5.0e-5, 0.08, 200.0).unwrap();
        let ratio = cored.inductance_henries() / air.inductance_henries();
        assert!((ratio - 200.0).abs() < 1.0e-7);
    }

    /// Zero turns yields exactly zero inductance.
    #[test]
    fn zero_turns_zero_inductance() {
        let coil = Solenoid::air_core(0.0, 1.0e-4, 0.1).unwrap();
        assert!(coil.inductance_henries().abs() < EPS);
    }

    /// Energy formula E = 0.5 L I^2, checked against a hand value.
    ///
    /// L = 2 H, I = 3 A -> E = 0.5 * 2 * 9 = 9 J.
    #[test]
    fn energy_ground_truth() {
        let e = energy_joules(2.0, 3.0);
        assert!((e - 9.0).abs() < EPS);
    }

    /// Energy is quadratic in current: doubling I quadruples E.
    #[test]
    fn energy_quadratic_in_current() {
        let l = 1.5e-3;
        let e1 = energy_joules(l, 2.0);
        let e2 = energy_joules(l, 4.0);
        assert!((e2 / e1 - 4.0).abs() < 1.0e-9);
    }

    /// Energy is even in the current sign: E(+I) == E(-I).
    #[test]
    fn energy_symmetric_in_sign() {
        let l = 0.7;
        let plus = energy_joules(l, 5.0);
        let minus = energy_joules(l, -5.0);
        assert!((plus - minus).abs() < EPS);
    }

    /// The struct getter agrees with the free `energy_joules` function.
    #[test]
    fn struct_energy_matches_free_fn() {
        let coil = Solenoid::new(200.0, 4.0e-4, 0.15, 50.0).unwrap();
        let l = coil.inductance_henries();
        let i = 1.2;
        assert!((coil.energy_joules(i) - energy_joules(l, i)).abs() < EPS);
    }

    /// Reactance formula X_L = 2 pi f L, checked against a hand value.
    ///
    /// L = 1e-3 H, f = 1000 Hz -> X_L = 2 pi * 1000 * 1e-3 = 2 pi ohms.
    #[test]
    fn reactance_ground_truth() {
        let x = reactance_ohms(1.0e-3, 1000.0).unwrap();
        let expected = 2.0 * std::f64::consts::PI;
        assert!((x - expected).abs() < 1.0e-9);
    }

    /// Reactance scales linearly with frequency: 2x f -> 2x X_L.
    #[test]
    fn reactance_scales_linearly_with_frequency() {
        let l = 2.2e-6;
        let x1 = reactance_ohms(l, 1.0e6).unwrap();
        let x2 = reactance_ohms(l, 2.0e6).unwrap();
        assert!((x2 / x1 - 2.0).abs() < 1.0e-9);
    }

    /// Reactance scales linearly with inductance: 2x L -> 2x X_L.
    #[test]
    fn reactance_scales_linearly_with_inductance() {
        let f = 60.0;
        let x1 = reactance_ohms(1.0e-3, f).unwrap();
        let x2 = reactance_ohms(2.0e-3, f).unwrap();
        assert!((x2 / x1 - 2.0).abs() < 1.0e-9);
    }

    /// Non-positive frequency is rejected.
    #[test]
    fn reactance_rejects_nonpositive_frequency() {
        let err = reactance_ohms(1.0e-3, 0.0).unwrap_err();
        assert_eq!(err.code(), "coil.non-positive");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    /// Time-constant formula tau = L / R, checked against a hand value.
    ///
    /// L = 10 mH, R = 5 ohms -> tau = 0.01 / 5 = 0.002 s.
    #[test]
    fn time_constant_ground_truth() {
        let tau = time_constant_seconds(10.0e-3, 5.0).unwrap();
        assert!((tau - 0.002).abs() < EPS);
    }

    /// tau scales linearly with L and inversely with R.
    #[test]
    fn time_constant_scaling() {
        let tau_base = time_constant_seconds(1.0e-3, 10.0).unwrap();
        let tau_more_l = time_constant_seconds(2.0e-3, 10.0).unwrap();
        let tau_more_r = time_constant_seconds(1.0e-3, 20.0).unwrap();
        assert!((tau_more_l / tau_base - 2.0).abs() < 1.0e-9);
        assert!((tau_more_r / tau_base - 0.5).abs() < 1.0e-9);
    }

    /// Non-positive resistance is rejected.
    #[test]
    fn time_constant_rejects_nonpositive_resistance() {
        let err = time_constant_seconds(1.0e-3, -2.0).unwrap_err();
        assert!(matches!(err, CoilError::NonPositive { .. }));
    }

    /// Constructor rejects a negative turn count.
    #[test]
    fn new_rejects_negative_turns() {
        let err = Solenoid::air_core(-5.0, 1.0e-4, 0.1).unwrap_err();
        assert_eq!(err.code(), "coil.negative-turns");
    }

    /// Constructor rejects a non-positive length.
    #[test]
    fn new_rejects_nonpositive_length() {
        let err = Solenoid::air_core(10.0, 1.0e-4, 0.0).unwrap_err();
        assert!(matches!(err, CoilError::NonPositive { name: "l", .. }));
    }

    /// Constructor rejects a non-positive area.
    #[test]
    fn new_rejects_nonpositive_area() {
        let err = Solenoid::air_core(10.0, -1.0e-4, 0.1).unwrap_err();
        assert!(matches!(err, CoilError::NonPositive { name: "A", .. }));
    }

    /// Constructor rejects relative permeability below 1.
    #[test]
    fn new_rejects_low_permeability() {
        let err = Solenoid::new(10.0, 1.0e-4, 0.1, 0.5).unwrap_err();
        assert_eq!(err.code(), "coil.bad-permeability");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    /// Constructor rejects a non-finite turn count.
    #[test]
    fn new_rejects_nonfinite_turns() {
        let err = Solenoid::air_core(f64::NAN, 1.0e-4, 0.1).unwrap_err();
        assert!(matches!(err, CoilError::NotFinite { name: "N", .. }));
    }

    /// A round-trip through serde JSON preserves the solenoid.
    #[test]
    fn serde_round_trip() {
        let coil = Solenoid::new(150.0, 2.5e-4, 0.09, 12.0).unwrap();
        let json = serde_json::to_string(&coil).unwrap();
        let back: Solenoid = serde_json::from_str(&json).unwrap();
        assert_eq!(coil, back);
    }

    /// Rising RL current reaches 63.2% of V/R after one time constant, with
    /// the exponent using tau = L/R, and starts at 0 / approaches V/R.
    #[test]
    fn rising_current_one_tau_and_limits() {
        let (v, r, l) = (10.0, 5.0, 0.01); // tau = L/R = 0.002 s, V/R = 2 A
        let tau = time_constant_seconds(l, r).unwrap();
        let at_tau = rising_current(v, r, l, tau).unwrap();
        let expected = (v / r) * (1.0 - (-1.0f64).exp());
        assert!((at_tau - expected).abs() < EPS, "i(tau) = {at_tau}");
        // Starts from zero.
        assert!(rising_current(v, r, l, 0.0).unwrap().abs() < EPS);
        // Approaches the steady current V/R: at 6 tau it is within ~0.5 %
        // and strictly below the asymptote (beyond ~37 tau the exponential
        // underflows to exactly V/R in f64).
        let near = rising_current(v, r, l, 6.0 * tau).unwrap();
        assert!(near < v / r, "must approach from below");
        assert!((near - v / r).abs() < 0.01, "within ~0.5% at 6 tau");
    }

    /// Decaying RL current falls to 1/e of I0 after one time constant, and
    /// starts at I0.
    #[test]
    fn decaying_current_one_tau_and_start() {
        let (i0, r, l) = (2.0, 5.0, 0.01);
        let tau = time_constant_seconds(l, r).unwrap();
        let at_tau = decaying_current(i0, r, l, tau).unwrap();
        assert!(
            (at_tau - i0 * (-1.0f64).exp()).abs() < EPS,
            "i(tau) = {at_tau}"
        );
        assert!((decaying_current(i0, r, l, 0.0).unwrap() - i0).abs() < EPS);
    }

    /// The rising and decaying currents are complementary about the steady
    /// value: (V/R - rising) equals a decay from V/R.
    #[test]
    fn rise_and_decay_are_complementary() {
        let (v, r, l) = (12.0, 4.0, 0.02);
        let i_inf = v / r;
        for &t in &[0.0, 1.0e-3, 5.0e-3, 2.0e-2] {
            let rise = rising_current(v, r, l, t).unwrap();
            let decay = decaying_current(i_inf, r, l, t).unwrap();
            assert!(
                (i_inf - rise - decay).abs() < 1.0e-12,
                "t={t}: {rise} {decay}"
            );
        }
    }

    /// RL transient rejects bad resistance / inductance / time / source.
    #[test]
    fn rl_transient_rejects_bad_inputs() {
        assert!(rising_current(10.0, 0.0, 0.01, 0.001).is_err()); // R = 0
        assert!(rising_current(10.0, 5.0, 0.0, 0.001).is_err()); // L = 0
        assert!(rising_current(10.0, 5.0, 0.01, -1.0).is_err()); // t < 0
        assert!(rising_current(f64::NAN, 5.0, 0.01, 0.001).is_err()); // V NaN
        assert!(decaying_current(2.0, -5.0, 0.01, 0.001).is_err()); // R < 0
        assert!(decaying_current(f64::INFINITY, 5.0, 0.01, 0.001).is_err()); // I0 inf
    }

    /// turns_for_inductance inverts inductance_henries both ways: size N for
    /// a target L and recover L, and recover a coil's own N from its L.
    #[test]
    fn turns_for_inductance_inverts_inductance() {
        let (a, l, mu_r) = (2.5e-4, 0.09, 12.0);
        let target = 3.0e-4;
        let n = Solenoid::turns_for_inductance(target, a, l, mu_r).unwrap();
        let coil = Solenoid::new(n, a, l, mu_r).unwrap();
        assert!((coil.inductance_henries() - target).abs() < 1.0e-9 * target);
        // Forward: a coil's own turn count is recovered from its inductance.
        let known = Solenoid::new(137.0, a, l, mu_r).unwrap();
        let n_back =
            Solenoid::turns_for_inductance(known.inductance_henries(), a, l, mu_r).unwrap();
        assert!((n_back - 137.0).abs() < 1.0e-9, "n_back = {n_back}");
    }

    /// Ground truth: air core A=1e-4, l=0.1, target L = 4 pi e-6 -> N = 100
    /// (the inverse of the inductance ground-truth case).
    #[test]
    fn turns_for_inductance_ground_truth() {
        let target = 4.0 * std::f64::consts::PI * 1.0e-6;
        let n = Solenoid::turns_for_inductance(target, 1.0e-4, 0.1, 1.0).unwrap();
        assert!((n - 100.0).abs() < 1.0e-9, "n = {n}");
    }

    /// Turns grow as sqrt(L) (4x L -> 2x turns); a permeable core needs fewer.
    #[test]
    fn turns_for_inductance_scales_and_permeability() {
        let (a, l) = (1.0e-4, 0.1);
        let lo = Solenoid::turns_for_inductance(1.0e-5, a, l, 1.0).unwrap();
        let hi = Solenoid::turns_for_inductance(4.0e-5, a, l, 1.0).unwrap();
        assert!((hi / lo - 2.0).abs() < 1.0e-9, "lo={lo} hi={hi}");
        let cored = Solenoid::turns_for_inductance(1.0e-5, a, l, 100.0).unwrap();
        assert!(cored < lo, "cored {cored} should be below air {lo}");
    }

    #[test]
    fn turns_for_inductance_rejects_bad_inputs() {
        assert!(Solenoid::turns_for_inductance(0.0, 1.0e-4, 0.1, 1.0).is_err()); // L
        assert!(Solenoid::turns_for_inductance(1.0e-5, 0.0, 0.1, 1.0).is_err()); // A
        assert!(Solenoid::turns_for_inductance(1.0e-5, 1.0e-4, 0.0, 1.0).is_err()); // l
        assert!(Solenoid::turns_for_inductance(1.0e-5, 1.0e-4, 0.1, 0.5).is_err()); // mu_r<1
        assert!(Solenoid::turns_for_inductance(f64::NAN, 1.0e-4, 0.1, 1.0).is_err());
    }
}
