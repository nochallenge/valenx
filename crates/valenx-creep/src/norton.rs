//! Norton-Bailey secondary (steady-state) creep law.
//!
//! ## Model
//!
//! In the secondary stage of creep the strain rate is approximately
//! constant in time and follows a power law in the applied stress —
//! Norton's law (also called the Norton-Bailey or power law):
//!
//! ```text
//!   epsilon_dot = A * sigma^n
//! ```
//!
//! where `epsilon_dot` is the steady-state (minimum) creep strain rate
//! (per hour, say), `sigma` is the applied stress (MPa), `A` is a
//! temperature-dependent material constant and `n` is the stress
//! exponent (the slope of `log(epsilon_dot)` versus `log(sigma)`). For
//! metals `n` is commonly in the range `~3` (diffusion / climb-glide
//! creep) to `~8` and beyond (dislocation creep).
//!
//! Temperature enters through `A`. A common engineering form makes the
//! Arrhenius dependence explicit,
//!
//! ```text
//!   A(T) = A0 * exp(-Q / (R * T))
//! ```
//!
//! so that for fixed `A0`, `Q` and `sigma`, the creep rate rises with
//! temperature `T`. [`NortonLaw::with_arrhenius`] builds that form and
//! [`NortonLaw::rate_at`] evaluates the resulting rate. The inverse
//! [`norton_stress_for_rate`] / [`NortonLaw::stress_for_rate`] solves
//! the same law for the stress `sigma = (epsilon_dot / A)^(1/n)` that
//! produces a target (e.g. allowable) creep rate.
//!
//! ## Honest scope
//!
//! This is the textbook closed-form Norton power law. The constants
//! `A`, `n`, `A0` and `Q` are empirical fits that you must supply from
//! qualified data; the law describes only the secondary (steady-state)
//! stage and ignores primary transients and tertiary acceleration to
//! rupture. Research / educational grade only — not a substitute for a
//! validated creep-life assessment.

use crate::error::{require_finite, require_non_negative, require_positive, CreepError};
use serde::{Deserialize, Serialize};

/// Universal gas constant in J / (mol K). Used by the optional
/// Arrhenius temperature dependence of the Norton coefficient.
pub const GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314_462_618_153_24;

/// Evaluate Norton's secondary-creep law `epsilon_dot = A * sigma^n`.
///
/// `coefficient` is `A`, `stress` is `sigma` (in whatever stress unit
/// `A` was calibrated for, e.g. MPa) and `exponent` is the stress
/// exponent `n`. The returned strain rate is in the reciprocal of the
/// time unit baked into `A`.
///
/// # Errors
///
/// Returns [`CreepError`] if `coefficient` is non-finite or negative,
/// if `stress` is non-finite or negative, or if `exponent` is
/// non-finite.
///
/// # Examples
///
/// ```
/// use valenx_creep::norton::norton_creep_rate;
///
/// // A = 1e-12, sigma = 100, n = 5: rate = 1e-12 * 100^5 = 1e-2.
/// let rate = norton_creep_rate(1e-12, 100.0, 5.0).unwrap();
/// assert!((rate - 1e-2).abs() < 1e-12);
/// ```
pub fn norton_creep_rate(coefficient: f64, stress: f64, exponent: f64) -> Result<f64, CreepError> {
    let coefficient = require_non_negative("coefficient", coefficient)?;
    let stress = require_non_negative("stress", stress)?;
    let exponent = require_finite("exponent", exponent)?;
    Ok(coefficient * stress.powf(exponent))
}

/// Invert Norton's law for the stress that produces a target
/// steady-state creep rate: `sigma = (epsilon_dot / A)^(1/n)`.
///
/// This is the inverse of [`norton_creep_rate`]: given an allowable
/// secondary-creep rate, it returns the stress at which the law predicts
/// exactly that rate. The coefficient `A` and the stress exponent `n`
/// must be strictly positive (a zero coefficient or exponent has no
/// invertible stress dependence); `target_rate` must be non-negative,
/// and a zero rate maps to zero stress.
///
/// # Errors
///
/// Returns [`CreepError`] if `coefficient` or `exponent` is non-finite
/// or not strictly positive, or if `target_rate` is non-finite or
/// negative.
///
/// # Examples
///
/// ```
/// use valenx_creep::norton::{norton_creep_rate, norton_stress_for_rate};
///
/// // Round-trip: stress -> rate -> stress.
/// let rate = norton_creep_rate(1e-12, 100.0, 5.0).unwrap();
/// let sigma = norton_stress_for_rate(1e-12, rate, 5.0).unwrap();
/// assert!((sigma - 100.0).abs() < 1e-6);
/// ```
pub fn norton_stress_for_rate(
    coefficient: f64,
    target_rate: f64,
    exponent: f64,
) -> Result<f64, CreepError> {
    let coefficient = require_positive("coefficient", coefficient)?;
    let target_rate = require_non_negative("target_rate", target_rate)?;
    let exponent = require_positive("exponent", exponent)?;
    Ok((target_rate / coefficient).powf(1.0 / exponent))
}

/// A calibrated Norton-Bailey secondary-creep law.
///
/// Holds the stress exponent `n` and the coefficient `A` (optionally
/// resolved from an Arrhenius `A0 * exp(-Q / (R T))` form), and
/// evaluates the steady-state creep rate for an applied stress.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NortonLaw {
    /// Norton coefficient `A` (per the calibration's time unit and the
    /// stress unit raised to `-n`).
    coefficient: f64,
    /// Stress exponent `n`.
    exponent: f64,
}

impl NortonLaw {
    /// Build a Norton law from an explicit coefficient `A` and stress
    /// exponent `n`.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if `coefficient` is non-finite or
    /// negative, or if `exponent` is non-finite.
    pub fn new(coefficient: f64, exponent: f64) -> Result<Self, CreepError> {
        Ok(Self {
            coefficient: require_non_negative("coefficient", coefficient)?,
            exponent: require_finite("exponent", exponent)?,
        })
    }

    /// Build a Norton law whose coefficient is resolved from an
    /// Arrhenius temperature dependence `A = A0 * exp(-Q / (R * T))`.
    ///
    /// `a0` is the pre-exponential factor `A0`, `activation_energy` is
    /// the creep activation energy `Q` in J/mol, and `temperature_k` is
    /// the absolute temperature `T` in kelvin.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if `a0` is non-finite or negative, if
    /// `activation_energy` is non-finite or negative, if
    /// `temperature_k` is non-finite or not strictly positive, or if
    /// `exponent` is non-finite.
    pub fn with_arrhenius(
        a0: f64,
        activation_energy: f64,
        temperature_k: f64,
        exponent: f64,
    ) -> Result<Self, CreepError> {
        let a0 = require_non_negative("a0", a0)?;
        let activation_energy = require_non_negative("activation_energy", activation_energy)?;
        let temperature_k = require_positive("temperature_k", temperature_k)?;
        let coefficient =
            a0 * (-activation_energy / (GAS_CONSTANT_J_PER_MOL_K * temperature_k)).exp();
        Self::new(coefficient, exponent)
    }

    /// The Norton coefficient `A`.
    pub fn coefficient(&self) -> f64 {
        self.coefficient
    }

    /// The stress exponent `n`.
    pub fn exponent(&self) -> f64 {
        self.exponent
    }

    /// Steady-state creep rate `A * sigma^n` for the applied `stress`.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if `stress` is non-finite or negative.
    pub fn rate_at(&self, stress: f64) -> Result<f64, CreepError> {
        norton_creep_rate(self.coefficient, stress, self.exponent)
    }

    /// The applied stress that yields a given steady-state creep rate,
    /// inverting `epsilon_dot = A * sigma^n` for `sigma`.
    ///
    /// The inverse of [`NortonLaw::rate_at`]; delegates to
    /// [`norton_stress_for_rate`]. Requires the law's coefficient `A`
    /// and exponent `n` to be strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if the coefficient or exponent is not
    /// strictly positive, or if `target_rate` is non-finite or negative.
    pub fn stress_for_rate(&self, target_rate: f64) -> Result<f64, CreepError> {
        norton_stress_for_rate(self.coefficient, target_rate, self.exponent)
    }

    /// The **accumulated secondary-creep strain** after a time `time` under
    /// a constant `stress`: the time integral of the (constant) steady-state
    /// rate, `epsilon = epsilon_dot * time = A * sigma^n * time`.
    ///
    /// This is the secondary-stage strain only — it ignores the primary
    /// transient and the tertiary acceleration (the crate's stated scope),
    /// so it is the linear-in-time accumulation along the minimum-rate line.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if `stress` is non-finite or negative, or if
    /// `time` is non-finite or negative.
    pub fn accumulated_strain(&self, stress: f64, time: f64) -> Result<f64, CreepError> {
        let time = require_non_negative("time", time)?;
        let rate = self.rate_at(stress)?;
        Ok(rate * time)
    }

    /// The time for the accumulated secondary-creep strain to reach
    /// `strain_limit` under a constant `stress` — the inverse of
    /// [`NortonLaw::accumulated_strain`], `t = strain_limit / (A sigma^n)`.
    ///
    /// The classic creep design question ("time to 1 % strain"). A zero
    /// strain limit returns `0`; a zero creep rate (zero stress, or a zero
    /// coefficient) returns `+infinity`, the correct "never reached" answer
    /// for a positive strain limit.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if `stress` is non-finite or negative, or if
    /// `strain_limit` is non-finite or negative.
    pub fn time_to_strain(&self, stress: f64, strain_limit: f64) -> Result<f64, CreepError> {
        let strain_limit = require_non_negative("strain_limit", strain_limit)?;
        let rate = self.rate_at(stress)?;
        if strain_limit == 0.0 {
            return Ok(0.0);
        }
        // rate == 0 (zero stress / coefficient) => strain never reached;
        // strain_limit / 0.0 == +infinity, the correct infinite life.
        Ok(strain_limit / rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn rate_matches_defining_formula() {
        // Ground truth: A * sigma^n with clean integer powers.
        // 2.0 * 10^3 = 2000.
        let rate = norton_creep_rate(2.0, 10.0, 3.0).unwrap();
        assert!((rate - 2000.0).abs() < 1e-9, "got {rate}");
    }

    #[test]
    fn rate_scales_as_stress_to_the_n() {
        // Doubling the stress must multiply the rate by 2^n exactly.
        let a = 3.5e-9;
        let n = 4.0;
        let base = norton_creep_rate(a, 50.0, n).unwrap();
        let doubled = norton_creep_rate(a, 100.0, n).unwrap();
        let ratio = doubled / base;
        assert!(
            (ratio - 2f64.powf(n)).abs() < 1e-6,
            "ratio {ratio} should equal 2^{n}"
        );
    }

    #[test]
    fn rate_rises_with_stress() {
        let law = NortonLaw::new(1e-10, 5.0).unwrap();
        let lo = law.rate_at(80.0).unwrap();
        let hi = law.rate_at(120.0).unwrap();
        assert!(hi > lo, "higher stress should creep faster: {hi} vs {lo}");
    }

    #[test]
    fn linear_creep_when_exponent_is_one() {
        // n = 1 makes the rate exactly proportional to stress.
        let law = NortonLaw::new(0.25, 1.0).unwrap();
        let r = law.rate_at(40.0).unwrap();
        assert!((r - 10.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn zero_stress_gives_zero_rate_for_positive_exponent() {
        let r = norton_creep_rate(5.0, 0.0, 3.0).unwrap();
        assert!(r.abs() < EPS, "expected 0, got {r}");
    }

    #[test]
    fn arrhenius_rate_rises_with_temperature() {
        // Fixed A0, Q, stress and n: hotter must creep faster because
        // A(T) = A0 exp(-Q/RT) increases with T.
        let a0 = 1.0e6;
        let q = 250_000.0; // J/mol, a representative creep activation.
        let n = 5.0;
        let sigma = 100.0;
        let cool = NortonLaw::with_arrhenius(a0, q, 900.0, n)
            .unwrap()
            .rate_at(sigma)
            .unwrap();
        let hot = NortonLaw::with_arrhenius(a0, q, 1000.0, n)
            .unwrap()
            .rate_at(sigma)
            .unwrap();
        assert!(hot > cool, "hotter should creep faster: {hot} vs {cool}");
    }

    #[test]
    fn arrhenius_coefficient_matches_closed_form() {
        // Verify A(T) against a hand-evaluated exponential.
        let a0 = 2.0;
        let q = 100_000.0;
        let t = 800.0;
        let law = NortonLaw::with_arrhenius(a0, q, t, 4.0).unwrap();
        let expected = a0 * (-q / (GAS_CONSTANT_J_PER_MOL_K * t)).exp();
        let rel = (law.coefficient() - expected).abs() / expected;
        assert!(rel < 1e-12, "coeff {} vs {expected}", law.coefficient());
    }

    #[test]
    fn accessors_round_trip_inputs() {
        let law = NortonLaw::new(7.0e-8, 6.5).unwrap();
        assert!((law.coefficient() - 7.0e-8).abs() < EPS);
        assert!((law.exponent() - 6.5).abs() < EPS);
    }

    #[test]
    fn rejects_negative_coefficient_and_stress() {
        assert!(norton_creep_rate(-1.0, 100.0, 5.0).is_err());
        assert!(norton_creep_rate(1.0, -100.0, 5.0).is_err());
        assert!(NortonLaw::new(-1.0e-9, 5.0).is_err());
    }

    #[test]
    fn stress_for_rate_inverts_the_forward_law() {
        // Round-trip: stress -> rate -> stress recovers the original.
        let (a, n) = (1e-12, 5.0);
        let sigma0 = 100.0;
        let rate = norton_creep_rate(a, sigma0, n).unwrap();
        let sigma = norton_stress_for_rate(a, rate, n).unwrap();
        assert!((sigma - sigma0).abs() / sigma0 < 1e-9, "got {sigma}");
    }

    #[test]
    fn stress_for_rate_matches_closed_form() {
        // sigma = (rate/A)^(1/n). A=2, rate=2000, n=3 -> (1000)^(1/3) = 10.
        let sigma = norton_stress_for_rate(2.0, 2000.0, 3.0).unwrap();
        assert!((sigma - 10.0).abs() < 1e-9, "got {sigma}");
    }

    #[test]
    fn stress_for_rate_via_law_round_trips() {
        let law = NortonLaw::new(3.5e-9, 4.0).unwrap();
        let rate = law.rate_at(75.0).unwrap();
        let sigma = law.stress_for_rate(rate).unwrap();
        assert!((sigma - 75.0).abs() / 75.0 < 1e-9, "got {sigma}");
    }

    #[test]
    fn zero_rate_maps_to_zero_stress() {
        let sigma = norton_stress_for_rate(5.0, 0.0, 3.0).unwrap();
        assert!(sigma.abs() < EPS, "got {sigma}");
    }

    #[test]
    fn stress_for_rate_is_monotonic() {
        let law = NortonLaw::new(1e-10, 5.0).unwrap();
        let lo = law.stress_for_rate(1e-3).unwrap();
        let hi = law.stress_for_rate(1e-1).unwrap();
        assert!(
            hi > lo,
            "higher target rate needs higher stress: {hi} vs {lo}"
        );
    }

    #[test]
    fn stress_for_rate_rejects_bad_domain() {
        // A zero/negative coefficient or exponent is not invertible.
        assert!(norton_stress_for_rate(0.0, 1.0, 5.0).is_err());
        assert!(norton_stress_for_rate(1e-10, 1.0, 0.0).is_err());
        assert!(norton_stress_for_rate(1e-10, 1.0, -2.0).is_err());
        // Negative or non-finite target rate.
        assert!(norton_stress_for_rate(1e-10, -1.0, 5.0).is_err());
        assert!(norton_stress_for_rate(1e-10, f64::NAN, 5.0).is_err());
        // A law with a zero coefficient cannot be inverted.
        assert!(NortonLaw::new(0.0, 5.0)
            .unwrap()
            .stress_for_rate(1.0)
            .is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(norton_creep_rate(f64::NAN, 100.0, 5.0).is_err());
        assert!(norton_creep_rate(1.0, 100.0, f64::INFINITY).is_err());
        assert!(NortonLaw::with_arrhenius(1.0, 1.0, 0.0, 5.0).is_err());
        assert!(NortonLaw::with_arrhenius(1.0, -1.0, 800.0, 5.0).is_err());
    }

    #[test]
    fn accumulated_strain_is_rate_times_time() {
        let law = NortonLaw::new(2.0, 3.0).unwrap();
        let stress = 10.0;
        let rate = law.rate_at(stress).unwrap(); // 2 * 10^3 = 2000
        let eps = law.accumulated_strain(stress, 5.0).unwrap();
        assert!((eps - rate * 5.0).abs() < 1e-6 * (rate * 5.0), "eps {eps}");
        // Linear in time: doubling the time doubles the strain; t=0 -> 0.
        let eps2 = law.accumulated_strain(stress, 10.0).unwrap();
        assert!((eps2 - 2.0 * eps).abs() < 1e-6 * eps);
        assert!(law.accumulated_strain(stress, 0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn time_to_strain_inverts_accumulated_strain() {
        let law = NortonLaw::new(1e-10, 5.0).unwrap();
        let stress = 120.0;
        // time -> strain -> time round-trips.
        for &t in &[1.0_f64, 100.0, 1.0e4] {
            let eps = law.accumulated_strain(stress, t).unwrap();
            let back = law.time_to_strain(stress, eps).unwrap();
            assert!((back - t).abs() < 1e-6 * t, "t {t} -> eps {eps} -> {back}");
        }
        // strain -> time -> strain round-trips at a 1% strain limit.
        let eps_limit = 0.01;
        let t = law.time_to_strain(stress, eps_limit).unwrap();
        assert!((law.accumulated_strain(stress, t).unwrap() - eps_limit).abs() < 1e-12);
    }

    #[test]
    fn time_to_strain_hand_value_and_limits() {
        // rate = 2 * 10^3 = 2000 per unit time; 1% strain in 0.01/2000.
        let law = NortonLaw::new(2.0, 3.0).unwrap();
        let t = law.time_to_strain(10.0, 0.01).unwrap();
        assert!((t - 0.01 / 2000.0).abs() < 1e-15, "t {t}");
        // Zero strain limit -> zero time.
        assert!(law.time_to_strain(10.0, 0.0).unwrap().abs() < EPS);
        // Zero stress -> zero rate -> infinite life for a positive strain.
        assert!(law.time_to_strain(0.0, 0.01).unwrap().is_infinite());
    }

    #[test]
    fn strain_and_time_reject_bad_inputs() {
        let law = NortonLaw::new(1e-10, 5.0).unwrap();
        assert!(law.accumulated_strain(100.0, -1.0).is_err()); // time < 0
        assert!(law.accumulated_strain(-1.0, 1.0).is_err()); // stress < 0
        assert!(law.time_to_strain(100.0, -0.01).is_err()); // strain < 0
        assert!(law.time_to_strain(100.0, f64::NAN).is_err());
    }
}
