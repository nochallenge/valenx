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
//! [`NortonLaw::rate_at`] evaluates the resulting rate.
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
    fn rejects_non_finite_inputs() {
        assert!(norton_creep_rate(f64::NAN, 100.0, 5.0).is_err());
        assert!(norton_creep_rate(1.0, 100.0, f64::INFINITY).is_err());
        assert!(NortonLaw::with_arrhenius(1.0, 1.0, 0.0, 5.0).is_err());
        assert!(NortonLaw::with_arrhenius(1.0, -1.0, 800.0, 5.0).is_err());
    }
}
