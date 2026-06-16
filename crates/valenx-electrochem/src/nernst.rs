//! The Nernst equation for a half-reaction's reduction potential.
//!
//! For the reduction half-reaction
//!
//! ```text
//! Ox + n e-  ->  Red
//! ```
//!
//! the electrode (reduction) potential at arbitrary composition and
//! temperature is
//!
//! ```text
//! E = E0 - (R T / (n F)) ln Q
//! ```
//!
//! where `E0` is the standard reduction potential (volts), `R` the molar
//! gas constant, `T` the absolute temperature (kelvin), `n` the number of
//! electrons transferred, `F` the Faraday constant, and `Q` the
//! dimensionless reaction quotient written reduced-over-oxidised in the
//! same direction as the half-reaction.
//!
//! The relation inverts cleanly: given a measured potential `E`, the
//! reaction quotient is `Q = exp((E0 - E) n F / (R T))`
//! ([`HalfReaction::quotient_from_potential`] / [`nernst_quotient`]) —
//! the basis of a potentiometric (ion-selective or pH) sensor.
//!
//! ## Honest scope
//!
//! `Q` is built from ideal activities (concentrations / partial pressures
//! relative to the standard state); activity coefficients, liquid-junction
//! potentials, and electrode kinetics are out of scope. The result is the
//! thermodynamic equilibrium electrode potential, not a measured cell
//! voltage under load.

use serde::{Deserialize, Serialize};

use crate::constants::{FARADAY_C_PER_MOL, GAS_CONSTANT_J_PER_MOL_K, STANDARD_TEMPERATURE_K};
use crate::error::{require_finite, require_positive, ElectrochemError};

/// The thermal (Nernst) prefactor `R T / F`, in volts.
///
/// This is the voltage scale that multiplies `(1 / n) ln Q` in the Nernst
/// equation. At the standard temperature 298.15 K it is approximately
/// `0.025693 V`; multiplying by `ln 10 ~ 2.302585` recovers the familiar
/// base-10 slope of about `0.0592 V` per decade of `Q` (per unit `n`).
///
/// # Errors
///
/// Returns [`ElectrochemError::BadParameter`] if `temperature_k <= 0`, or
/// [`ElectrochemError::NonFinite`] if it is not finite.
pub fn thermal_voltage(temperature_k: f64) -> Result<f64, ElectrochemError> {
    let t = require_positive(
        temperature_k,
        "temperature_k",
        "absolute temperature must be strictly positive (kelvin)",
    )?;
    Ok(GAS_CONSTANT_J_PER_MOL_K * t / FARADAY_C_PER_MOL)
}

/// The thermal voltage `R T / F` at the standard temperature (298.15 K).
///
/// Provided as a convenience for the common "at 25 C" case. Equal to
/// `thermal_voltage(STANDARD_TEMPERATURE_K)` and approximately `0.025693 V`.
pub fn thermal_voltage_standard() -> f64 {
    GAS_CONSTANT_J_PER_MOL_K * STANDARD_TEMPERATURE_K / FARADAY_C_PER_MOL
}

/// The base-10 Nernst slope `(R T ln 10) / (n F)`, in volts per decade.
///
/// This is the change in electrode potential per ten-fold change in the
/// reaction quotient `Q`. For `n = 1` at 298.15 K it is approximately
/// `0.05916 V` (the textbook "59 mV per decade"); for general `n` it is
/// that value divided by `n`.
///
/// # Errors
///
/// Returns [`ElectrochemError::BadParameter`] if `temperature_k <= 0` or
/// `n <= 0`, or [`ElectrochemError::NonFinite`] for non-finite inputs.
pub fn nernst_slope_per_decade(n: f64, temperature_k: f64) -> Result<f64, ElectrochemError> {
    let n = require_positive(n, "n", "electrons transferred must be strictly positive")?;
    let vt = thermal_voltage(temperature_k)?;
    Ok(vt * std::f64::consts::LN_10 / n)
}

/// A reduction half-reaction parameterised for the Nernst equation.
///
/// Holds the standard reduction potential `E0`, the number of electrons
/// transferred `n`, and the absolute temperature `T`. Construct it through
/// [`HalfReaction::new`], which validates every field, then evaluate the
/// electrode potential at a given reaction quotient with
/// [`HalfReaction::potential`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct HalfReaction {
    /// Standard reduction potential `E0`, in volts (versus SHE).
    pub e_standard_v: f64,
    /// Number of electrons transferred, `n` (strictly positive).
    pub electrons: f64,
    /// Absolute temperature `T`, in kelvin (strictly positive).
    pub temperature_k: f64,
}

impl HalfReaction {
    /// Build and validate a [`HalfReaction`].
    ///
    /// `e_standard_v` is the standard reduction potential (any finite real
    /// value, including negative for easily-oxidised couples). `electrons`
    /// is `n` and must be strictly positive. `temperature_k` is the
    /// absolute temperature and must be strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::NonFinite`] if `e_standard_v` is not
    /// finite, or [`ElectrochemError::BadParameter`] if `electrons <= 0` or
    /// `temperature_k <= 0`.
    pub fn new(
        e_standard_v: f64,
        electrons: f64,
        temperature_k: f64,
    ) -> Result<Self, ElectrochemError> {
        let e_standard_v = require_finite(e_standard_v, "e_standard_v")?;
        let electrons = require_positive(
            electrons,
            "electrons",
            "electrons transferred must be strictly positive",
        )?;
        let temperature_k = require_positive(
            temperature_k,
            "temperature_k",
            "absolute temperature must be strictly positive (kelvin)",
        )?;
        Ok(Self {
            e_standard_v,
            electrons,
            temperature_k,
        })
    }

    /// Build a [`HalfReaction`] at the standard temperature (298.15 K).
    ///
    /// Equivalent to [`HalfReaction::new`] with
    /// `temperature_k = STANDARD_TEMPERATURE_K`.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::NonFinite`] if `e_standard_v` is not
    /// finite, or [`ElectrochemError::BadParameter`] if `electrons <= 0`.
    pub fn at_standard_temperature(
        e_standard_v: f64,
        electrons: f64,
    ) -> Result<Self, ElectrochemError> {
        Self::new(e_standard_v, electrons, STANDARD_TEMPERATURE_K)
    }

    /// The thermal voltage `R T / F` for this half-reaction's temperature.
    ///
    /// Always succeeds because the temperature was validated at
    /// construction time.
    pub fn thermal_voltage(&self) -> f64 {
        GAS_CONSTANT_J_PER_MOL_K * self.temperature_k / FARADAY_C_PER_MOL
    }

    /// The electrode (reduction) potential at reaction quotient `q`.
    ///
    /// Evaluates `E = E0 - (R T / (n F)) ln q`. At `q = 1` the logarithm
    /// vanishes and `E == E0`. With `n > 0`, increasing `q` (more product,
    /// less reactant) lowers `E`; decreasing `q` raises it.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::BadParameter`] if `q <= 0` (the natural
    /// log is undefined there), or [`ElectrochemError::NonFinite`] if `q` is
    /// not finite.
    pub fn potential(&self, q: f64) -> Result<f64, ElectrochemError> {
        let q = require_positive(
            q,
            "q",
            "reaction quotient must be strictly positive for ln Q",
        )?;
        let prefactor = self.thermal_voltage() / self.electrons;
        Ok(self.e_standard_v - prefactor * q.ln())
    }

    /// The reaction quotient `Q` that produces a measured electrode
    /// potential `E`, inverting [`HalfReaction::potential`]:
    ///
    /// ```text
    /// Q = exp((E0 - E) * n F / (R T)) = exp((E0 - E) / (R T / (n F)))
    /// ```
    ///
    /// This is how a potentiometric sensor (an ion-selective or pH
    /// electrode) reads a concentration ratio back from a voltage. At
    /// `E == E0` it returns `Q = 1`; a potential below `E0` gives `Q > 1`
    /// and one above gives `Q < 1`. The result is always strictly
    /// positive.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::NonFinite`] if `potential_v` is not
    /// finite.
    pub fn quotient_from_potential(&self, potential_v: f64) -> Result<f64, ElectrochemError> {
        let e = require_finite(potential_v, "potential_v")?;
        let prefactor = self.thermal_voltage() / self.electrons;
        Ok(((self.e_standard_v - e) / prefactor).exp())
    }
}

/// Free-function form of the Nernst equation.
///
/// Computes `E = E0 - (R T / (n F)) ln Q` directly from its parameters,
/// validating each. Equivalent to building a [`HalfReaction`] and calling
/// [`HalfReaction::potential`], but convenient for one-off evaluations.
///
/// # Errors
///
/// Returns [`ElectrochemError::NonFinite`] if `e_standard_v` is not finite,
/// or [`ElectrochemError::BadParameter`] if `electrons <= 0`,
/// `temperature_k <= 0`, or `q <= 0`.
pub fn nernst_potential(
    e_standard_v: f64,
    electrons: f64,
    temperature_k: f64,
    q: f64,
) -> Result<f64, ElectrochemError> {
    HalfReaction::new(e_standard_v, electrons, temperature_k)?.potential(q)
}

/// Free-function inverse Nernst: the reaction quotient `Q` for a measured
/// electrode potential `E`.
///
/// Inverts [`nernst_potential`]: `Q = exp((E0 - E) n F / (R T))`.
/// Equivalent to building a [`HalfReaction`] and calling
/// [`HalfReaction::quotient_from_potential`].
///
/// # Errors
///
/// Returns [`ElectrochemError::NonFinite`] if `e_standard_v` or
/// `potential_v` is not finite, or [`ElectrochemError::BadParameter`] if
/// `electrons <= 0` or `temperature_k <= 0`.
pub fn nernst_quotient(
    e_standard_v: f64,
    electrons: f64,
    temperature_k: f64,
    potential_v: f64,
) -> Result<f64, ElectrochemError> {
    HalfReaction::new(e_standard_v, electrons, temperature_k)?.quotient_from_potential(potential_v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    /// Loose tolerance for derived voltages (millivolt scale).
    const EPS: f64 = 1e-9;

    #[test]
    fn thermal_voltage_at_298k_is_about_25_7_mv() {
        // R T / F at 298.15 K. Textbook value ~0.02569 V.
        let vt = thermal_voltage(298.15).unwrap();
        assert!(
            (vt - 0.025_692_5).abs() < 1e-6,
            "RT/F at 298.15 K should be ~0.025693 V, got {vt}"
        );
        // The convenience wrapper must agree exactly with the explicit call.
        let vt_std = thermal_voltage_standard();
        assert!(
            (vt - vt_std).abs() < EPS,
            "standard-temperature wrapper disagreed: {vt} vs {vt_std}"
        );
    }

    #[test]
    fn nernst_slope_is_about_59_mv_per_decade_for_n1_at_25c() {
        // (R T ln10) / (n F) at 298.15 K, n = 1 -> ~0.05916 V/decade.
        let slope = nernst_slope_per_decade(1.0, 298.15).unwrap();
        assert!(
            (slope - 0.059_16).abs() < 1e-4,
            "n=1 slope should be ~0.0592 V/decade, got {slope}"
        );
        // For n = 2 the slope must be exactly half the n = 1 slope.
        let slope2 = nernst_slope_per_decade(2.0, 298.15).unwrap();
        assert!(
            (slope2 - slope / 2.0).abs() < EPS,
            "n=2 slope should be half of n=1: {slope2} vs {}",
            slope / 2.0
        );
    }

    #[test]
    fn potential_equals_e0_at_q_equals_one() {
        // ln(1) = 0, so E must equal E0 regardless of n, T, E0.
        for &(e0, n, t) in &[
            (0.34_f64, 2.0_f64, 298.15_f64),
            (-0.76, 2.0, 310.0),
            (1.23, 4.0, 273.15),
        ] {
            let hr = HalfReaction::new(e0, n, t).unwrap();
            let e = hr.potential(1.0).unwrap();
            assert!(
                (e - e0).abs() < EPS,
                "at Q=1 expected E==E0 ({e0}), got {e} for n={n}, T={t}"
            );
        }
    }

    #[test]
    fn potential_falls_as_q_rises_for_positive_n() {
        // With n > 0, dE/dQ < 0: increasing Q strictly lowers E.
        let hr = HalfReaction::at_standard_temperature(0.34, 2.0).unwrap();
        let e_low = hr.potential(0.1).unwrap();
        let e_mid = hr.potential(1.0).unwrap();
        let e_high = hr.potential(10.0).unwrap();
        assert!(
            e_low > e_mid && e_mid > e_high,
            "E must strictly fall as Q rises: {e_low} -> {e_mid} -> {e_high}"
        );
        // Symmetry: a decade below and a decade above Q=1 are mirror images
        // about E0, each offset by exactly the per-decade slope.
        let slope = nernst_slope_per_decade(2.0, STANDARD_TEMPERATURE_K).unwrap();
        assert!(
            (e_low - (e_mid + slope)).abs() < EPS,
            "Q=0.1 should sit one slope above E0: {e_low} vs {}",
            e_mid + slope
        );
        assert!(
            (e_high - (e_mid - slope)).abs() < EPS,
            "Q=10 should sit one slope below E0: {e_high} vs {}",
            e_mid - slope
        );
    }

    #[test]
    fn potential_matches_explicit_decade_offset() {
        // E(Q=10) - E0 should equal -(RT ln10)/(nF) exactly.
        let e0 = 0.799; // Ag+/Ag standard reduction potential, n = 1.
        let hr = HalfReaction::at_standard_temperature(e0, 1.0).unwrap();
        let e10 = hr.potential(10.0).unwrap();
        let expected = e0 - thermal_voltage_standard() * std::f64::consts::LN_10;
        assert!(
            (e10 - expected).abs() < EPS,
            "one-decade Nernst shift mismatch: {e10} vs {expected}"
        );
    }

    #[test]
    fn larger_n_compresses_the_swing() {
        // The prefactor RT/(nF) shrinks with n, so a fixed Q shifts E less.
        let q = 100.0;
        let small_n = HalfReaction::at_standard_temperature(0.0, 1.0)
            .unwrap()
            .potential(q)
            .unwrap();
        let large_n = HalfReaction::at_standard_temperature(0.0, 4.0)
            .unwrap()
            .potential(q)
            .unwrap();
        // Both negative (Q > 1), but the n=4 shift is 1/4 the n=1 shift.
        assert!(
            (small_n / large_n - 4.0).abs() < 1e-6,
            "n=1 shift should be 4x the n=4 shift: ratio {}",
            small_n / large_n
        );
    }

    #[test]
    fn free_function_agrees_with_struct() {
        let e_struct = HalfReaction::new(0.34, 2.0, 298.15)
            .unwrap()
            .potential(5.0)
            .unwrap();
        let e_free = nernst_potential(0.34, 2.0, 298.15, 5.0).unwrap();
        assert!(
            (e_struct - e_free).abs() < EPS,
            "free fn vs struct mismatch: {e_free} vs {e_struct}"
        );
    }

    #[test]
    fn rejects_non_positive_q() {
        let hr = HalfReaction::at_standard_temperature(0.0, 1.0).unwrap();
        let err = hr.potential(0.0).unwrap_err();
        assert_eq!(err.code(), "electrochem.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);
        assert!(hr.potential(-1.0).is_err());
    }

    #[test]
    fn rejects_bad_construction_inputs() {
        assert!(HalfReaction::new(0.0, 0.0, 298.15).is_err()); // n = 0
        assert!(HalfReaction::new(0.0, -2.0, 298.15).is_err()); // n < 0
        assert!(HalfReaction::new(0.0, 2.0, 0.0).is_err()); // T = 0
        assert!(HalfReaction::new(0.0, 2.0, -10.0).is_err()); // T < 0
        let nf = HalfReaction::new(f64::NAN, 2.0, 298.15).unwrap_err();
        assert_eq!(nf.code(), "electrochem.non_finite");
        assert_eq!(nf.category(), ErrorCategory::Numeric);
    }

    #[test]
    fn quotient_inverts_potential() {
        // Round-trip both directions at a non-trivial operating point.
        let hr = HalfReaction::new(0.34, 2.0, 298.15).unwrap();
        // Q -> E -> Q.
        let q0 = 37.0;
        let e = hr.potential(q0).unwrap();
        let q_back = hr.quotient_from_potential(e).unwrap();
        assert!((q_back - q0).abs() / q0 < 1e-9, "Q round-trip: {q_back}");
        // E -> Q -> E.
        let e_meas = 0.21;
        let q = hr.quotient_from_potential(e_meas).unwrap();
        let e_back = hr.potential(q).unwrap();
        assert!((e_back - e_meas).abs() < EPS, "E round-trip: {e_back}");
    }

    #[test]
    fn quotient_is_one_at_standard_potential() {
        let hr = HalfReaction::at_standard_temperature(0.799, 1.0).unwrap();
        let q = hr.quotient_from_potential(0.799).unwrap();
        assert!((q - 1.0).abs() < 1e-9, "Q at E0 should be 1, got {q}");
    }

    #[test]
    fn one_decade_offsets_give_quotient_ten_and_tenth() {
        // E = E0 - slope_per_decade => Q = 10 exactly; +slope => Q = 0.1.
        let hr = HalfReaction::at_standard_temperature(0.0, 1.0).unwrap();
        let slope = nernst_slope_per_decade(1.0, STANDARD_TEMPERATURE_K).unwrap();
        let q_lo = hr.quotient_from_potential(-slope).unwrap();
        assert!(
            (q_lo - 10.0).abs() < 1e-6,
            "one decade below E0 -> Q=10: {q_lo}"
        );
        let q_hi = hr.quotient_from_potential(slope).unwrap();
        assert!(
            (q_hi - 0.1).abs() < 1e-7,
            "one decade above E0 -> Q=0.1: {q_hi}"
        );
    }

    #[test]
    fn potential_below_e0_gives_quotient_above_one() {
        let hr = HalfReaction::at_standard_temperature(0.34, 2.0).unwrap();
        assert!(hr.quotient_from_potential(0.30).unwrap() > 1.0);
        assert!(hr.quotient_from_potential(0.40).unwrap() < 1.0);
    }

    #[test]
    fn quotient_free_function_agrees_with_struct() {
        let via_struct = HalfReaction::new(0.34, 2.0, 298.15)
            .unwrap()
            .quotient_from_potential(0.25)
            .unwrap();
        let via_free = nernst_quotient(0.34, 2.0, 298.15, 0.25).unwrap();
        assert!((via_struct - via_free).abs() < EPS);
    }

    #[test]
    fn quotient_rejects_non_finite_and_bad_construction() {
        let hr = HalfReaction::at_standard_temperature(0.0, 1.0).unwrap();
        assert!(hr.quotient_from_potential(f64::NAN).is_err());
        assert!(hr.quotient_from_potential(f64::INFINITY).is_err());
        // The free function also validates the construction inputs.
        assert!(nernst_quotient(0.0, 0.0, 298.15, 0.1).is_err()); // n = 0
        assert!(nernst_quotient(f64::NAN, 1.0, 298.15, 0.1).is_err()); // bad E0
    }

    #[test]
    fn temperature_scales_the_slope_linearly() {
        // RT/F is linear in T, so doubling (absolute) T doubles the slope.
        let s1 = nernst_slope_per_decade(1.0, 150.0).unwrap();
        let s2 = nernst_slope_per_decade(1.0, 300.0).unwrap();
        assert!(
            (s2 - 2.0 * s1).abs() < EPS,
            "slope should double with T: {s2} vs {}",
            2.0 * s1
        );
    }
}
