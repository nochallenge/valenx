//! Van't Hoff osmotic pressure — the colligative core.
//!
//! For a dilute solution the osmotic pressure is
//!
//! ```text
//! Pi = i * c * R * T
//! ```
//!
//! where
//!
//! - `i` is the **van't Hoff factor** (the number of dissolved
//!   particles each formula unit dissociates into — `1` for a non-
//!   electrolyte like glucose, `2` for NaCl, `3` for CaCl2, ...),
//! - `c` is the **molar concentration** of the solute formula unit
//!   in mol per litre (`mol/L` = `mol/dm^3`),
//! - `R` is the **ideal-gas constant**, and
//! - `T` is the **absolute temperature** in kelvin.
//!
//! The product `i * c` is the **osmolarity** (osmol/L) — the
//! concentration of osmotically active particles, which is what the
//! pressure actually responds to.
//!
//! ## Units
//!
//! Osmotic pressure follows the ideal-gas analogy `Pi*V = n*R*T`, so the
//! unit of `Pi` is fixed by the unit chosen for `R`:
//!
//! - [`R_L_ATM`] (`L*atm/(mol*K)`) with `c` in `mol/L` gives `Pi` in
//!   **atmospheres**.
//! - [`R_SI`] (`J/(mol*K)` = `Pa*m^3/(mol*K)`) requires `c` in
//!   `mol/m^3` (note: `1 mol/L = 1000 mol/m^3`) and yields `Pi` in
//!   **pascals**.
//!
//! [`Solution::osmotic_pressure_atm`] and
//! [`Solution::osmotic_pressure_pa`] handle the bookkeeping for you.

use crate::error::OsmosisError;
use serde::{Deserialize, Serialize};

/// Ideal-gas constant in `L*atm/(mol*K)` (CODATA-derived).
///
/// Use with molar concentration in `mol/L` to obtain osmotic pressure
/// in atmospheres.
pub const R_L_ATM: f64 = 0.082_057_366_080_96;

/// Ideal-gas constant in SI units, `J/(mol*K)` = `Pa*m^3/(mol*K)`
/// (CODATA 2018 exact value).
///
/// Use with molar concentration in `mol/m^3` to obtain osmotic pressure
/// in pascals.
pub const R_SI: f64 = 8.314_462_618;

/// Conversion factor: `1 atm` expressed in pascals (exact, by
/// definition of the standard atmosphere).
pub const ATM_IN_PA: f64 = 101_325.0;

/// Absolute zero offset: `0 C` in kelvin (exact, by definition).
pub const CELSIUS_ZERO_K: f64 = 273.15;

/// A dilute aqueous solution characterized by what matters for the
/// colligative (van't Hoff) osmotic pressure: the solute concentration,
/// its dissociation factor, and the temperature.
///
/// Construct via [`Solution::new`], which validates every field.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Solution {
    /// Molar concentration of the solute **formula unit**, `mol/L`.
    concentration_mol_per_l: f64,
    /// Van't Hoff dissociation factor `i` (dimensionless, `>= 1`).
    vant_hoff_i: f64,
    /// Absolute temperature, kelvin (`> 0`).
    temperature_k: f64,
}

impl Solution {
    /// Build a validated [`Solution`].
    ///
    /// # Parameters
    ///
    /// - `concentration_mol_per_l` — solute formula-unit molarity in
    ///   `mol/L`; must be finite and `>= 0` (a zero-concentration
    ///   solution is pure solvent and is allowed — it simply has zero
    ///   osmotic pressure).
    /// - `vant_hoff_i` — dissociation factor; must be finite and
    ///   `>= 1.0` (each formula unit yields at least one particle).
    /// - `temperature_k` — absolute temperature; must be finite and
    ///   strictly `> 0` (Kelvin is an absolute scale).
    ///
    /// # Errors
    ///
    /// Returns [`OsmosisError::InvalidParameter`] naming the first field
    /// that fails its constraint.
    pub fn new(
        concentration_mol_per_l: f64,
        vant_hoff_i: f64,
        temperature_k: f64,
    ) -> Result<Self, OsmosisError> {
        if !concentration_mol_per_l.is_finite() || concentration_mol_per_l < 0.0 {
            return Err(OsmosisError::invalid(
                "concentration_mol_per_l",
                concentration_mol_per_l,
                "must be finite and >= 0",
            ));
        }
        if !vant_hoff_i.is_finite() || vant_hoff_i < 1.0 {
            return Err(OsmosisError::invalid(
                "vant_hoff_i",
                vant_hoff_i,
                "must be finite and >= 1",
            ));
        }
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(OsmosisError::invalid(
                "temperature_k",
                temperature_k,
                "must be finite and > 0",
            ));
        }
        Ok(Self {
            concentration_mol_per_l,
            vant_hoff_i,
            temperature_k,
        })
    }

    /// Build a [`Solution`] taking the temperature in degrees Celsius.
    ///
    /// Convenience wrapper around [`Solution::new`]; converts via
    /// [`CELSIUS_ZERO_K`] and applies the same validation (so e.g.
    /// `-300 C`, which is below absolute zero, is rejected).
    ///
    /// # Errors
    ///
    /// Propagates [`OsmosisError::InvalidParameter`] from
    /// [`Solution::new`].
    pub fn new_celsius(
        concentration_mol_per_l: f64,
        vant_hoff_i: f64,
        temperature_c: f64,
    ) -> Result<Self, OsmosisError> {
        Self::new(
            concentration_mol_per_l,
            vant_hoff_i,
            temperature_c + CELSIUS_ZERO_K,
        )
    }

    /// Molar concentration of the solute formula unit, `mol/L`.
    pub fn concentration_mol_per_l(&self) -> f64 {
        self.concentration_mol_per_l
    }

    /// Van't Hoff dissociation factor `i`.
    pub fn vant_hoff_i(&self) -> f64 {
        self.vant_hoff_i
    }

    /// Absolute temperature, kelvin.
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Osmolarity `i * c`, in osmol/L — the concentration of
    /// osmotically active particles.
    ///
    /// This is the quantity tonicity comparisons (see
    /// [`crate::tonicity`]) are made on.
    pub fn osmolarity_osmol_per_l(&self) -> f64 {
        self.vant_hoff_i * self.concentration_mol_per_l
    }

    /// Osmotic pressure in **atmospheres**, `Pi = i * c * R * T` using
    /// [`R_L_ATM`].
    pub fn osmotic_pressure_atm(&self) -> f64 {
        self.osmolarity_osmol_per_l() * R_L_ATM * self.temperature_k
    }

    /// Osmotic pressure in **pascals**.
    ///
    /// Computed in SI directly: the osmolarity is converted from
    /// `osmol/L` to `osmol/m^3` (`* 1000`) and multiplied by
    /// [`R_SI`] and `T`. (Equivalently `osmotic_pressure_atm * ATM_IN_PA`
    /// to within floating-point rounding.)
    pub fn osmotic_pressure_pa(&self) -> f64 {
        let osmolarity_per_m3 = self.osmolarity_osmol_per_l() * 1000.0;
        osmolarity_per_m3 * R_SI * self.temperature_k
    }
}

/// Free-function form of the van't Hoff law in atmospheres.
///
/// `Pi = i * c * R_L_ATM * T`, with `c` in `mol/L` and `T` in kelvin.
///
/// # Errors
///
/// Validates the inputs exactly as [`Solution::new`] does and returns
/// [`OsmosisError::InvalidParameter`] on the first violation.
pub fn osmotic_pressure_atm(
    vant_hoff_i: f64,
    concentration_mol_per_l: f64,
    temperature_k: f64,
) -> Result<f64, OsmosisError> {
    Ok(Solution::new(concentration_mol_per_l, vant_hoff_i, temperature_k)?.osmotic_pressure_atm())
}

/// Invert the van't Hoff law: the osmolarity (`i * c`, osmol/L) that
/// produces a given osmotic pressure at temperature `T`.
///
/// `i * c = Pi / (R_L_ATM * T)`.
///
/// # Errors
///
/// Returns [`OsmosisError::InvalidParameter`] if `pressure_atm` is not
/// finite or negative, or if `temperature_k` is not strictly positive
/// and finite.
pub fn osmolarity_from_pressure_atm(
    pressure_atm: f64,
    temperature_k: f64,
) -> Result<f64, OsmosisError> {
    if !pressure_atm.is_finite() || pressure_atm < 0.0 {
        return Err(OsmosisError::invalid(
            "pressure_atm",
            pressure_atm,
            "must be finite and >= 0",
        ));
    }
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(OsmosisError::invalid(
            "temperature_k",
            temperature_k,
            "must be finite and > 0",
        ));
    }
    Ok(pressure_atm / (R_L_ATM * temperature_k))
}

/// Determine a solute's **molar mass** from a measured osmotic pressure
/// — the classic membrane-osmometry method for sizing a macromolecule
/// (polymer, protein) too large to weigh out mole-by-mole.
///
/// Writing the molar concentration as the mass concentration over the
/// molar mass, `c = rho / M`, in the van't Hoff law `Pi = i * c * R * T`
/// and solving for `M` gives
///
/// ```text
/// M = i * rho * R_L_ATM * T / Pi
/// ```
///
/// with `rho` the mass concentration in `g/L`, `Pi` in atmospheres, `T`
/// in kelvin, and the result in `g/mol`. For a non-dissociating
/// macromolecule pass `vant_hoff_i = 1`. This is the exact inverse of
/// the forward law solved for `M`: feeding the result back as
/// `c = rho / M` reproduces `Pi`.
///
/// # Errors
///
/// Returns [`OsmosisError::InvalidParameter`] if `pressure_atm` is not
/// finite and strictly positive (a molar mass cannot be inferred from a
/// zero or absent pressure), if `mass_concentration_g_per_l` is not
/// finite and strictly positive, if `vant_hoff_i` is not finite and
/// `>= 1`, or if `temperature_k` is not finite and strictly positive.
pub fn molar_mass_from_pressure_atm(
    pressure_atm: f64,
    mass_concentration_g_per_l: f64,
    vant_hoff_i: f64,
    temperature_k: f64,
) -> Result<f64, OsmosisError> {
    if !pressure_atm.is_finite() || pressure_atm <= 0.0 {
        return Err(OsmosisError::invalid(
            "pressure_atm",
            pressure_atm,
            "must be finite and > 0",
        ));
    }
    if !mass_concentration_g_per_l.is_finite() || mass_concentration_g_per_l <= 0.0 {
        return Err(OsmosisError::invalid(
            "mass_concentration_g_per_l",
            mass_concentration_g_per_l,
            "must be finite and > 0",
        ));
    }
    if !vant_hoff_i.is_finite() || vant_hoff_i < 1.0 {
        return Err(OsmosisError::invalid(
            "vant_hoff_i",
            vant_hoff_i,
            "must be finite and >= 1",
        ));
    }
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(OsmosisError::invalid(
            "temperature_k",
            temperature_k,
            "must be finite and > 0",
        ));
    }
    Ok(vant_hoff_i * mass_concentration_g_per_l * R_L_ATM * temperature_k / pressure_atm)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_negative_concentration() {
        let err = Solution::new(-0.1, 1.0, 300.0).unwrap_err();
        assert_eq!(err.code(), "osmosis.invalid_parameter");
        match err {
            OsmosisError::InvalidParameter { name, .. } => {
                assert_eq!(name, "concentration_mol_per_l")
            }
        }
    }

    #[test]
    fn rejects_i_below_one() {
        let err = Solution::new(0.1, 0.5, 300.0).unwrap_err();
        match err {
            OsmosisError::InvalidParameter { name, .. } => assert_eq!(name, "vant_hoff_i"),
        }
    }

    #[test]
    fn rejects_nonpositive_temperature() {
        assert!(Solution::new(0.1, 1.0, 0.0).is_err());
        assert!(Solution::new(0.1, 1.0, -10.0).is_err());
    }

    #[test]
    fn celsius_below_absolute_zero_is_rejected() {
        // -300 C is below absolute zero (~-273.15 C).
        assert!(Solution::new_celsius(0.1, 1.0, -300.0).is_err());
        // 0 C is fine and equals 273.15 K.
        let s = Solution::new_celsius(0.1, 1.0, 0.0).unwrap();
        assert!((s.temperature_k() - 273.15).abs() < EPS);
    }

    #[test]
    fn zero_concentration_gives_zero_pressure() {
        let s = Solution::new(0.0, 2.0, 310.15).unwrap();
        assert!(s.osmotic_pressure_atm().abs() < EPS);
        assert!(s.osmotic_pressure_pa().abs() < EPS);
        assert!(s.osmolarity_osmol_per_l().abs() < EPS);
    }

    // ---- Ground-truth numeric checks -------------------------------

    #[test]
    fn one_molar_ideal_solute_at_stp_matches_known_value() {
        // A 1 mol/L non-electrolyte at 273.15 K: Pi = 1*1*R*T.
        // R*T = 0.0820574 * 273.15 = 22.414... L*atm/mol — the molar
        // volume of an ideal gas at STP. So Pi ~= 22.414 atm.
        let s = Solution::new(1.0, 1.0, CELSIUS_ZERO_K).unwrap();
        let pi = s.osmotic_pressure_atm();
        assert!((pi - 22.414).abs() < 1e-3, "expected ~22.414 atm, got {pi}");
    }

    #[test]
    fn physiological_saline_is_about_seven_point_seven_atm() {
        // 0.9% NaCl ~= 0.154 mol/L, i ~= 2 (full dissociation),
        // body temperature 37 C = 310.15 K.
        // Pi = 2 * 0.154 * 0.0820574 * 310.15 ~= 7.84 atm.
        // Physiology textbooks quote plasma osmotic pressure ~= 7.7 atm
        // (i for NaCl is closer to 1.85 in vivo, which lowers it a hair);
        // we check the ideal i=2 prediction lands in the right band.
        let s = Solution::new(0.154, 2.0, 310.15).unwrap();
        let pi = s.osmotic_pressure_atm();
        assert!((7.0..8.5).contains(&pi), "expected ~7.7-7.8 atm, got {pi}");
    }

    #[test]
    fn nacl_gives_double_the_pressure_of_glucose_at_equal_molarity() {
        // i=2 (NaCl) vs i=1 (glucose) at the same c and T.
        let nacl = Solution::new(0.1, 2.0, 298.15).unwrap();
        let glucose = Solution::new(0.1, 1.0, 298.15).unwrap();
        let ratio = nacl.osmotic_pressure_atm() / glucose.osmotic_pressure_atm();
        assert!((ratio - 2.0).abs() < EPS, "ratio = {ratio}");
    }

    // ---- The required VALIDATE properties --------------------------

    #[test]
    fn pressure_is_linear_in_concentration() {
        let t = 298.15;
        let base = Solution::new(0.05, 1.0, t).unwrap().osmotic_pressure_atm();
        // Doubling and tripling c must scale Pi by exactly 2 and 3.
        let double = Solution::new(0.10, 1.0, t).unwrap().osmotic_pressure_atm();
        let triple = Solution::new(0.15, 1.0, t).unwrap().osmotic_pressure_atm();
        assert!((double - 2.0 * base).abs() < EPS, "double = {double}");
        assert!((triple - 3.0 * base).abs() < EPS, "triple = {triple}");
    }

    #[test]
    fn pressure_is_linear_in_absolute_temperature() {
        let c = 0.2;
        let p200 = Solution::new(c, 1.0, 200.0).unwrap().osmotic_pressure_atm();
        let p400 = Solution::new(c, 1.0, 400.0).unwrap().osmotic_pressure_atm();
        // Pi(400 K) must be exactly twice Pi(200 K) — proportional to T
        // on the ABSOLUTE (kelvin) scale, the key teaching point.
        assert!(
            (p400 - 2.0 * p200).abs() < EPS,
            "p400 = {p400}, p200 = {p200}"
        );
    }

    #[test]
    fn pressure_extrapolates_to_zero_at_zero_kelvin() {
        // Linearity through the origin: as T -> 0+, Pi -> 0.
        let c = 0.2;
        let tiny = Solution::new(c, 1.0, 1e-6).unwrap().osmotic_pressure_atm();
        assert!(tiny.abs() < 1e-6, "Pi near 0 K should be ~0, got {tiny}");
    }

    #[test]
    fn atm_and_pa_results_are_consistent() {
        let s = Solution::new(0.3, 2.0, 305.0).unwrap();
        let from_atm = s.osmotic_pressure_atm() * ATM_IN_PA;
        let direct_pa = s.osmotic_pressure_pa();
        // The two unit paths must agree to a tight relative tolerance.
        let rel = ((from_atm - direct_pa) / direct_pa).abs();
        assert!(rel < 1e-4, "atm-path {from_atm} vs SI-path {direct_pa}");
    }

    #[test]
    fn free_function_matches_method() {
        let pi_fn = osmotic_pressure_atm(2.0, 0.15, 310.0).unwrap();
        let pi_method = Solution::new(0.15, 2.0, 310.0)
            .unwrap()
            .osmotic_pressure_atm();
        assert!((pi_fn - pi_method).abs() < EPS);
    }

    #[test]
    fn inversion_round_trips() {
        let s = Solution::new(0.12, 2.0, 300.0).unwrap();
        let pi = s.osmotic_pressure_atm();
        let recovered = osmolarity_from_pressure_atm(pi, 300.0).unwrap();
        // i*c = 2 * 0.12 = 0.24 osmol/L.
        assert!((recovered - 0.24).abs() < EPS, "recovered = {recovered}");
        assert!((recovered - s.osmolarity_osmol_per_l()).abs() < EPS);
    }

    #[test]
    fn inversion_rejects_bad_inputs() {
        assert!(osmolarity_from_pressure_atm(-1.0, 300.0).is_err());
        assert!(osmolarity_from_pressure_atm(1.0, 0.0).is_err());
        assert!(osmolarity_from_pressure_atm(f64::NAN, 300.0).is_err());
    }

    // ---- Molar mass from osmotic pressure (membrane osmometry) ------

    #[test]
    fn molar_mass_round_trips_through_the_forward_law() {
        // Pick a molar mass, make a solution at a known mass
        // concentration, compute its osmotic pressure, then recover the
        // mass. The last case exercises the dissociation factor i > 1.
        for &(m_true, rho, i) in &[
            (50_000.0, 10.0, 1.0),
            (1_800.0, 25.0, 1.0),
            (12_000.0, 8.0, 2.0),
        ] {
            let t = 298.15;
            let c = rho / m_true; // molar concentration, mol/L
            let pi = Solution::new(c, i, t).unwrap().osmotic_pressure_atm();
            let m = molar_mass_from_pressure_atm(pi, rho, i, t).unwrap();
            assert!(
                (m - m_true).abs() / m_true < 1e-9,
                "recovered M = {m}, expected {m_true}"
            );
        }
    }

    #[test]
    fn molar_mass_matches_hand_formula() {
        // M = i * rho * R * T / Pi, recomputed independently.
        let m = molar_mass_from_pressure_atm(0.5, 20.0, 1.0, 300.0).unwrap();
        let expected = 1.0 * 20.0 * R_L_ATM * 300.0 / 0.5;
        assert!((m - expected).abs() < EPS, "M = {m} vs {expected}");
    }

    #[test]
    fn molar_mass_equals_mass_concentration_over_molarity() {
        // Derivable identity: M = i * rho / (i*c), where the osmolarity
        // i*c is exactly what osmolarity_from_pressure_atm recovers.
        let (rho, i, t, pi) = (12.0, 1.0, 310.0, 0.03);
        let m = molar_mass_from_pressure_atm(pi, rho, i, t).unwrap();
        let osmolarity = osmolarity_from_pressure_atm(pi, t).unwrap();
        let via_osmolarity = i * rho / osmolarity;
        assert!(
            (m - via_osmolarity).abs() / m < 1e-12,
            "M {m} vs {via_osmolarity}"
        );
    }

    #[test]
    fn molar_mass_scales_linearly_with_mass_concentration() {
        // At fixed Pi, T, i: M is proportional to the mass concentration.
        let base = molar_mass_from_pressure_atm(0.2, 5.0, 1.0, 300.0).unwrap();
        let double = molar_mass_from_pressure_atm(0.2, 10.0, 1.0, 300.0).unwrap();
        assert!(
            (double - 2.0 * base).abs() / base < 1e-12,
            "double = {double}"
        );
    }

    #[test]
    fn molar_mass_scales_inversely_with_pressure() {
        // At fixed rho, T, i: a higher osmotic pressure implies a smaller
        // (more numerous) solute, so M is proportional to 1/Pi.
        let low = molar_mass_from_pressure_atm(0.1, 8.0, 1.0, 300.0).unwrap();
        let high = molar_mass_from_pressure_atm(0.2, 8.0, 1.0, 300.0).unwrap();
        assert!(
            (high - low / 2.0).abs() / low < 1e-12,
            "high = {high}, low = {low}"
        );
    }

    #[test]
    fn molar_mass_rejects_bad_inputs() {
        assert!(molar_mass_from_pressure_atm(0.0, 10.0, 1.0, 300.0).is_err()); // Pi = 0
        assert!(molar_mass_from_pressure_atm(-1.0, 10.0, 1.0, 300.0).is_err());
        assert!(molar_mass_from_pressure_atm(0.5, 0.0, 1.0, 300.0).is_err()); // rho = 0
        assert!(molar_mass_from_pressure_atm(0.5, 10.0, 0.5, 300.0).is_err()); // i < 1
        assert!(molar_mass_from_pressure_atm(0.5, 10.0, 1.0, 0.0).is_err()); // T = 0
        assert!(molar_mass_from_pressure_atm(f64::NAN, 10.0, 1.0, 300.0).is_err());
    }
}
