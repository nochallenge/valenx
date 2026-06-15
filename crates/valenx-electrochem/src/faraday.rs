//! Faraday's law of electrolysis: mass from charge.
//!
//! The mass of a substance deposited at (or dissolved from) an electrode by
//! passing a charge `Q` is
//!
//! ```text
//! m = (Q M) / (n F)
//! ```
//!
//! where `M` is the molar mass of the substance (grams per mole), `n` is
//! the number of electrons transferred per formula unit, and `F` is the
//! Faraday constant. The charge itself, for a constant current `I` flowing
//! for a time `t`, is `Q = I t`.
//!
//! Equivalently, the moles of substance are `Q / (n F)` and the mass is
//! that times `M`. The quantity `Q / F` is the number of moles of
//! electrons; dividing by `n` gives moles of product.
//!
//! ## Honest scope
//!
//! This assumes 100 % current efficiency (every electron drives the target
//! half-reaction). Side reactions, dendrite formation, and mass-transport
//! limits are not modelled.

use serde::{Deserialize, Serialize};

use crate::constants::FARADAY_C_PER_MOL;
use crate::error::{require_non_negative, require_positive, ElectrochemError};

/// Charge transported by a constant current, `Q = I t`, in coulombs.
///
/// `current_a` is the current in amperes and `seconds` the elapsed time in
/// seconds. Both must be non-negative (a zero current or zero time gives
/// zero charge).
///
/// # Errors
///
/// Returns [`ElectrochemError::BadParameter`] if either argument is
/// negative, or [`ElectrochemError::NonFinite`] for non-finite inputs.
pub fn charge_from_current(current_a: f64, seconds: f64) -> Result<f64, ElectrochemError> {
    let current_a = require_non_negative(
        current_a,
        "current_a",
        "current must be non-negative (amperes)",
    )?;
    let seconds = require_non_negative(seconds, "seconds", "elapsed time must be non-negative")?;
    Ok(current_a * seconds)
}

/// Moles of substance produced by a charge `Q`, equal to `Q / (n F)`.
///
/// `charge_c` is the charge in coulombs and `electrons` is `n`, the number
/// of electrons transferred per formula unit (strictly positive). The
/// result is in moles.
///
/// # Errors
///
/// Returns [`ElectrochemError::BadParameter`] if `charge_c < 0` or
/// `electrons <= 0`, or [`ElectrochemError::NonFinite`] for non-finite
/// inputs.
pub fn moles_from_charge(charge_c: f64, electrons: f64) -> Result<f64, ElectrochemError> {
    let charge_c = require_non_negative(
        charge_c,
        "charge_c",
        "charge must be non-negative (coulombs)",
    )?;
    let electrons = require_positive(
        electrons,
        "electrons",
        "electrons transferred must be strictly positive",
    )?;
    Ok(charge_c / (electrons * FARADAY_C_PER_MOL))
}

/// Mass deposited / dissolved by a charge `Q`, equal to `(Q M) / (n F)`.
///
/// `charge_c` is the charge in coulombs, `molar_mass_g_per_mol` is `M` (the
/// molar mass in grams per mole, strictly positive), and `electrons` is `n`
/// (strictly positive). The result is in grams.
///
/// The mass scales linearly with the charge and inversely with `n`: halving
/// `n` (other inputs fixed) doubles the deposited mass.
///
/// # Errors
///
/// Returns [`ElectrochemError::BadParameter`] if `charge_c < 0`,
/// `molar_mass_g_per_mol <= 0`, or `electrons <= 0`, or
/// [`ElectrochemError::NonFinite`] for non-finite inputs.
pub fn mass_from_charge(
    charge_c: f64,
    molar_mass_g_per_mol: f64,
    electrons: f64,
) -> Result<f64, ElectrochemError> {
    let molar_mass_g_per_mol = require_positive(
        molar_mass_g_per_mol,
        "molar_mass_g_per_mol",
        "molar mass must be strictly positive (g/mol)",
    )?;
    let moles = moles_from_charge(charge_c, electrons)?;
    Ok(moles * molar_mass_g_per_mol)
}

/// Mass deposited by a constant current over a time, in grams.
///
/// A convenience that chains [`charge_from_current`] into
/// [`mass_from_charge`]: it computes `Q = I t` and then
/// `m = (Q M) / (n F)`.
///
/// # Errors
///
/// Propagates any error from [`charge_from_current`] or
/// [`mass_from_charge`].
pub fn mass_from_current(
    current_a: f64,
    seconds: f64,
    molar_mass_g_per_mol: f64,
    electrons: f64,
) -> Result<f64, ElectrochemError> {
    let charge_c = charge_from_current(current_a, seconds)?;
    mass_from_charge(charge_c, molar_mass_g_per_mol, electrons)
}

/// An electrolysis deposition described by its substance and electron count.
///
/// Construct through [`Electrolysis::new`] (which validates the molar mass
/// and electron count), then ask for the [`Electrolysis::mass`] or
/// [`Electrolysis::moles`] produced by a given charge.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Electrolysis {
    /// Molar mass `M` of the deposited substance, in grams per mole.
    pub molar_mass_g_per_mol: f64,
    /// Electrons transferred per formula unit, `n` (strictly positive).
    pub electrons: f64,
}

impl Electrolysis {
    /// Build and validate an [`Electrolysis`] descriptor.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::BadParameter`] if
    /// `molar_mass_g_per_mol <= 0` or `electrons <= 0`, or
    /// [`ElectrochemError::NonFinite`] for non-finite inputs.
    pub fn new(molar_mass_g_per_mol: f64, electrons: f64) -> Result<Self, ElectrochemError> {
        let molar_mass_g_per_mol = require_positive(
            molar_mass_g_per_mol,
            "molar_mass_g_per_mol",
            "molar mass must be strictly positive (g/mol)",
        )?;
        let electrons = require_positive(
            electrons,
            "electrons",
            "electrons transferred must be strictly positive",
        )?;
        Ok(Self {
            molar_mass_g_per_mol,
            electrons,
        })
    }

    /// Moles produced by a charge `charge_c` (coulombs), `Q / (n F)`.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::BadParameter`] if `charge_c < 0`, or
    /// [`ElectrochemError::NonFinite`] if it is not finite.
    pub fn moles(&self, charge_c: f64) -> Result<f64, ElectrochemError> {
        moles_from_charge(charge_c, self.electrons)
    }

    /// Mass produced by a charge `charge_c` (coulombs), in grams.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::BadParameter`] if `charge_c < 0`, or
    /// [`ElectrochemError::NonFinite`] if it is not finite.
    pub fn mass(&self, charge_c: f64) -> Result<f64, ElectrochemError> {
        mass_from_charge(charge_c, self.molar_mass_g_per_mol, self.electrons)
    }

    /// Mass produced by a constant current `current_a` over `seconds`.
    ///
    /// Chains `Q = I t` into [`Electrolysis::mass`]; result in grams.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::BadParameter`] if `current_a < 0` or
    /// `seconds < 0`, or [`ElectrochemError::NonFinite`] for non-finite
    /// inputs.
    pub fn mass_from_current(&self, current_a: f64, seconds: f64) -> Result<f64, ElectrochemError> {
        let charge_c = charge_from_current(current_a, seconds)?;
        self.mass(charge_c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn one_faraday_deposits_one_equivalent() {
        // Passing exactly F coulombs through an n=1 substance deposits
        // exactly one mole. For silver (M = 107.868 g/mol, n = 1) that is
        // 107.868 g.
        let m = mass_from_charge(FARADAY_C_PER_MOL, 107.868, 1.0).unwrap();
        assert!(
            (m - 107.868).abs() < 1e-6,
            "1 F should deposit 1 mol Ag = 107.868 g, got {m}"
        );
        // And exactly one mole regardless of M.
        let mol = moles_from_charge(FARADAY_C_PER_MOL, 1.0).unwrap();
        assert!(
            (mol - 1.0).abs() < EPS,
            "1 F / (1 * F) should be 1 mol, got {mol}"
        );
    }

    #[test]
    fn copper_electroplating_textbook_value() {
        // 2.00 A for 1.00 hour through Cu2+ (M = 63.546 g/mol, n = 2).
        // Q = 2.00 * 3600 = 7200 C; m = 7200*63.546/(2*96485.33) ~= 2.371 g.
        let m = mass_from_current(2.00, 3600.0, 63.546, 2.0).unwrap();
        assert!(
            (m - 2.371).abs() < 1e-3,
            "Cu plating should be ~2.371 g, got {m}"
        );
    }

    #[test]
    fn charge_is_current_times_time() {
        let q = charge_from_current(2.5, 120.0).unwrap();
        assert!((q - 300.0).abs() < EPS, "Q = I t should be 300 C, got {q}");
    }

    #[test]
    fn mass_scales_linearly_with_charge() {
        // Doubling the charge doubles the mass (M, n fixed).
        let m1 = mass_from_charge(1000.0, 58.69, 2.0).unwrap();
        let m2 = mass_from_charge(2000.0, 58.69, 2.0).unwrap();
        assert!(
            (m2 - 2.0 * m1).abs() < EPS,
            "mass should double with charge: {m2} vs {}",
            2.0 * m1
        );
    }

    #[test]
    fn mass_scales_inversely_with_n() {
        // Halving n (Q, M fixed) doubles the deposited mass.
        let m_n1 = mass_from_charge(5000.0, 100.0, 1.0).unwrap();
        let m_n2 = mass_from_charge(5000.0, 100.0, 2.0).unwrap();
        assert!(
            (m_n1 - 2.0 * m_n2).abs() < EPS,
            "n=1 mass should be 2x n=2 mass: {m_n1} vs {}",
            2.0 * m_n2
        );
        let m_n4 = mass_from_charge(5000.0, 100.0, 4.0).unwrap();
        assert!(
            (m_n1 - 4.0 * m_n4).abs() < EPS,
            "n=1 mass should be 4x n=4 mass: {m_n1} vs {}",
            4.0 * m_n4
        );
    }

    #[test]
    fn moles_times_molar_mass_equals_mass() {
        // Internal consistency: mass = moles * M for the same inputs.
        let q = 4321.0;
        let m_substance = 26.98; // aluminium.
        let n = 3.0;
        let moles = moles_from_charge(q, n).unwrap();
        let mass = mass_from_charge(q, m_substance, n).unwrap();
        assert!(
            (mass - moles * m_substance).abs() < EPS,
            "mass != moles*M: {mass} vs {}",
            moles * m_substance
        );
    }

    #[test]
    fn zero_charge_gives_zero_mass() {
        let m = mass_from_charge(0.0, 63.546, 2.0).unwrap();
        assert!(m.abs() < EPS, "zero charge should give zero mass, got {m}");
        let m2 = mass_from_current(0.0, 3600.0, 63.546, 2.0).unwrap();
        assert!(
            m2.abs() < EPS,
            "zero current should give zero mass, got {m2}"
        );
    }

    #[test]
    fn struct_matches_free_functions() {
        let cell = Electrolysis::new(63.546, 2.0).unwrap();
        let q = 7200.0;
        let m_struct = cell.mass(q).unwrap();
        let m_free = mass_from_charge(q, 63.546, 2.0).unwrap();
        assert!(
            (m_struct - m_free).abs() < EPS,
            "Electrolysis::mass disagreed: {m_struct} vs {m_free}"
        );
        let mol_struct = cell.moles(q).unwrap();
        let mol_free = moles_from_charge(q, 2.0).unwrap();
        assert!((mol_struct - mol_free).abs() < EPS);
        let m_curr = cell.mass_from_current(2.0, 3600.0).unwrap();
        assert!(
            (m_curr - m_struct).abs() < EPS,
            "current path disagreed with charge path: {m_curr} vs {m_struct}"
        );
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(mass_from_charge(-1.0, 10.0, 1.0).is_err()); // Q < 0
        assert!(mass_from_charge(10.0, 0.0, 1.0).is_err()); // M = 0
        assert!(mass_from_charge(10.0, 10.0, 0.0).is_err()); // n = 0
        assert!(charge_from_current(-1.0, 10.0).is_err()); // I < 0
        assert!(charge_from_current(1.0, -10.0).is_err()); // t < 0
        assert!(Electrolysis::new(-1.0, 1.0).is_err());
        assert!(Electrolysis::new(10.0, 0.0).is_err());
        let err = moles_from_charge(10.0, f64::INFINITY).unwrap_err();
        assert_eq!(err.code(), "electrochem.non_finite");
    }
}
