//! Cell potential from cathode and anode reduction potentials.
//!
//! A galvanic / electrolytic cell is built from two half-reactions, each
//! written as a reduction and each carrying its own reduction potential.
//! The overall cell potential is
//!
//! ```text
//! E_cell = E_cathode - E_anode
//! ```
//!
//! where both `E_cathode` and `E_anode` are *reduction* potentials (the
//! subtraction converts the anode's tabulated reduction potential into the
//! oxidation that actually happens there). A positive `E_cell` indicates a
//! thermodynamically spontaneous (galvanic) reaction as written; a negative
//! `E_cell` indicates a non-spontaneous reaction that requires an external
//! driving voltage (electrolysis).
//!
//! ## Honest scope
//!
//! This is the equilibrium (open-circuit, zero-current) cell potential from
//! ideal thermodynamics. Ohmic losses, activation / concentration
//! overpotentials, and junction potentials are not modelled, so this is not
//! a terminal voltage under load.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, ElectrochemError};

/// Whether a cell reaction (as written) is thermodynamically spontaneous.
///
/// Determined by the sign of `E_cell`: a strictly positive cell potential
/// is [`Spontaneity::Spontaneous`], a strictly negative one is
/// [`Spontaneity::NonSpontaneous`], and exactly zero is
/// [`Spontaneity::Equilibrium`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Spontaneity {
    /// `E_cell > 0`: the forward reaction releases energy (galvanic).
    Spontaneous,
    /// `E_cell < 0`: an external voltage is required (electrolytic).
    NonSpontaneous,
    /// `E_cell == 0`: the cell is at equilibrium, no net driving force.
    Equilibrium,
}

/// Cell potential `E_cell = E_cathode - E_anode`, in volts.
///
/// Both arguments are *reduction* potentials in volts. The result is the
/// open-circuit cell potential; subtract the anode's reduction potential to
/// account for oxidation occurring there.
///
/// # Errors
///
/// Returns [`ElectrochemError::NonFinite`] if either potential is not
/// finite.
pub fn cell_potential(e_cathode_v: f64, e_anode_v: f64) -> Result<f64, ElectrochemError> {
    let e_cathode_v = require_finite(e_cathode_v, "e_cathode_v")?;
    let e_anode_v = require_finite(e_anode_v, "e_anode_v")?;
    Ok(e_cathode_v - e_anode_v)
}

/// Classify a cell potential by its [`Spontaneity`].
///
/// Maps the sign of `e_cell_v` to the corresponding variant. A non-finite
/// input is rejected so that `NaN` never silently classifies as
/// equilibrium.
///
/// # Errors
///
/// Returns [`ElectrochemError::NonFinite`] if `e_cell_v` is not finite.
pub fn spontaneity(e_cell_v: f64) -> Result<Spontaneity, ElectrochemError> {
    let e_cell_v = require_finite(e_cell_v, "e_cell_v")?;
    Ok(if e_cell_v > 0.0 {
        Spontaneity::Spontaneous
    } else if e_cell_v < 0.0 {
        Spontaneity::NonSpontaneous
    } else {
        Spontaneity::Equilibrium
    })
}

/// A two-electrode cell described by its cathode and anode reduction
/// potentials.
///
/// Construct through [`Cell::new`] (which validates both potentials), then
/// query [`Cell::potential`] for `E_cell` and [`Cell::spontaneity`] for the
/// spontaneity classification.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    /// Cathode standard reduction potential, in volts (reduction site).
    pub e_cathode_v: f64,
    /// Anode standard reduction potential, in volts (oxidation site).
    pub e_anode_v: f64,
}

impl Cell {
    /// Build and validate a [`Cell`] from two reduction potentials.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::NonFinite`] if either potential is not
    /// finite.
    pub fn new(e_cathode_v: f64, e_anode_v: f64) -> Result<Self, ElectrochemError> {
        let e_cathode_v = require_finite(e_cathode_v, "e_cathode_v")?;
        let e_anode_v = require_finite(e_anode_v, "e_anode_v")?;
        Ok(Self {
            e_cathode_v,
            e_anode_v,
        })
    }

    /// The cell potential `E_cell = E_cathode - E_anode`, in volts.
    ///
    /// Always succeeds because both potentials were validated at
    /// construction time.
    pub fn potential(&self) -> f64 {
        self.e_cathode_v - self.e_anode_v
    }

    /// The [`Spontaneity`] of this cell's reaction as written.
    ///
    /// Equivalent to classifying [`Cell::potential`] by sign. Always
    /// succeeds.
    pub fn spontaneity(&self) -> Spontaneity {
        let e = self.potential();
        if e > 0.0 {
            Spontaneity::Spontaneous
        } else if e < 0.0 {
            Spontaneity::NonSpontaneous
        } else {
            Spontaneity::Equilibrium
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn daniell_cell_potential_is_about_1_1_v() {
        // Daniell cell: cathode Cu2+/Cu (+0.34 V), anode Zn2+/Zn (-0.76 V).
        // E_cell = 0.34 - (-0.76) = 1.10 V.
        let e = cell_potential(0.34, -0.76).unwrap();
        assert!(
            (e - 1.10).abs() < 1e-9,
            "Daniell cell should be ~1.10 V, got {e}"
        );
        assert_eq!(spontaneity(e).unwrap(), Spontaneity::Spontaneous);
    }

    #[test]
    fn cell_is_exactly_cathode_minus_anode() {
        // Check against a spread of arbitrary potential pairs.
        for &(cath, an) in &[
            (1.23_f64, 0.0_f64),
            (-0.13, -0.76),
            (0.80, 0.34),
            (0.0, 0.0),
        ] {
            let e = cell_potential(cath, an).unwrap();
            assert!(
                (e - (cath - an)).abs() < EPS,
                "E_cell != cathode-anode for ({cath}, {an}): got {e}"
            );
        }
    }

    #[test]
    fn negative_cell_potential_is_non_spontaneous() {
        // Reversing the Daniell electrodes gives a negative, electrolytic cell.
        let e = cell_potential(-0.76, 0.34).unwrap();
        assert!(e < 0.0, "reversed Daniell should be negative, got {e}");
        assert!((e + 1.10).abs() < 1e-9, "should be -1.10 V, got {e}");
        assert_eq!(spontaneity(e).unwrap(), Spontaneity::NonSpontaneous);
    }

    #[test]
    fn equal_potentials_give_equilibrium() {
        let e = cell_potential(0.34, 0.34).unwrap();
        assert!(e.abs() < EPS, "equal electrodes should give 0 V, got {e}");
        assert_eq!(spontaneity(e).unwrap(), Spontaneity::Equilibrium);
    }

    #[test]
    fn struct_matches_free_functions() {
        let cell = Cell::new(0.34, -0.76).unwrap();
        let e_free = cell_potential(0.34, -0.76).unwrap();
        assert!(
            (cell.potential() - e_free).abs() < EPS,
            "Cell::potential disagreed with free fn: {} vs {e_free}",
            cell.potential()
        );
        assert_eq!(cell.spontaneity(), spontaneity(e_free).unwrap());
    }

    #[test]
    fn antisymmetry_of_swapping_electrodes() {
        // Swapping cathode/anode negates the cell potential exactly.
        let forward = cell_potential(0.80, 0.15).unwrap();
        let reverse = cell_potential(0.15, 0.80).unwrap();
        assert!(
            (forward + reverse).abs() < EPS,
            "swap should negate: {forward} vs {reverse}"
        );
    }

    #[test]
    fn rejects_non_finite_potentials() {
        assert!(cell_potential(f64::NAN, 0.0).is_err());
        assert!(cell_potential(0.0, f64::INFINITY).is_err());
        assert!(Cell::new(f64::NAN, 0.0).is_err());
        let err = spontaneity(f64::NAN).unwrap_err();
        assert_eq!(err.code(), "electrochem.non_finite");
    }
}
