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
//! The same standard cell potential fixes the reaction's equilibrium
//! constant through `ΔG° = -n F E°_cell = -R T ln K`, i.e.
//! `K = exp(n F E°_cell / (R T))` ([`equilibrium_constant`]). A positive
//! `E°_cell` (`K > 1`) favours products; a negative one (`K < 1`) favours
//! reactants. At equilibrium the reaction quotient `Q` has climbed to `K`
//! and the Nernst cell potential has fallen to zero.
//!
//! ## Honest scope
//!
//! This is the equilibrium (open-circuit, zero-current) cell potential from
//! ideal thermodynamics. Ohmic losses, activation / concentration
//! overpotentials, and junction potentials are not modelled, so this is not
//! a terminal voltage under load.

use serde::{Deserialize, Serialize};

use crate::constants::{FARADAY_C_PER_MOL, GAS_CONSTANT_J_PER_MOL_K};
use crate::error::{require_finite, require_positive, ElectrochemError};

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

/// The thermodynamic equilibrium constant `K` implied by a standard cell
/// potential `E°_cell`.
///
/// Combining the two standard-state free-energy relations
/// `ΔG° = -n F E°_cell` (electrochemistry) and `ΔG° = -R T ln K`
/// (thermodynamics) eliminates `ΔG°` and gives
///
/// ```text
/// ln K = n F E°_cell / (R T)   =>   K = exp(n F E°_cell / (R T))
/// ```
///
/// where `n` is the number of electrons transferred in the balanced overall
/// reaction, `F` the Faraday constant, `R` the molar gas constant, and `T`
/// the absolute temperature. A spontaneous cell (`E°_cell > 0`) gives
/// `K > 1` (products favoured); a non-spontaneous one (`E°_cell < 0`) gives
/// `K < 1`; and `E°_cell = 0` gives `K = 1`. At `Q = K` the Nernst cell
/// potential is exactly zero, the defining condition of equilibrium.
///
/// Like the rest of the crate this uses ideal activities, so `K` is the
/// thermodynamic (activity-based) equilibrium constant.
///
/// # Errors
///
/// Returns [`ElectrochemError::NonFinite`] if `e_cell_standard_v` is not
/// finite, or [`ElectrochemError::BadParameter`] if `electrons <= 0` or
/// `temperature_k <= 0`.
pub fn equilibrium_constant(
    e_cell_standard_v: f64,
    electrons: f64,
    temperature_k: f64,
) -> Result<f64, ElectrochemError> {
    let e = require_finite(e_cell_standard_v, "e_cell_standard_v")?;
    let n = require_positive(
        electrons,
        "electrons",
        "electrons transferred must be strictly positive",
    )?;
    let t = require_positive(
        temperature_k,
        "temperature_k",
        "absolute temperature must be strictly positive (kelvin)",
    )?;
    Ok((n * FARADAY_C_PER_MOL * e / (GAS_CONSTANT_J_PER_MOL_K * t)).exp())
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

    /// The equilibrium constant `K = exp(n F E°_cell / (R T))` for this
    /// cell's reaction, treating the stored potentials as standard ones.
    ///
    /// Equivalent to [`equilibrium_constant`] applied to
    /// [`Cell::potential`]. `electrons` is the number transferred in the
    /// balanced overall reaction and `temperature_k` the absolute
    /// temperature.
    ///
    /// # Errors
    ///
    /// Returns [`ElectrochemError::BadParameter`] if `electrons <= 0` or
    /// `temperature_k <= 0`. (The potentials were validated at construction,
    /// so their difference is always finite.)
    pub fn equilibrium_constant(
        &self,
        electrons: f64,
        temperature_k: f64,
    ) -> Result<f64, ElectrochemError> {
        equilibrium_constant(self.potential(), electrons, temperature_k)
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

    // --- Equilibrium constant -------------------------------------------

    #[test]
    fn daniell_cell_equilibrium_constant_is_about_1e37() {
        // Daniell cell: E°_cell = 1.10 V, n = 2 electrons, 298.15 K.
        // ln K = nFE°/RT = 85.63 => log10 K = 37.19 => K ~ 1.5e37, matching
        // the textbook "K ~ 10^37" for Zn + Cu2+ -> Zn2+ + Cu.
        let k = equilibrium_constant(1.10, 2.0, 298.15).unwrap();
        assert!(k > 0.0, "K must be positive, got {k}");
        let log10_k = k.log10();
        assert!(
            (log10_k - 37.19).abs() < 0.05,
            "Daniell K should be ~10^37.2, got 10^{log10_k}"
        );
    }

    #[test]
    fn at_q_equals_k_the_nernst_cell_potential_is_zero() {
        // The defining property of K: when the reaction quotient climbs to
        // K, the Nernst cell potential collapses to zero. Feed K straight
        // back into the (independently implemented) Nernst equation as Q and
        // confirm E -> 0, across spontaneous and non-spontaneous cells.
        for &(e0, n, t) in &[
            (1.10_f64, 2.0_f64, 298.15_f64),
            (0.46, 1.0, 310.0),
            (-0.25, 3.0, 273.15),
        ] {
            let k = equilibrium_constant(e0, n, t).unwrap();
            let e_at_equilibrium = crate::nernst::nernst_potential(e0, n, t, k).unwrap();
            assert!(
                e_at_equilibrium.abs() < 1e-9,
                "E at Q=K should be 0 (e0={e0}, n={n}, T={t}), got {e_at_equilibrium}"
            );
        }
    }

    #[test]
    fn zero_cell_potential_gives_unit_equilibrium_constant() {
        // E°_cell = 0 => ln K = 0 => K = 1 exactly, for any valid n, T.
        for &(n, t) in &[(1.0_f64, 298.15_f64), (4.0, 350.0)] {
            let k = equilibrium_constant(0.0, n, t).unwrap();
            assert!((k - 1.0).abs() < 1e-12, "K should be 1 at E°=0, got {k}");
        }
    }

    #[test]
    fn positive_potential_favours_products_negative_favours_reactants() {
        let t = 298.15;
        let k_pos = equilibrium_constant(0.5, 2.0, t).unwrap();
        let k_neg = equilibrium_constant(-0.5, 2.0, t).unwrap();
        assert!(k_pos > 1.0, "spontaneous cell should give K>1, got {k_pos}");
        assert!(
            k_neg < 1.0,
            "non-spontaneous cell should give K<1, got {k_neg}"
        );
        // Reversing the sign of E° inverts K (ln K changes sign), so the
        // product of the two is exactly one.
        assert!(
            (k_pos * k_neg - 1.0).abs() < 1e-9,
            "K(+E°)*K(-E°) should be 1: {k_pos} * {k_neg}"
        );
    }

    #[test]
    fn doubling_electrons_squares_the_constant() {
        // ln K is linear in n, so K(2n) = K(n)^2 at fixed E°, T.
        let (e0, t) = (0.1, 298.15);
        let k1 = equilibrium_constant(e0, 1.0, t).unwrap();
        let k2 = equilibrium_constant(e0, 2.0, t).unwrap();
        assert!(
            (k2 / (k1 * k1) - 1.0).abs() < 1e-9,
            "K(n=2) should equal K(n=1)^2: {k2} vs {}",
            k1 * k1
        );
    }

    #[test]
    fn cell_method_matches_free_function() {
        // Cell::equilibrium_constant routes through the free function with
        // E°_cell = cathode - anode.
        let cell = Cell::new(0.34, -0.76).unwrap(); // Daniell, E° = 1.10 V
        let via_method = cell.equilibrium_constant(2.0, 298.15).unwrap();
        let via_free = equilibrium_constant(1.10, 2.0, 298.15).unwrap();
        assert!(
            (via_method / via_free - 1.0).abs() < 1e-12,
            "method vs free fn: {via_method} vs {via_free}"
        );
    }

    #[test]
    fn equilibrium_constant_rejects_bad_inputs() {
        assert!(equilibrium_constant(f64::NAN, 2.0, 298.15).is_err()); // bad E°
        assert!(equilibrium_constant(f64::INFINITY, 2.0, 298.15).is_err());
        assert!(equilibrium_constant(1.0, 0.0, 298.15).is_err()); // n = 0
        assert!(equilibrium_constant(1.0, -2.0, 298.15).is_err()); // n < 0
        assert!(equilibrium_constant(1.0, 2.0, 0.0).is_err()); // T = 0
        assert!(equilibrium_constant(1.0, 2.0, -5.0).is_err()); // T < 0
        let err = equilibrium_constant(f64::NAN, 2.0, 298.15).unwrap_err();
        assert_eq!(err.code(), "electrochem.non_finite");
    }
}
