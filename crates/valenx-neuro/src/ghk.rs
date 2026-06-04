//! Goldman–Hodgkin–Katz (GHK) voltage equation — the steady resting membrane
//! potential set by *several* permeant ions at once.
//!
//! Where the Nernst potential gives the equilibrium of a single ion, a real
//! membrane sits at the weighted compromise of every ion that can cross it,
//! each weighted by its permeability `P`:
//!
//! ```text
//!         R·T      P_K[K]_o + P_Na[Na]_o + … + P_Cl[Cl]_i + …
//! V_m  =  ───· ln ───────────────────────────────────────────
//!          F       P_K[K]_i + P_Na[Na]_i + … + P_Cl[Cl]_o + …
//! ```
//!
//! Cations contribute their *outside* concentration to the numerator and
//! *inside* to the denominator; anions (e.g. Cl⁻) flip, reflecting their
//! opposite charge. Only monovalent ions are handled — the classic GHK form.
//! Permeabilities are relative (only ratios matter) and concentrations may be
//! in any shared unit. With a single permeant ion the equation collapses back
//! to that ion's Nernst potential.

use crate::nernst::thermal_voltage_mv;

/// A monovalent permeant ion for the [`ghk_potential_mv`] equation.
///
/// `permeability` is the relative membrane permeability `P` (only ratios
/// matter); `outside` and `inside` are the extra- and intracellular
/// concentrations in any shared unit.
#[derive(Clone, Copy, Debug)]
pub struct GhkIon {
    /// Relative membrane permeability `P` (arbitrary units; only ratios matter).
    pub permeability: f64,
    /// External (extracellular) concentration.
    pub outside: f64,
    /// Internal (intracellular) concentration.
    pub inside: f64,
}

impl GhkIon {
    /// Construct a permeant ion from its permeability and its outside / inside
    /// concentrations.
    pub fn new(permeability: f64, outside: f64, inside: f64) -> GhkIon {
        GhkIon {
            permeability,
            outside,
            inside,
        }
    }
}

/// The GHK resting membrane potential in **millivolts** at absolute
/// temperature `temp_k` (K), from the permeant monovalent `cations` and
/// `anions`.
///
/// Each cation adds `P·[outside]` to the numerator and `P·[inside]` to the
/// denominator; each anion flips (`P·[inside]` up top, `P·[outside]` below).
/// At least one ion with a positive permeability and concentration must be
/// supplied — an empty mixture yields a non-finite result the caller should
/// guard.
pub fn ghk_potential_mv(temp_k: f64, cations: &[GhkIon], anions: &[GhkIon]) -> f64 {
    let numerator: f64 = cations.iter().map(|c| c.permeability * c.outside).sum::<f64>()
        + anions.iter().map(|a| a.permeability * a.inside).sum::<f64>();
    let denominator: f64 = cations.iter().map(|c| c.permeability * c.inside).sum::<f64>()
        + anions.iter().map(|a| a.permeability * a.outside).sum::<f64>();
    thermal_voltage_mv(temp_k) * (numerator / denominator).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nernst::{nernst_potential_mv, BODY_TEMPERATURE_K};

    #[test]
    fn ghk_reduces_to_nernst_for_a_single_cation() {
        // One permeant cation ⇒ the GHK potential is exactly its Nernst
        // potential, independent of the (lone) permeability value.
        let k = GhkIon::new(2.7, 5.0, 140.0);
        let ghk = ghk_potential_mv(BODY_TEMPERATURE_K, &[k], &[]);
        let nernst = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 5.0, 140.0);
        assert!((ghk - nernst).abs() < 1e-9, "GHK {ghk} vs Nernst {nernst}");
    }

    #[test]
    fn ghk_reduces_to_nernst_for_a_single_anion() {
        // One permeant anion ⇒ its Nernst potential with valence −1 (the
        // numerator/denominator flip exactly mirrors the z = −1 sign).
        let cl = GhkIon::new(1.0, 110.0, 10.0);
        let ghk = ghk_potential_mv(BODY_TEMPERATURE_K, &[], &[cl]);
        let nernst = nernst_potential_mv(BODY_TEMPERATURE_K, -1.0, 110.0, 10.0);
        assert!((ghk - nernst).abs() < 1e-9, "GHK {ghk} vs Nernst {nernst}");
    }

    #[test]
    fn resting_potential_is_physiological_and_near_e_k() {
        // A mammalian-ish mix: K dominates permeability, Na leaks a little,
        // Cl is moderately permeant. P_K : P_Na : P_Cl = 1 : 0.05 : 0.45.
        let k = GhkIon::new(1.0, 5.0, 140.0);
        let na = GhkIon::new(0.05, 145.0, 15.0);
        let cl = GhkIon::new(0.45, 110.0, 10.0);
        let vm = ghk_potential_mv(BODY_TEMPERATURE_K, &[k, na], &[cl]);
        // A resting cell sits well on the hyperpolarized side.
        assert!((-85.0..=-45.0).contains(&vm), "resting Vm {vm} mV");
        // It is bracketed by the K and Na equilibrium potentials, and sits
        // much closer to E_K because potassium dominates the permeability.
        let e_k = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 5.0, 140.0);
        let e_na = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 145.0, 15.0);
        assert!(e_k < vm && vm < e_na, "Vm {vm} must lie between E_K {e_k} and E_Na {e_na}");
        assert!((vm - e_k).abs() < (vm - e_na).abs(), "resting Vm should hug E_K");
    }

    #[test]
    fn raising_sodium_permeability_depolarises() {
        // The action-potential upstroke: when Na permeability shoots up, the
        // membrane swings from rest toward E_Na (depolarizes).
        let k = GhkIon::new(1.0, 5.0, 140.0);
        let cl = GhkIon::new(0.45, 110.0, 10.0);
        let rest = ghk_potential_mv(
            BODY_TEMPERATURE_K,
            &[k, GhkIon::new(0.05, 145.0, 15.0)],
            &[cl],
        );
        let firing = ghk_potential_mv(
            BODY_TEMPERATURE_K,
            &[k, GhkIon::new(5.0, 145.0, 15.0)],
            &[cl],
        );
        assert!(firing > rest, "more Na permeability depolarizes: {rest} -> {firing}");
        assert!(firing > 0.0, "a Na-dominated membrane overshoots positive: {firing}");
    }
}
