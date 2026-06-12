//! Scaled 1-4 nonbonded pairs.
//!
//! AMBER / OPLS-AA force fields evaluate the nonbonded (Lennard-Jones +
//! Coulomb) interaction between the two **end atoms of each proper
//! dihedral** — the *1-4 pairs* — at a reduced strength, because the
//! bonded dihedral term already accounts for much of that interaction and
//! the original parameterisations were fitted with the reduction in place.
//! The reduction factors are [`crate::forcefield::ForceField::lj_14_scale`]
//! / [`crate::forcefield::ForceField::coulomb_14_scale`] (AMBER: 0.5 /
//! 0.8333; OPLS-AA: 0.5 / 0.5).
//!
//! These pairs are **excluded** from the ordinary cut-off LJ/Coulomb sum
//! (the [`crate::sim::Simulation`] builder adds them to the exclusion set)
//! and re-added here at scaled strength, using the **direct**
//! (un-reaction-field, un-shifted) `1/r` Coulomb and the plain 12-6 LJ —
//! the standard through-bond convention. A 1-4 pair that is *also* a 1-2 or
//! 1-3 neighbour (e.g. in a small ring) is dropped, since it is fully
//! excluded.
//!
//! Before this term existed the scale factors were dead: 1-4 pairs were
//! evaluated at *full* strength in the normal nonbonded loop, over-counting
//! 1-4 LJ ~2× and 1-4 Coulomb ~1.2-2×, which distorts torsional barriers
//! and rotamer energetics.

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::forcefield::{CombiningRule, ForceField, LjParam};
use crate::system::System;
use crate::units::COULOMB;

/// Scaled 1-4 Lennard-Jones + Coulomb pair term.
#[derive(Clone, Debug, PartialEq)]
pub struct ScaledPairs14 {
    /// The 1-4 atom-index pairs (dihedral end atoms, minus any that are
    /// also 1-2 / 1-3 neighbours).
    pairs: Vec<(usize, usize)>,
    /// Per-atom LJ parameters, indexed by atom index.
    per_atom: Vec<LjParam>,
    /// Per-atom partial charges (e), indexed by atom index.
    charges: Vec<f64>,
    /// LJ combining rule.
    combining: CombiningRule,
    /// 1-4 Lennard-Jones scale factor.
    lj_scale: f64,
    /// 1-4 Coulomb scale factor.
    coulomb_scale: f64,
}

impl ScaledPairs14 {
    /// Builds the term from a system + force field. The 1-4 pairs come
    /// from [`crate::system::Topology::one_four_pairs`].
    ///
    /// # Errors
    /// [`MdError::Invalid`] if an atom type is missing from the force field.
    pub fn from_system(system: &System, ff: &ForceField) -> Result<Self> {
        let mut per_atom = Vec::with_capacity(system.len());
        for atom in &system.topology.atoms {
            per_atom.push(ff.lj(&atom.type_name)?);
        }
        Ok(ScaledPairs14 {
            pairs: system.topology.one_four_pairs(),
            per_atom,
            charges: system.topology.atoms.iter().map(|a| a.charge).collect(),
            combining: ff.combining_rule,
            lj_scale: ff.lj_14_scale,
            coulomb_scale: ff.coulomb_14_scale,
        })
    }

    /// True when there are no 1-4 pairs to evaluate.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Number of 1-4 pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }
}

impl ForceTerm for ScaledPairs14 {
    fn name(&self) -> &str {
        "scaled-1-4-pairs"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        if out.forces.len() != n {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        for &(i, j) in &self.pairs {
            if i >= n || j >= n {
                return Err(MdError::dimension("1-4 pair index out of range"));
            }
            let d = system
                .cell
                .min_image(system.positions[i] - system.positions[j]);
            let r2 = d.norm_squared();
            if r2 < 1e-24 {
                continue;
            }
            let inv_r2 = 1.0 / r2;

            // Scaled 12-6 Lennard-Jones — direct (no cut-off shift), since
            // 1-4 pairs are intramolecular and always well within a cutoff.
            let combined = self.combining.combine(self.per_atom[i], self.per_atom[j]);
            if self.lj_scale != 0.0 && combined.epsilon > 0.0 {
                let sig2 = combined.sigma * combined.sigma;
                let sr2 = sig2 * inv_r2;
                let sr6 = sr2 * sr2 * sr2;
                let sr12 = sr6 * sr6;
                out.energy += self.lj_scale * 4.0 * combined.epsilon * (sr12 - sr6);
                let fscalar = self.lj_scale * 24.0 * combined.epsilon * (2.0 * sr12 - sr6) * inv_r2;
                let fij = fscalar * d;
                out.forces[i] += fij;
                out.forces[j] -= fij;
                out.virial += d.dot(&fij);
            }

            // Scaled direct Coulomb V = f·qᵢqⱼ/r (the reaction-field tail is
            // not applied to through-bond 1-4 pairs).
            let qq = self.charges[i] * self.charges[j];
            if self.coulomb_scale != 0.0 && qq != 0.0 {
                let inv_r = inv_r2.sqrt();
                let energy = self.coulomb_scale * COULOMB * qq * inv_r;
                out.energy += energy;
                // F = f·qᵢqⱼ/r² along d̂; f_vec = (f·qᵢqⱼ/r³)·d = energy·inv_r²·d.
                let fscalar = energy * inv_r2;
                let fij = fscalar * d;
                out.forces[i] += fij;
                out.forces[j] -= fij;
                out.virial += d.dot(&fij);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nonbonded::lj::pair_energy;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// A 0-1-2-3 chain with a single proper dihedral; the 1-4 pair is (0,3).
    fn four_atom_chain() -> (System, ForceField) {
        let mut top = Topology::new();
        for _ in 0..4 {
            top.push_atom(Atom::new("CT", 12.0, 0.2).unwrap()); // charge +0.2 each
        }
        top.add_bond(0, 1).unwrap();
        top.add_bond(1, 2).unwrap();
        top.add_bond(2, 3).unwrap();
        top.add_angle(0, 1, 2).unwrap();
        top.add_angle(1, 2, 3).unwrap();
        top.add_dihedral(0, 1, 2, 3).unwrap();
        let pos = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.15, 0.0, 0.0),
            Vector3::new(0.30, 0.10, 0.0),
            Vector3::new(0.45, 0.10, 0.0),
        ];
        let sys = System::new(top, pos)
            .unwrap()
            .with_cell(SimBox::cubic(10.0).unwrap());
        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("CT", LjParam::new(0.35, 0.3).unwrap());
        ff.lj_14_scale = 0.5;
        ff.coulomb_14_scale = 0.5;
        (sys, ff)
    }

    #[test]
    fn one_four_pair_is_the_dihedral_ends_only() {
        let (sys, _) = four_atom_chain();
        // (0,3) is the only 1-4 pair; the 1-2 (0,1),(1,2),(2,3) and 1-3
        // (0,2),(1,3) neighbours are excluded, not reported here.
        assert_eq!(sys.topology.one_four_pairs(), vec![(0, 3)]);
    }

    #[test]
    fn energy_is_the_scaled_direct_lj_plus_coulomb() {
        let (sys, ff) = four_atom_chain();
        let term = ScaledPairs14::from_system(&sys, &ff).unwrap();
        assert_eq!(term.len(), 1);
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&sys, &mut ef).unwrap();

        // Expected: 0.5 × (unscaled direct LJ + Coulomb) for the (0,3) pair.
        let r = (sys.positions[0] - sys.positions[3]).norm();
        let lj = pair_energy(0.35, 0.3, r);
        let coul = COULOMB * 0.2 * 0.2 / r;
        let expected = 0.5 * lj + 0.5 * coul;
        assert!(
            (ef.energy - expected).abs() < 1e-9,
            "got {}, expected {expected}",
            ef.energy
        );
    }

    #[test]
    fn unscaled_factors_recover_the_full_interaction() {
        // With both scale factors = 1, the 1-4 energy is the full direct
        // LJ + Coulomb — establishing that 0.5/0.5 above really is half.
        let (sys, mut ff) = four_atom_chain();
        ff.lj_14_scale = 1.0;
        ff.coulomb_14_scale = 1.0;
        let term = ScaledPairs14::from_system(&sys, &ff).unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&sys, &mut ef).unwrap();
        let r = (sys.positions[0] - sys.positions[3]).norm();
        let full = pair_energy(0.35, 0.3, r) + COULOMB * 0.2 * 0.2 / r;
        assert!(
            (ef.energy - full).abs() < 1e-9,
            "got {}, full {full}",
            ef.energy
        );
    }

    #[test]
    fn force_matches_finite_difference() {
        let (base, ff) = four_atom_chain();
        let term = ScaledPairs14::from_system(&base, &ff).unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&base, &mut ef).unwrap();
        let h = 1e-7;
        let energy_at = |dx: f64| {
            let mut s = base.clone();
            s.positions[0].x += dx;
            let mut e = EnergyForce::zeros(4);
            term.accumulate(&s, &mut e).unwrap();
            e.energy
        };
        let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!(
            (ef.forces[0].x - fd).abs() < 1e-4,
            "analytic {} vs finite-difference {fd}",
            ef.forces[0].x
        );
    }
}
