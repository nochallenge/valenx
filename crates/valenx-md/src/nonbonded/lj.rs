//! Lennard-Jones interactions — **roadmap feature 7**.
//!
//! The 12-6 van der Waals potential between an unbonded pair:
//!
//! ```text
//! V_LJ(r) = 4ε · [ (σ/r)¹² − (σ/r)⁶ ]
//! ```
//!
//! Beyond a cutoff radius `r_c` the interaction is dropped. To keep
//! the energy continuous at `r_c` (a discontinuity would inject heat
//! every time a pair crosses the cutoff) the **shifted-potential**
//! form is used:
//!
//! ```text
//! V_shift(r) = V_LJ(r) − V_LJ(r_c)   for r < r_c,   else 0
//! ```
//!
//! The shift is a constant, so it changes the *energy* but not the
//! *force* — the force is `−dV_LJ/dr` for `r < r_c` and zero beyond,
//! exactly as for the plain truncated potential.
//!
//! Per-pair σ/ε come from the [`crate::forcefield::ForceField`]
//! combining rule; 1-2 / 1-3 excluded pairs are skipped via an
//! [`ExclusionSet`].

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::forcefield::{CombiningRule, ForceField, LjParam};
use crate::nonbonded::ExclusionSet;
use crate::system::System;

/// The Lennard-Jones nonbonded force term.
///
/// Holds a precomputed per-atom [`LjParam`] table and the global
/// cutoff. Pairs are taken from a caller-supplied neighbour list (or
/// all pairs when none is given).
#[derive(Clone, Debug, PartialEq)]
pub struct LennardJones {
    /// Per-atom LJ parameters, indexed by atom index.
    per_atom: Vec<LjParam>,
    /// Cutoff radius (nm).
    cutoff: f64,
    /// Whether to apply the energy shift at the cutoff.
    shifted: bool,
    /// How unlike atom types are combined into a pair interaction.
    combining: CombiningRule,
    /// Excluded 1-2 / 1-3 pairs.
    exclusions: ExclusionSet,
}

impl LennardJones {
    /// Builds the term from a force field and a system.
    ///
    /// Every atom's type is looked up in `ff` and its single-type LJ
    /// parameters cached; the standard 1-2/1-3 exclusions are derived
    /// from the topology.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if an atom type is missing from the force
    /// field or `cutoff` is non-positive.
    pub fn from_system(system: &System, ff: &ForceField, cutoff: f64) -> Result<Self> {
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(MdError::invalid("cutoff", "must be finite and positive"));
        }
        let mut per_atom = Vec::with_capacity(system.len());
        for atom in &system.topology.atoms {
            per_atom.push(ff.lj(&atom.type_name)?);
        }
        Ok(LennardJones {
            per_atom,
            cutoff,
            shifted: true,
            combining: ff.combining_rule,
            exclusions: ExclusionSet::from_topology(&system.topology),
        })
    }

    /// Builder-style switch for the cutoff energy shift (on by
    /// default).
    pub fn with_shift(mut self, shifted: bool) -> Self {
        self.shifted = shifted;
        self
    }

    /// Replaces the exclusion set.
    pub fn with_exclusions(mut self, exclusions: ExclusionSet) -> Self {
        self.exclusions = exclusions;
        self
    }

    /// The cutoff radius.
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Evaluates the LJ energy + forces over an explicit pair list.
    ///
    /// This is the workhorse: [`accumulate`](ForceTerm::accumulate)
    /// calls it with the all-pairs list; callers with a
    /// [`crate::nonbonded::neighbor::NeighborList`] pass its pairs
    /// directly to skip the `O(N²)` enumeration.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if `out` is the wrong size or a
    /// pair index is out of range.
    pub fn accumulate_pairs(
        &self,
        system: &System,
        pairs: &[(usize, usize)],
        out: &mut EnergyForce,
    ) -> Result<()> {
        let n = system.len();
        if out.forces.len() != n {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        let rc2 = self.cutoff * self.cutoff;
        for &(i, j) in pairs {
            if i >= n || j >= n {
                return Err(MdError::dimension("pair index out of range"));
            }
            if self.exclusions.contains(i, j) {
                continue;
            }
            let p = self
                .per_atom
                .get(i)
                .copied()
                .zip(self.per_atom.get(j).copied());
            let Some((pi, pj)) = p else {
                return Err(MdError::dimension("LJ table index out of range"));
            };
            let combined = self.combining.combine(pi, pj);
            let (sigma, epsilon) = (combined.sigma, combined.epsilon);
            if epsilon <= 0.0 {
                continue;
            }
            let d = system.cell.min_image(system.positions[i] - system.positions[j]);
            let r2 = d.norm_squared();
            if r2 >= rc2 || r2 < 1e-24 {
                continue;
            }
            let inv_r2 = 1.0 / r2;
            let sig2 = sigma * sigma;
            let sr2 = sig2 * inv_r2;
            let sr6 = sr2 * sr2 * sr2;
            let sr12 = sr6 * sr6;
            let mut energy = 4.0 * epsilon * (sr12 - sr6);
            if self.shifted {
                let sr2c = sig2 / rc2;
                let sr6c = sr2c * sr2c * sr2c;
                energy -= 4.0 * epsilon * (sr6c * sr6c - sr6c);
            }
            out.energy += energy;
            // Force magnitude / r: dV/dr = 4eps(-12 sr12 + 6 sr6)/r,
            // so f_vec = 24 eps (2 sr12 - sr6)/r^2 * d.
            let fscalar = 24.0 * epsilon * (2.0 * sr12 - sr6) * inv_r2;
            let fij = fscalar * d;
            out.forces[i] += fij;
            out.forces[j] -= fij;
            out.virial += d.dot(&fij);
        }
        Ok(())
    }
}

impl ForceTerm for LennardJones {
    fn name(&self) -> &str {
        "lennard-jones"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        // Default path: every pair. The Simulation driver uses a
        // neighbour list and calls `accumulate_pairs` directly.
        let n = system.len();
        let mut pairs = Vec::with_capacity(n * n.saturating_sub(1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                pairs.push((i, j));
            }
        }
        self.accumulate_pairs(system, &pairs, out)
    }
}

/// The single-pair Lennard-Jones energy (no cutoff) — handy for tests
/// and tabulation.
pub fn pair_energy(sigma: f64, epsilon: f64, r: f64) -> f64 {
    if r < 1e-24 {
        return f64::INFINITY;
    }
    let sr6 = (sigma / r).powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    fn two_atom(sep: f64) -> (System, ForceField) {
        let mut top = Topology::new();
        top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
        top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
        let sys = System::new(top, vec![Vector3::zeros(), Vector3::new(sep, 0.0, 0.0)])
            .unwrap()
            .with_cell(SimBox::cubic(10.0).unwrap());
        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("Ar", LjParam::new(0.34, 0.996).unwrap());
        (sys, ff)
    }

    #[test]
    fn minimum_is_near_two_to_the_sixth_sigma() {
        // V_LJ minimum is at r = 2^(1/6) sigma.
        let sigma = 0.34;
        let rmin = sigma * 2f64.powf(1.0 / 6.0);
        let e_min = pair_energy(sigma, 0.996, rmin);
        // At the minimum V = -epsilon.
        assert!((e_min - (-0.996)).abs() < 1e-6, "e_min = {e_min}");
    }

    #[test]
    fn repulsive_at_short_range_attractive_at_long_range() {
        let (sys_close, ff) = two_atom(0.30); // < 2^(1/6) sigma
        let lj = LennardJones::from_system(&sys_close, &ff, 1.0).unwrap();
        let mut ef = EnergyForce::zeros(2);
        lj.accumulate_pairs(&sys_close, &[(0, 1)], &mut ef).unwrap();
        // Close: atoms repel -> atom 0 pushed to -x.
        assert!(ef.forces[0].x < 0.0);

        let (sys_far, ff2) = two_atom(0.45); // > 2^(1/6) sigma
        let lj2 = LennardJones::from_system(&sys_far, &ff2, 1.0).unwrap();
        let mut ef2 = EnergyForce::zeros(2);
        lj2.accumulate_pairs(&sys_far, &[(0, 1)], &mut ef2).unwrap();
        // Far: atoms attract -> atom 0 pulled to +x.
        assert!(ef2.forces[0].x > 0.0);
    }

    #[test]
    fn shift_makes_energy_continuous_at_cutoff() {
        let (sys, ff) = two_atom(0.799); // just inside a 0.8 cutoff
        let lj = LennardJones::from_system(&sys, &ff, 0.8).unwrap();
        let mut ef = EnergyForce::zeros(2);
        lj.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        // With the shift, the energy at r just below r_c is ~0.
        assert!(ef.energy.abs() < 1e-3, "energy near cutoff = {}", ef.energy);
    }

    #[test]
    fn cutoff_drops_distant_pairs() {
        let (sys, ff) = two_atom(2.0); // far beyond a 1.0 cutoff
        let lj = LennardJones::from_system(&sys, &ff, 1.0).unwrap();
        let mut ef = EnergyForce::zeros(2);
        lj.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        assert_eq!(ef.energy, 0.0);
        assert!(ef.max_force() < 1e-12);
    }

    #[test]
    fn force_matches_finite_difference() {
        let (base, ff) = two_atom(0.38);
        let lj = LennardJones::from_system(&base, &ff, 1.0).unwrap();
        let mut ef = EnergyForce::zeros(2);
        lj.accumulate_pairs(&base, &[(0, 1)], &mut ef).unwrap();

        let h = 1e-7;
        let energy_at = |dx: f64| {
            let mut s = base.clone();
            s.positions[0].x += dx;
            let mut e = EnergyForce::zeros(2);
            lj.accumulate_pairs(&s, &[(0, 1)], &mut e).unwrap();
            e.energy
        };
        let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!((ef.forces[0].x - fd).abs() < 1e-2, "{} vs {}", ef.forces[0].x, fd);
    }

    #[test]
    fn excluded_pairs_are_skipped() {
        let (sys, ff) = two_atom(0.34);
        let mut ex = ExclusionSet::none();
        ex.insert(0, 1);
        let lj = LennardJones::from_system(&sys, &ff, 1.0)
            .unwrap()
            .with_exclusions(ex);
        let mut ef = EnergyForce::zeros(2);
        lj.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        assert_eq!(ef.energy, 0.0);
    }
}
