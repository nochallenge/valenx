//! Wolf damped-shifted-force electrostatics — a reciprocal-space-free
//! alternative to [`super::pme`].
//!
//! Ewald summation ([`super::pme`]) is *exact* for a periodic
//! point-charge system but pays for a reciprocal-space sum over
//! `k`-vectors. The **Wolf method** (Wolf, Keefer, Asthagiri & Madden,
//! *J. Chem. Phys.* **110**, 8254 (1999)) observes that the slow
//! convergence of a bare cut-off Coulomb sum is an artefact of the cell
//! not being charge-neutral *within the cutoff sphere*: subtract the
//! net charge inside the sphere, damp the interaction with the same
//! `erfc(αr)` complementary-error screening Ewald uses, and the
//! pairwise sum alone converges to the Madelung energy — no
//! reciprocal-space term at all. The cost is `O(N·N_neighbours)`, the
//! same as a plain cut-off sum.
//!
//! This module implements the **damped shifted force (DSF)** refinement
//! of Fennell & Gezelter (*J. Chem. Phys.* **124**, 234104 (2006)),
//! which is the form everyone uses in practice because it makes *both*
//! the energy and the force continuous (and the force zero) at the
//! cutoff. The pair energy for `r < r_c` is
//!
//! ```text
//! V(r) = q_i q_j [  erfc(αr)/r
//!                 − erfc(αr_c)/r_c
//!                 + ( erfc(αr_c)/r_c²
//!                     + 2α/√π · e^{−α²r_c²}/r_c ) · (r − r_c) ]
//! ```
//!
//! The first term is the damped interaction; the second is the value
//! shift (so `V(r_c)=0`); the third is the force shift (so
//! `−dV/dr|_{r_c}=0`). The corresponding force magnitude along the pair
//! vector is
//!
//! ```text
//! −dV/dr = q_i q_j [ erfc(αr)/r² + 2α/√π · e^{−α²r²}/r
//!                    − ( erfc(αr_c)/r_c² + 2α/√π · e^{−α²r_c²}/r_c ) ]
//! ```
//!
//! which vanishes at `r = r_c` by construction.
//!
//! ## Accuracy vs Ewald
//!
//! Wolf/DSF is an **approximation**, not the exact lattice sum: it
//! trades the reciprocal-space term for a cheap pairwise shift and is
//! exact only in the limit of a large cutoff. With a sensible damping
//! `α ≈ 0.2–0.3 / nm`-scaled to the system (here in `1/nm`) and a
//! cutoff of order 1 nm it reproduces Ewald energies for ionic and
//! aqueous systems to well under a percent — the crate's tests drive
//! the NaCl lattice energy toward the Madelung value as the cutoff
//! grows. For a reference-quality energy, use [`super::pme`]; for a
//! fast, mesh-free, constant-pressure-friendly electrostatics in a
//! large run, Wolf is the pragmatic choice.
//!
//! Excluded 1-2 / 1-3 pairs are skipped through an [`ExclusionSet`],
//! exactly as in the other nonbonded terms. Forces are analytic and are
//! cross-checked against an energy finite difference in the tests.

use std::f64::consts::PI;

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::nonbonded::erfc;
use crate::nonbonded::ExclusionSet;
use crate::system::System;
use crate::units::COULOMB;

/// The Wolf damped-shifted-force electrostatics term.
#[derive(Clone, Debug, PartialEq)]
pub struct Wolf {
    /// Per-atom partial charges (e), indexed by atom index.
    charges: Vec<f64>,
    /// Cutoff radius `r_c` (nm).
    cutoff: f64,
    /// Damping parameter α (1/nm). Larger α localises the interaction
    /// faster but needs a matching cutoff to stay accurate.
    alpha: f64,
    /// Excluded 1-2 / 1-3 pairs.
    exclusions: ExclusionSet,
}

impl Wolf {
    /// Builds a Wolf term from a system, taking partial charges off the
    /// atoms and the 1-2/1-3 exclusions off the topology.
    ///
    /// * `cutoff` — real-space cutoff (nm); must be finite and positive.
    /// * `alpha` — damping parameter (1/nm); must be finite and
    ///   positive. A common rule of thumb is `α ≈ 3 / r_c`, i.e. the
    ///   damping that makes `erfc(α·r_c)` small at the cutoff.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `cutoff` or `alpha` is non-positive or
    /// non-finite.
    pub fn from_system(system: &System, cutoff: f64, alpha: f64) -> Result<Self> {
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(MdError::invalid("cutoff", "must be finite and positive"));
        }
        if !(alpha.is_finite() && alpha > 0.0) {
            return Err(MdError::invalid("alpha", "must be finite and positive"));
        }
        Ok(Wolf {
            charges: system.topology.atoms.iter().map(|a| a.charge).collect(),
            cutoff,
            alpha,
            exclusions: ExclusionSet::from_topology(&system.topology),
        })
    }

    /// Builds a Wolf term choosing the damping `α` automatically from
    /// the cutoff so that `erfc(α·r_c) ≈ accuracy`, i.e.
    /// `α = erfc⁻¹(accuracy) / r_c`. This is the same convergence
    /// criterion the Ewald real-space split uses.
    ///
    /// # Errors
    /// [`MdError::Invalid`] for a non-positive cutoff or an `accuracy`
    /// outside `(0, 1)`.
    pub fn with_accuracy(system: &System, cutoff: f64, accuracy: f64) -> Result<Self> {
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(MdError::invalid("cutoff", "must be finite and positive"));
        }
        if !(accuracy.is_finite() && accuracy > 0.0 && accuracy < 1.0) {
            return Err(MdError::invalid("accuracy", "must lie in (0, 1)"));
        }
        // erfc(x) ≈ accuracy  ⇒  x ≈ √(−ln accuracy) is the same recipe
        // used by `Pme::from_system`; α = x / r_c.
        let alpha = (-accuracy.ln()).sqrt() / cutoff;
        Self::from_system(system, cutoff, alpha)
    }

    /// Replaces the exclusion set.
    pub fn with_exclusions(mut self, exclusions: ExclusionSet) -> Self {
        self.exclusions = exclusions;
        self
    }

    /// The cutoff radius (nm).
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// The damping parameter α (1/nm).
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Evaluates the Wolf energy + forces over an explicit pair list.
    ///
    /// `pairs` is typically a neighbour list;
    /// [`accumulate`](ForceTerm::accumulate) supplies the all-pairs
    /// list. The self-energy correction (each charge's
    /// interaction with the neutralising background) is added once over
    /// all atoms regardless of the pair list.
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
        if self.charges.len() != n {
            return Err(MdError::dimension(
                "charge table size does not match the system",
            ));
        }

        let rc = self.cutoff;
        let rc2 = rc * rc;
        let alpha = self.alpha;
        let two_alpha_over_sqrtpi = 2.0 * alpha / PI.sqrt();

        // Shift constants evaluated at the cutoff (Fennell-Gezelter DSF).
        let erfc_rc = erfc(alpha * rc);
        // Value shift: erfc(αr_c)/r_c.
        let v_shift = erfc_rc / rc;
        // Force-shift magnitude: erfc(αr_c)/r_c² + 2α/√π·e^{−α²r_c²}/r_c.
        let f_shift = erfc_rc / rc2 + two_alpha_over_sqrtpi * (-(alpha * rc).powi(2)).exp() / rc;

        for &(i, j) in pairs {
            if i >= n || j >= n {
                return Err(MdError::dimension("pair index out of range"));
            }
            if self.exclusions.contains(i, j) {
                continue;
            }
            let qq = COULOMB * self.charges[i] * self.charges[j];
            if qq == 0.0 {
                continue;
            }
            let d = system
                .cell
                .min_image(system.positions[i] - system.positions[j]);
            let r2 = d.norm_squared();
            if r2 >= rc2 || r2 < 1e-24 {
                continue;
            }
            let r = r2.sqrt();
            let inv_r = 1.0 / r;
            let ar = alpha * r;
            let erfc_ar = erfc(ar);

            // DSF pair energy:
            //   V = qq [ erfc(αr)/r − erfc(αr_c)/r_c + f_shift·(r − r_c) ].
            let energy = qq * (erfc_ar * inv_r - v_shift + f_shift * (r - rc));
            out.energy += energy;

            // Force magnitude along d / r:
            //   −dV/dr = qq [ erfc(αr)/r² + 2α/√π·e^{−α²r²}/r − f_shift ].
            let neg_dv_dr = qq
                * (erfc_ar * inv_r * inv_r + two_alpha_over_sqrtpi * (-ar * ar).exp() * inv_r
                    - f_shift);
            // f_vec = (−dV/dr) · d/r.
            let fij = (neg_dv_dr * inv_r) * d;
            out.forces[i] += fij;
            out.forces[j] -= fij;
            out.virial += d.dot(&fij);
        }

        // Self-energy correction. Each charge interacts with its own
        // neutralising background; the Wolf self term (with the damping
        // limit erfc→1 at r→0 handled analytically) is
        //   E_self = −( erfc(αr_c)/(2 r_c) + α/√π ) · Σ q² · COULOMB.
        let self_q2: f64 = self.charges.iter().map(|q| q * q).sum();
        out.energy -= (erfc_rc / (2.0 * rc) + alpha / PI.sqrt()) * self_q2 * COULOMB;

        Ok(())
    }
}

impl ForceTerm for Wolf {
    fn name(&self) -> &str {
        "wolf-dsf"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// Builds a cubic NaCl rock-salt crystal of `2·half` ions per side
    /// (an even number, charge-neutral) with the given nearest-neighbour
    /// spacing, replicated so the box is periodic.
    fn nacl_crystal(half: usize, spacing: f64) -> System {
        let per_side = 2 * half;
        let edge = per_side as f64 * spacing;
        let mut top = Topology::new();
        let mut pos = Vec::new();
        for i in 0..per_side {
            for j in 0..per_side {
                for k in 0..per_side {
                    // Alternating-charge rock salt: +1 on even sublattice.
                    let charge = if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 };
                    let (name, mass) = if charge > 0.0 {
                        ("Na", 23.0)
                    } else {
                        ("Cl", 35.45)
                    };
                    top.push_atom(Atom::new(name, mass, charge).unwrap());
                    pos.push(Vector3::new(
                        i as f64 * spacing,
                        j as f64 * spacing,
                        k as f64 * spacing,
                    ));
                }
            }
        }
        System::new(top, pos)
            .unwrap()
            .with_cell(SimBox::cubic(edge).unwrap())
    }

    fn two_charges(q1: f64, q2: f64, sep: f64) -> System {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, q1).unwrap());
        top.push_atom(Atom::new("B", 1.0, q2).unwrap());
        System::new(top, vec![Vector3::zeros(), Vector3::new(sep, 0.0, 0.0)])
            .unwrap()
            .with_cell(SimBox::cubic(10.0).unwrap())
    }

    #[test]
    fn rejects_bad_parameters() {
        let sys = two_charges(1.0, -1.0, 0.3);
        assert!(Wolf::from_system(&sys, 0.0, 1.0).is_err());
        assert!(Wolf::from_system(&sys, -1.0, 1.0).is_err());
        assert!(Wolf::from_system(&sys, 1.0, 0.0).is_err());
        assert!(Wolf::from_system(&sys, 1.0, f64::NAN).is_err());
        assert!(Wolf::with_accuracy(&sys, 1.0, 0.0).is_err());
        assert!(Wolf::with_accuracy(&sys, 1.0, 1.0).is_err());
        assert!(Wolf::with_accuracy(&sys, -1.0, 1e-5).is_err());
    }

    #[test]
    fn energy_zero_at_cutoff() {
        // A single pair right at the cutoff contributes ~0 energy (the
        // value shift makes V(r_c) = 0). Add back the self term to
        // isolate the pair contribution.
        let sys = two_charges(1.0, -1.0, 1.4999);
        let wolf = Wolf::from_system(&sys, 1.5, 3.0).unwrap();
        let mut ef = EnergyForce::zeros(2);
        wolf.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        // Self energy alone:
        let mut self_only = EnergyForce::zeros(2);
        wolf.accumulate_pairs(&sys, &[], &mut self_only).unwrap();
        let pair_energy = ef.energy - self_only.energy;
        assert!(
            pair_energy.abs() < 1e-2,
            "pair energy near cutoff = {pair_energy}"
        );
    }

    #[test]
    fn opposite_charges_attract() {
        let sys = two_charges(1.0, -1.0, 0.4);
        let wolf = Wolf::from_system(&sys, 1.2, 2.5).unwrap();
        let mut ef = EnergyForce::zeros(2);
        wolf.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        // +/- : atom 0 pulled toward atom 1 (at +x), so force.x > 0.
        assert!(ef.forces[0].x > 0.0, "force = {:?}", ef.forces[0]);
    }

    #[test]
    fn force_matches_finite_difference() {
        let base = two_charges(0.8, -0.6, 0.45);
        let wolf = Wolf::from_system(&base, 1.2, 2.8).unwrap();
        let mut ef = EnergyForce::zeros(2);
        wolf.accumulate_pairs(&base, &[(0, 1)], &mut ef).unwrap();

        let h = 1e-7;
        for comp in 0..3 {
            let energy_at = |delta: f64| {
                let mut s = base.clone();
                s.positions[0][comp] += delta;
                let mut e = EnergyForce::zeros(2);
                wolf.accumulate_pairs(&s, &[(0, 1)], &mut e).unwrap();
                e.energy
            };
            let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
            assert!(
                (ef.forces[0][comp] - fd).abs() < 1e-2,
                "comp {comp}: {} vs {}",
                ef.forces[0][comp],
                fd
            );
        }
    }

    #[test]
    fn excluded_pairs_skipped() {
        let sys = two_charges(1.0, -1.0, 0.4);
        let mut ex = ExclusionSet::none();
        ex.insert(0, 1);
        let wolf = Wolf::from_system(&sys, 1.2, 2.5)
            .unwrap()
            .with_exclusions(ex);
        let mut ef = EnergyForce::zeros(2);
        wolf.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        // Only the (negative) self term remains; no pair contribution.
        let mut self_only = EnergyForce::zeros(2);
        Wolf::from_system(&sys, 1.2, 2.5)
            .unwrap()
            .with_exclusions({
                let mut e = ExclusionSet::none();
                e.insert(0, 1);
                e
            })
            .accumulate_pairs(&sys, &[], &mut self_only)
            .unwrap();
        assert!((ef.energy - self_only.energy).abs() < 1e-12);
        assert!(ef.forces[0].norm() < 1e-12);
    }

    /// **Required test 2.** The Wolf/DSF Madelung energy of a NaCl
    /// lattice converges toward the analytic Madelung value
    /// `M ≈ 1.747565` (per ion pair) as the cutoff grows. We report the
    /// per-ion-pair electrostatic energy in units of `COULOMB/spacing`
    /// (which is `−M` for the attractive ground state) and check it
    /// approaches `−M` within a documented tolerance, and that a larger
    /// cutoff is no worse than a small one.
    #[test]
    fn wolf_converges_toward_madelung() {
        const MADELUNG_NACL: f64 = 1.747_564_594_633;
        let spacing = 0.282; // ~NaCl nearest-neighbour distance (nm).

        // Reduced Madelung energy per ion pair for a cutoff = `half`
        // unit cells: E_pair / (COULOMB / spacing) → −M.
        let reduced_madelung = |half: usize, cutoff_cells: f64| -> f64 {
            let sys = nacl_crystal(half, spacing);
            let cutoff = (cutoff_cells * spacing).min(sys.cell.max_cutoff() * 0.999);
            // Damping chosen by the accuracy recipe at this cutoff.
            let wolf = Wolf::with_accuracy(&sys, cutoff, 1e-5).unwrap();
            let mut ef = EnergyForce::zeros(sys.len());
            wolf.accumulate(&sys, &mut ef).unwrap();
            let n_pairs = sys.len() as f64 / 2.0;
            let e_per_pair = ef.energy / n_pairs;
            // Divide out COULOMB/spacing to get the dimensionless value.
            e_per_pair / (COULOMB / spacing)
        };

        // Small cutoff, then a larger one in a bigger crystal.
        let small = reduced_madelung(3, 2.5);
        let large = reduced_madelung(4, 3.5);

        // Both should be near −M; the larger cutoff should be at least
        // as close. Wolf/DSF reaches ~1% at a ~1 nm cutoff for NaCl.
        assert!(
            (small + MADELUNG_NACL).abs() < 0.1,
            "small-cutoff reduced Madelung = {small} (target {})",
            -MADELUNG_NACL
        );
        assert!(
            (large + MADELUNG_NACL).abs() < 0.03,
            "large-cutoff reduced Madelung = {large} (target {})",
            -MADELUNG_NACL
        );
        assert!(
            (large + MADELUNG_NACL).abs() <= (small + MADELUNG_NACL).abs() + 1e-9,
            "larger cutoff did not improve: small={small}, large={large}"
        );
    }
}
