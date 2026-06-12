//! Coulomb electrostatics — **roadmap feature 8**.
//!
//! Two cut-off electrostatics schemes share this module:
//!
//! - **Direct (bare) Coulomb** — `V = f·qᵢqⱼ/r` for `r < r_c`, dropped
//!   beyond. `f` is the [`crate::units::COULOMB`] prefactor. Simple,
//!   but the abrupt truncation of a `1/r` tail is physically harsh.
//!
//! - **Reaction field** — the standard cure for a cut-off Coulomb
//!   sum. The medium beyond `r_c` is modelled as a uniform dielectric
//!   continuum of permittivity `ε_rf`; its polarisation response adds
//!   a term that grows as `r²`:
//!
//!   ```text
//!   V_RF(r) = f·qᵢqⱼ · [ 1/r + k_rf·r² − c_rf ]
//!   k_rf = (ε_rf − ε₁) / [ r_c³·(2·ε_rf + ε₁) ]
//!   c_rf = 1/r_c + k_rf·r_c²
//!   ```
//!
//!   `ε₁` is the permittivity *inside* the cutoff (1 for an explicit-
//!   solvent simulation). The `−c_rf` shift makes `V_RF(r_c) = 0`, so
//!   the energy is continuous; the reaction-field force
//!   `f·qᵢqⱼ·(1/r² − 2·k_rf·r)` is likewise continuous at `r_c`.
//!   `ε_rf → ∞` (a conductor) is the common "tin-foil" choice.
//!
//! Excluded 1-2 / 1-3 pairs are skipped through an [`ExclusionSet`].

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::nonbonded::ExclusionSet;
use crate::system::System;
use crate::units::COULOMB;

/// Which cut-off electrostatics model to use.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CoulombMethod {
    /// Bare truncated `1/r` Coulomb.
    Direct,
    /// Reaction-field-corrected Coulomb.
    ///
    /// `epsilon_rf` is the dielectric constant assigned to the
    /// continuum beyond the cutoff; `epsilon_inner` is the
    /// permittivity inside (usually 1). `epsilon_rf = f64::INFINITY`
    /// selects the conductor / tin-foil limit.
    ReactionField {
        /// Dielectric constant of the reaction-field continuum.
        epsilon_rf: f64,
        /// Permittivity inside the cutoff.
        epsilon_inner: f64,
    },
}

impl CoulombMethod {
    /// The conductor ("tin-foil") reaction field — the most common
    /// choice, `ε_rf → ∞`, inner permittivity 1.
    pub fn conductor_reaction_field() -> Self {
        CoulombMethod::ReactionField {
            epsilon_rf: f64::INFINITY,
            epsilon_inner: 1.0,
        }
    }
}

/// The Coulomb electrostatics force term.
#[derive(Clone, Debug, PartialEq)]
pub struct Coulomb {
    /// Per-atom partial charges (e), indexed by atom index.
    charges: Vec<f64>,
    /// Cutoff radius (nm).
    cutoff: f64,
    /// Electrostatics model.
    method: CoulombMethod,
    /// Excluded 1-2 / 1-3 pairs.
    exclusions: ExclusionSet,
}

impl Coulomb {
    /// Builds the term from a system, taking partial charges off the
    /// atoms and the 1-2/1-3 exclusions off the topology.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `cutoff` is non-positive, or the
    /// reaction-field parameters are invalid.
    pub fn from_system(system: &System, cutoff: f64, method: CoulombMethod) -> Result<Self> {
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(MdError::invalid("cutoff", "must be finite and positive"));
        }
        if let CoulombMethod::ReactionField {
            epsilon_rf,
            epsilon_inner,
        } = method
        {
            if !(epsilon_inner.is_finite() && epsilon_inner > 0.0) {
                return Err(MdError::invalid(
                    "epsilon_inner",
                    "must be finite and positive",
                ));
            }
            // epsilon_rf == INFINITY is allowed (conductor limit).
            if epsilon_rf <= 0.0 || epsilon_rf.is_nan() {
                return Err(MdError::invalid(
                    "epsilon_rf",
                    "must be positive (or +inf for the conductor limit)",
                ));
            }
        }
        Ok(Coulomb {
            charges: system.topology.atoms.iter().map(|a| a.charge).collect(),
            cutoff,
            method,
            exclusions: ExclusionSet::from_topology(&system.topology),
        })
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

    /// The reaction-field `k_rf` coefficient for this term's cutoff
    /// and method. Zero for the direct method.
    fn k_rf(&self) -> f64 {
        match self.method {
            CoulombMethod::Direct => 0.0,
            CoulombMethod::ReactionField {
                epsilon_rf,
                epsilon_inner,
            } => {
                let rc3 = self.cutoff.powi(3);
                if epsilon_rf.is_infinite() {
                    // ε_rf → ∞: k_rf -> 1 / (2 r_c³).
                    1.0 / (2.0 * rc3)
                } else {
                    (epsilon_rf - epsilon_inner) / (rc3 * (2.0 * epsilon_rf + epsilon_inner))
                }
            }
        }
    }

    /// Evaluates the Coulomb energy + forces over an explicit pair
    /// list.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] for a size / index mismatch.
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
        let rc = self.cutoff;
        let rc2 = rc * rc;
        let k_rf = self.k_rf();
        let c_rf = 1.0 / rc + k_rf * rc2;
        for &(i, j) in pairs {
            if i >= n || j >= n {
                return Err(MdError::dimension("pair index out of range"));
            }
            if self.exclusions.contains(i, j) {
                continue;
            }
            let qi = self.charges[i];
            let qj = self.charges[j];
            let qq = COULOMB * qi * qj;
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
            // Energy.
            let energy = match self.method {
                CoulombMethod::Direct => qq * inv_r,
                CoulombMethod::ReactionField { .. } => qq * (inv_r + k_rf * r2 - c_rf),
            };
            out.energy += energy;
            // Force / r along d: dV/dr = qq(-1/r^2 + 2 k_rf r),
            // f_vec = -dV/dr * d/r = qq(1/r^2 - 2 k_rf r) * d/r.
            let dv_dr = qq * (-inv_r * inv_r + 2.0 * k_rf * r);
            let fij = -(dv_dr * inv_r) * d;
            out.forces[i] += fij;
            out.forces[j] -= fij;
            out.virial += d.dot(&fij);
        }
        Ok(())
    }
}

impl ForceTerm for Coulomb {
    fn name(&self) -> &str {
        "coulomb"
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

    fn two_charges(q1: f64, q2: f64, sep: f64) -> System {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, q1).unwrap());
        top.push_atom(Atom::new("B", 1.0, q2).unwrap());
        System::new(top, vec![Vector3::zeros(), Vector3::new(sep, 0.0, 0.0)])
            .unwrap()
            .with_cell(SimBox::cubic(10.0).unwrap())
    }

    #[test]
    fn direct_coulomb_energy_is_textbook() {
        // Two +1 charges 1 nm apart: V = COULOMB ~= 139 kJ/mol.
        let sys = two_charges(1.0, 1.0, 1.0);
        let c = Coulomb::from_system(&sys, 2.0, CoulombMethod::Direct).unwrap();
        let mut ef = EnergyForce::zeros(2);
        c.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        assert!((ef.energy - COULOMB).abs() < 1e-6);
    }

    #[test]
    fn like_charges_repel_unlike_attract() {
        let like = two_charges(1.0, 1.0, 0.5);
        let c = Coulomb::from_system(&like, 2.0, CoulombMethod::Direct).unwrap();
        let mut ef = EnergyForce::zeros(2);
        c.accumulate_pairs(&like, &[(0, 1)], &mut ef).unwrap();
        // +/+ : atom 0 pushed to -x.
        assert!(ef.forces[0].x < 0.0);

        let unlike = two_charges(1.0, -1.0, 0.5);
        let c2 = Coulomb::from_system(&unlike, 2.0, CoulombMethod::Direct).unwrap();
        let mut ef2 = EnergyForce::zeros(2);
        c2.accumulate_pairs(&unlike, &[(0, 1)], &mut ef2).unwrap();
        // +/- : atom 0 pulled to +x.
        assert!(ef2.forces[0].x > 0.0);
    }

    #[test]
    fn reaction_field_energy_is_zero_at_cutoff() {
        let sys = two_charges(1.0, -1.0, 1.499);
        let c = Coulomb::from_system(&sys, 1.5, CoulombMethod::conductor_reaction_field()).unwrap();
        let mut ef = EnergyForce::zeros(2);
        c.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        assert!(
            ef.energy.abs() < 1e-2,
            "RF energy near cutoff = {}",
            ef.energy
        );
    }

    #[test]
    fn cutoff_drops_distant_pairs() {
        let sys = two_charges(1.0, 1.0, 3.0);
        let c = Coulomb::from_system(&sys, 1.0, CoulombMethod::Direct).unwrap();
        let mut ef = EnergyForce::zeros(2);
        c.accumulate_pairs(&sys, &[(0, 1)], &mut ef).unwrap();
        assert_eq!(ef.energy, 0.0);
    }

    #[test]
    fn force_matches_finite_difference_direct() {
        let base = two_charges(0.8, -0.6, 0.4);
        let c = Coulomb::from_system(&base, 2.0, CoulombMethod::Direct).unwrap();
        let mut ef = EnergyForce::zeros(2);
        c.accumulate_pairs(&base, &[(0, 1)], &mut ef).unwrap();

        let h = 1e-7;
        let energy_at = |dx: f64| {
            let mut s = base.clone();
            s.positions[0].x += dx;
            let mut e = EnergyForce::zeros(2);
            c.accumulate_pairs(&s, &[(0, 1)], &mut e).unwrap();
            e.energy
        };
        let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!(
            (ef.forces[0].x - fd).abs() < 1e-2,
            "{} vs {}",
            ef.forces[0].x,
            fd
        );
    }

    #[test]
    fn force_matches_finite_difference_reaction_field() {
        let base = two_charges(1.0, -1.0, 0.7);
        let c = Coulomb::from_system(
            &base,
            1.2,
            CoulombMethod::ReactionField {
                epsilon_rf: 78.0,
                epsilon_inner: 1.0,
            },
        )
        .unwrap();
        let mut ef = EnergyForce::zeros(2);
        c.accumulate_pairs(&base, &[(0, 1)], &mut ef).unwrap();

        let h = 1e-7;
        let energy_at = |dx: f64| {
            let mut s = base.clone();
            s.positions[0].x += dx;
            let mut e = EnergyForce::zeros(2);
            c.accumulate_pairs(&s, &[(0, 1)], &mut e).unwrap();
            e.energy
        };
        let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!(
            (ef.forces[0].x - fd).abs() < 1e-2,
            "{} vs {}",
            ef.forces[0].x,
            fd
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let sys = two_charges(1.0, 1.0, 0.5);
        assert!(Coulomb::from_system(&sys, 0.0, CoulombMethod::Direct).is_err());
        assert!(Coulomb::from_system(
            &sys,
            1.0,
            CoulombMethod::ReactionField {
                epsilon_rf: -1.0,
                epsilon_inner: 1.0
            }
        )
        .is_err());
    }
}
