//! **Feature 6 — spatial-restraint model assembly (Modeller's core).**
//!
//! Modeller's defining idea is *satisfaction of spatial restraints*:
//! the template structure is not copied directly — instead, distances
//! and angles measured from the template become **restraints**, and
//! the model is the conformation that best satisfies all of them at
//! once. This naturally blends information from one (or several)
//! templates and tolerates small alignment errors.
//!
//! This module implements that idea at Cα resolution:
//!
//! 1. [`derive_restraints`] measures Cα–Cα distance restraints from
//!    the template for every equivalenced residue pair (each restraint
//!    is a harmonic well around the template distance).
//! 2. [`satisfy_restraints`] moves the model's Cα atoms — by gradient
//!    descent on the summed harmonic restraint energy — until the
//!    restraints are jointly satisfied as well as possible.
//!
//! It is a real, working restraint-satisfaction optimiser. A full
//! Modeller run additionally restrains dihedrals, uses a statistical
//! pdf rather than a plain harmonic well, and optimises all atoms; at
//! Cα level this captures the algorithm's essence.

use nalgebra::{Point3, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::{Result, StructPredictError};
use crate::homology::align::TargetTemplateAlignment;
use crate::model::ProteinModel;

/// A harmonic spatial restraint between two model residues' Cα atoms.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SpatialRestraint {
    /// Index of the first restrained residue.
    pub i: usize,
    /// Index of the second restrained residue.
    pub j: usize,
    /// Target Cα–Cα distance, ångström (measured from the template).
    pub target: f64,
    /// Harmonic force constant — a stiffer restraint when the two
    /// residues are close in the template (a contact) than when far.
    pub weight: f64,
}

impl SpatialRestraint {
    /// The harmonic energy `½·k·(d − d₀)²` of this restraint at the
    /// current Cα–Cα distance `d`.
    pub fn energy(&self, d: f64) -> f64 {
        let dx = d - self.target;
        0.5 * self.weight * dx * dx
    }
}

/// Derives Cα–Cα distance restraints from a template.
///
/// For every pair of equivalenced residues whose template Cα atoms
/// are within `contact_cutoff` ångström, a restraint is created at the
/// template distance. The weight rises as the template distance
/// falls (a tight contact is a more confident, stiffer restraint).
///
/// `template_model` is the model that carries the *template's* Cα
/// coordinates (e.g. produced by
/// [`crate::homology::transfer_backbone`] from the template onto its
/// own sequence, or the template structure itself). `alignment`
/// supplies the target↔template equivalence map.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a non-positive cutoff.
pub fn derive_restraints(
    template_model: &ProteinModel,
    alignment: &TargetTemplateAlignment,
    contact_cutoff: f64,
) -> Result<Vec<SpatialRestraint>> {
    if !(contact_cutoff.is_finite() && contact_cutoff > 0.0) {
        return Err(StructPredictError::invalid(
            "contact_cutoff",
            "must be finite and positive",
        ));
    }
    // Map template residue index → its Cα coordinate.
    let tmpl_ca: Vec<Option<Point3<f64>>> =
        template_model.residues.iter().map(|r| r.ca).collect();
    // The equivalence map gives (target_index, template_index).
    let eqs = &alignment.equivalences;
    let mut restraints = Vec::new();
    for a in 0..eqs.len() {
        for b in (a + 1)..eqs.len() {
            let (ti_a, pi_a) = eqs[a];
            let (ti_b, pi_b) = eqs[b];
            let (Some(Some(ca_a)), Some(Some(ca_b))) =
                (tmpl_ca.get(pi_a), tmpl_ca.get(pi_b))
            else {
                continue;
            };
            let d = (ca_a - ca_b).norm();
            if d > contact_cutoff {
                continue;
            }
            // Stiffer for closer contacts; never zero.
            let weight = (1.0 + (contact_cutoff - d) / contact_cutoff).max(0.1);
            restraints.push(SpatialRestraint {
                i: ti_a,
                j: ti_b,
                target: d,
                weight,
            });
        }
    }
    Ok(restraints)
}

/// The total restraint energy of a model under a restraint set.
pub fn total_restraint_energy(model: &ProteinModel, restraints: &[SpatialRestraint]) -> f64 {
    let mut e = 0.0;
    for r in restraints {
        if let (Some(Some(ci)), Some(Some(cj))) = (
            model.residues.get(r.i).map(|x| x.ca),
            model.residues.get(r.j).map(|x| x.ca),
        ) {
            e += r.energy((ci - cj).norm());
        }
    }
    e
}

/// Refines a model by satisfying a set of spatial restraints.
///
/// Moves the model's Cα atoms down the gradient of the summed
/// harmonic restraint energy. Each restraint pulls or pushes its two
/// Cα atoms toward the template distance; the optimiser iterates
/// gradient steps with an adaptive step size until the energy stops
/// decreasing or `max_iterations` is reached.
///
/// Returns the number of iterations performed. The model's Cα
/// coordinates are updated in place; non-Cα backbone atoms are *not*
/// moved here — call [`crate::refine`] for a full backbone relax.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty restraint set;
/// [`StructPredictError::NotConverged`] is *not* raised — a
/// non-converged run still returns its iteration count (the model is
/// improved, just not to tolerance).
pub fn satisfy_restraints(
    model: &mut ProteinModel,
    restraints: &[SpatialRestraint],
    max_iterations: usize,
) -> Result<usize> {
    if restraints.is_empty() {
        return Err(StructPredictError::invalid(
            "restraints",
            "no restraints to satisfy",
        ));
    }
    let n = model.residues.len();
    let mut step = 0.05;
    let mut energy = total_restraint_energy(model, restraints);
    let mut iterations = 0;
    for _ in 0..max_iterations {
        iterations += 1;
        // Gradient of the restraint energy w.r.t. each Cα.
        let mut grad = vec![Vector3::zeros(); n];
        for r in restraints {
            let (Some(ci), Some(cj)) = (
                model.residues.get(r.i).and_then(|x| x.ca),
                model.residues.get(r.j).and_then(|x| x.ca),
            ) else {
                continue;
            };
            let diff = ci - cj;
            let d = diff.norm();
            if d < 1e-9 {
                continue;
            }
            // dE/dci = k·(d - d0)·(diff/d)
            let f = r.weight * (d - r.target);
            let dir = diff / d;
            grad[r.i] += f * dir;
            grad[r.j] -= f * dir;
        }
        // Trial step downhill.
        let mut trial = model.clone();
        for (k, g) in grad.iter().enumerate() {
            if let Some(ca) = trial.residues[k].ca.as_mut() {
                *ca -= step * g;
            }
        }
        let new_energy = total_restraint_energy(&trial, restraints);
        if new_energy < energy {
            *model = trial;
            energy = new_energy;
            step *= 1.2;
        } else {
            step *= 0.5;
            if step < 1e-6 {
                break;
            }
        }
    }
    Ok(iterations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::homology::align::target_template_alignment;

    fn linear_template(n: usize, spacing: f64) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new(i as f64 * spacing, 0.0, 0.0));
        }
        m
    }

    #[test]
    fn derives_contact_restraints() {
        let tmpl = linear_template(6, 3.8);
        let aln = target_template_alignment("AAAAAA", "AAAAAA", &[]).expect("align");
        let restraints = derive_restraints(&tmpl, &aln, 8.0).expect("derive");
        assert!(!restraints.is_empty());
        // Adjacent Cα pairs are ~3.8 Å.
        let adjacent = restraints.iter().find(|r| r.j == r.i + 1).expect("pair");
        assert!((adjacent.target - 3.8).abs() < 1e-6);
    }

    #[test]
    fn satisfy_lowers_restraint_energy() {
        let tmpl = linear_template(6, 3.8);
        let aln = target_template_alignment("AAAAAA", "AAAAAA", &[]).expect("align");
        let restraints = derive_restraints(&tmpl, &aln, 10.0).expect("derive");
        // Perturbed model: Cα atoms jittered off the template line.
        let mut model = linear_template(6, 3.8);
        for (i, r) in model.residues.iter_mut().enumerate() {
            if let Some(ca) = r.ca.as_mut() {
                ca.y += if i % 2 == 0 { 1.5 } else { -1.5 };
            }
        }
        let before = total_restraint_energy(&model, &restraints);
        satisfy_restraints(&mut model, &restraints, 300).expect("satisfy");
        let after = total_restraint_energy(&model, &restraints);
        assert!(after < before, "restraint energy {before} -> {after}");
    }

    #[test]
    fn empty_restraints_rejected() {
        let mut model = linear_template(3, 3.8);
        assert!(satisfy_restraints(&mut model, &[], 10).is_err());
    }

    #[test]
    fn bad_cutoff_rejected() {
        let tmpl = linear_template(3, 3.8);
        let aln = target_template_alignment("AAA", "AAA", &[]).expect("align");
        assert!(derive_restraints(&tmpl, &aln, -1.0).is_err());
    }
}
