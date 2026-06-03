//! **Feature 14 — energy minimisation / relaxation.**
//!
//! A built model — even after loop closure and repacking — sits in a
//! strained, slightly-clashing geometry. **Energy minimisation**
//! relaxes it: it walks the structure downhill on a potential-energy
//! surface until the forces are small, removing the worst strain
//! without changing the fold.
//!
//! This module drives the minimisation with the **[`valenx_md`]
//! force-field minimisers** — the same callback-based steepest-descent
//! / conjugate-gradient / L-BFGS optimisers `valenx-md` uses to
//! prepare an MD system. The model's Cα atoms become a `valenx_md`
//! `System`; the energy + force callback evaluates a **structure
//! restraint potential**:
//!
//! - a **chain term** — a harmonic well on every consecutive Cα–Cα
//!   distance about the ideal ~3.8 Å virtual-bond length (keeps the
//!   chain connected);
//! - a **non-bonded term** — a soft repulsion between non-adjacent
//!   Cα atoms that are too close (removes clashes).
//!
//! The minimiser itself is `valenx-md`'s, exercised through its
//! public `EnergyForce`-callback interface — a real reuse of the
//! force-field engine, not a re-implementation.

use nalgebra::Vector3;
use valenx_md::bonded::EnergyForce;
use valenx_md::minimize::{steepest_descent, MinimizeOptions};
use valenx_md::system::{Atom, Topology};
use valenx_md::System;

use crate::error::{Result, StructPredictError};
use crate::model::{ideal, ProteinModel};

/// The outcome of a relaxation run.
#[derive(Clone, Debug, PartialEq)]
pub struct RelaxResult {
    /// Structure-restraint energy before relaxation.
    pub initial_energy: f64,
    /// Structure-restraint energy after relaxation.
    pub final_energy: f64,
    /// Maximum force on any Cα at the end.
    pub final_max_force: f64,
    /// Minimiser iterations performed.
    pub iterations: usize,
    /// Whether the force tolerance was reached.
    pub converged: bool,
}

/// Ideal virtual Cα–Cα bond length (ångström) — the chain restraint
/// target.
const CA_BOND: f64 = ideal::CA_CA;
/// Harmonic force constant for the Cα–Cα chain restraint.
const CHAIN_K: f64 = 50.0;
/// Soft-clash radius for non-adjacent Cα atoms (ångström).
const CLASH_R: f64 = 4.5;
/// Soft-clash force constant.
const CLASH_K: f64 = 20.0;

/// Evaluates the structure-restraint energy and forces on a Cα set.
///
/// This is the energy + force callback handed to the `valenx_md`
/// minimiser. `positions` are the current Cα coordinates.
fn structure_energy_force(positions: &[Vector3<f64>]) -> EnergyForce {
    let n = positions.len();
    let mut ef = EnergyForce::zeros(n);
    // Chain term: harmonic on consecutive Cα–Cα distances.
    for i in 0..n.saturating_sub(1) {
        let d = positions[i + 1] - positions[i];
        let len = d.norm();
        if len < 1e-9 {
            continue;
        }
        let dx = len - CA_BOND;
        ef.energy += 0.5 * CHAIN_K * dx * dx;
        // Force on atom i: +k·dx·(d/len); equal and opposite on i+1.
        let f = CHAIN_K * dx * (d / len);
        ef.forces[i] += f;
        ef.forces[i + 1] -= f;
    }
    // Non-bonded soft clash for non-adjacent pairs.
    for i in 0..n {
        for j in (i + 2)..n {
            let d = positions[j] - positions[i];
            let len = d.norm();
            if !(1e-9..CLASH_R).contains(&len) {
                continue;
            }
            let pen = CLASH_R - len;
            ef.energy += 0.5 * CLASH_K * pen * pen;
            // Repulsive: push i and j apart.
            let f = CLASH_K * pen * (d / len);
            ef.forces[i] -= f;
            ef.forces[j] += f;
        }
    }
    ef
}

/// Builds a bond-free `valenx_md` `System` whose atoms are the
/// model's Cα atoms.
fn ca_system(model: &ProteinModel) -> Result<(System, Vec<usize>)> {
    let mut topology = Topology::new();
    let mut positions = Vec::new();
    let mut residue_of_atom = Vec::new();
    for (i, r) in model.residues.iter().enumerate() {
        if let Some(ca) = r.ca {
            // A 12 g/mol pseudo-atom — mass is irrelevant to a
            // minimiser, which only needs energy + force.
            let atom = Atom::new("CA", 12.0, 0.0)
                .map_err(|e| StructPredictError::invalid("md_atom", e.to_string()))?;
            topology.push_atom(atom);
            positions.push(Vector3::new(ca.x, ca.y, ca.z));
            residue_of_atom.push(i);
        }
    }
    if positions.len() < 2 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 2 Cα atoms to relax",
        ));
    }
    let system = System::new(topology, positions)
        .map_err(|e| StructPredictError::invalid("md_system", e.to_string()))?;
    Ok((system, residue_of_atom))
}

/// Relaxes a model by gradient energy minimisation.
///
/// Builds a `valenx_md` `System` from the Cα trace, runs
/// `valenx-md`'s steepest-descent minimiser against the
/// structure-restraint potential, and writes the relaxed Cα
/// coordinates back into the model. Non-Cα backbone atoms are *not*
/// moved — call this after the backbone is built but treat it as a
/// Cα-level relax.
///
/// `max_iterations` caps the minimiser; `force_tolerance` is the
/// stopping force.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a model with fewer than 2 Cα
/// atoms or bad options.
pub fn relax_model(
    model: &mut ProteinModel,
    max_iterations: usize,
    force_tolerance: f64,
) -> Result<RelaxResult> {
    if max_iterations == 0 {
        return Err(StructPredictError::invalid(
            "max_iterations",
            "must be at least 1",
        ));
    }
    if !(force_tolerance.is_finite() && force_tolerance > 0.0) {
        return Err(StructPredictError::invalid(
            "force_tolerance",
            "must be finite and positive",
        ));
    }
    let (mut system, residue_of_atom) = ca_system(model)?;

    let initial_energy = structure_energy_force(&system.positions).energy;

    let options = MinimizeOptions {
        force_tolerance,
        max_iterations,
        initial_step: 0.01,
    };
    // The callback `valenx-md` drives — note `valenx-md` passes a
    // whole `System`; we score its positions.
    let mut force_fn = |sys: &System| -> valenx_md::Result<EnergyForce> {
        Ok(structure_energy_force(&sys.positions))
    };
    let result = steepest_descent(&mut system, options, &mut force_fn)
        .map_err(|e| StructPredictError::invalid("minimize", e.to_string()))?;

    // Write the relaxed Cα coordinates back into the model.
    for (atom_idx, &res_idx) in residue_of_atom.iter().enumerate() {
        let p = system.positions[atom_idx];
        if let Some(ca) = model.residues[res_idx].ca.as_mut() {
            *ca = nalgebra::Point3::new(p.x, p.y, p.z);
        }
    }

    Ok(RelaxResult {
        initial_energy,
        final_energy: result.final_energy,
        final_max_force: result.final_max_force,
        iterations: result.iterations,
        converged: result.converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    /// A chain with one bad Cα–Cα contact (residue 2 sits on top of
    /// residue 1).
    fn strained_model() -> ProteinModel {
        let mut m = ProteinModel::from_sequence("AAAAAA").expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new(i as f64 * 3.8, 0.0, 0.0));
        }
        // Strain: pull residue 3 far off the chain.
        m.residues[3].ca = Some(Point3::new(3.0 * 3.8, 9.0, 0.0));
        m
    }

    #[test]
    fn relaxation_lowers_the_energy() {
        let mut m = strained_model();
        let res = relax_model(&mut m, 800, 1.0).expect("relax");
        assert!(
            res.final_energy <= res.initial_energy,
            "energy {} -> {}",
            res.initial_energy,
            res.final_energy
        );
    }

    #[test]
    fn relaxation_pulls_the_chain_toward_ideal_spacing() {
        let mut m = strained_model();
        relax_model(&mut m, 1500, 0.5).expect("relax");
        // After relaxing, the residue-2→3 virtual bond is closer to
        // the ideal ~3.8 Å than the 9-Å-off starting geometry.
        let d = (m.residues[3].ca.unwrap() - m.residues[2].ca.unwrap()).norm();
        assert!(d < 7.0, "Cα-Cα spacing relaxed toward ideal: {d}");
    }

    #[test]
    fn too_short_model_rejected() {
        let mut m = ProteinModel::from_sequence("A").expect("model");
        m.residues[0].ca = Some(Point3::origin());
        assert!(relax_model(&mut m, 100, 1.0).is_err());
    }

    #[test]
    fn bad_options_rejected() {
        let mut m = strained_model();
        assert!(relax_model(&mut m, 0, 1.0).is_err());
        assert!(relax_model(&mut m, 10, -1.0).is_err());
    }
}
