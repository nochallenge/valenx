//! Holonomic constraints — SHAKE and RATTLE (**roadmap feature 23**).
//!
//! Bond stretching is the fastest motion in a molecule. Constraining
//! the stiffest bonds (typically every X-H bond) to a *fixed length*
//! removes that motion entirely and lets the integrator take a larger
//! time step. Two classic iterative constraint solvers do this:
//!
//! - **SHAKE** — works on the *positions*. After the unconstrained
//!   integrator step has moved the atoms, SHAKE iteratively adjusts
//!   them, one constraint at a time, until every constrained distance
//!   `|rᵢ − rⱼ|` is back to its target `d₀` (within tolerance). Each
//!   correction is applied along the *old* bond direction so it
//!   conserves linear and angular momentum.
//!
//! - **RATTLE** — the velocity-Verlet-compatible companion. It does
//!   the SHAKE position correction *and* a second pass that removes
//!   the component of the relative velocity along each constrained
//!   bond, so the velocities also satisfy the constraint
//!   (`(vᵢ − vⱼ)·(rᵢ − rⱼ) = 0`). RATTLE is what a velocity-Verlet
//!   run uses.
//!
//! Both are iterative Gauss-Seidel solvers; they converge
//! geometrically for the well-separated constraints of a normal
//! biomolecule. A [`MdError::NotConverged`] is returned if the
//! iteration cap is hit.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::system::System;

/// One distance constraint: hold `|r_i − r_j|` at `length`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DistanceConstraint {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Target distance (nm).
    pub length: f64,
}

impl DistanceConstraint {
    /// Builds a distance constraint.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if the two indices coincide or `length` is
    /// not finite and positive.
    pub fn new(i: usize, j: usize, length: f64) -> Result<Self> {
        if i == j {
            return Err(MdError::invalid("constraint", "indices must differ"));
        }
        if !(length.is_finite() && length > 0.0) {
            return Err(MdError::invalid(
                "constraint.length",
                "must be finite and positive",
            ));
        }
        Ok(DistanceConstraint { i, j, length })
    }
}

/// A set of distance constraints with a SHAKE / RATTLE solver.
#[derive(Clone, Debug, PartialEq)]
pub struct Constraints {
    constraints: Vec<DistanceConstraint>,
    /// Relative distance tolerance for convergence.
    tolerance: f64,
    /// Maximum solver iterations.
    max_iterations: usize,
}

impl Constraints {
    /// Builds a constraint set with the default tolerance (`1e-8`) and
    /// iteration cap (`500`).
    pub fn new(constraints: Vec<DistanceConstraint>) -> Self {
        Constraints {
            constraints,
            tolerance: 1e-8,
            max_iterations: 500,
        }
    }

    /// Builder-style override of the relative tolerance.
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance.abs().max(1e-15);
        self
    }

    /// Builder-style override of the iteration cap.
    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations.max(1);
        self
    }

    /// Number of constraints (use this as the `constraints` argument
    /// to [`System::degrees_of_freedom`](crate::system::System::degrees_of_freedom)).
    pub fn len(&self) -> usize {
        self.constraints.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.constraints.is_empty()
    }

    /// The largest relative distance violation over all constraints at
    /// the system's current positions.
    pub fn max_violation(&self, system: &System) -> f64 {
        let mut worst: f64 = 0.0;
        for c in &self.constraints {
            let d = system
                .cell
                .min_image(system.positions[c.i] - system.positions[c.j]);
            let r = d.norm();
            worst = worst.max(((r - c.length) / c.length).abs());
        }
        worst
    }

    /// **SHAKE** — corrects the *positions* so every constrained
    /// distance is satisfied.
    ///
    /// `reference` are the positions *before* the unconstrained step
    /// (SHAKE applies each correction along the old bond direction).
    /// If `velocities` is `Some`, the position shift is also folded
    /// into the velocities divided by `dt` — handy after a plain
    /// position move.
    ///
    /// # Errors
    /// [`MdError::NotConverged`] if the iteration cap is hit;
    /// [`MdError::Invalid`] on an out-of-range index.
    pub fn shake(
        &self,
        system: &mut System,
        reference: &[Vector3<f64>],
        velocities: Option<f64>,
    ) -> Result<usize> {
        let n = system.len();
        if reference.len() != n {
            return Err(MdError::dimension(
                "reference position count does not match the system",
            ));
        }
        for c in &self.constraints {
            if c.i >= n || c.j >= n {
                return Err(MdError::invalid("constraint", "atom index out of range"));
            }
        }
        for iteration in 0..self.max_iterations {
            let mut max_rel: f64 = 0.0;
            for c in &self.constraints {
                let mass_i = system.topology.atoms[c.i].mass;
                let mass_j = system.topology.atoms[c.j].mass;
                let inv_mass_i = 1.0 / mass_i;
                let inv_mass_j = 1.0 / mass_j;
                let r = system
                    .cell
                    .min_image(system.positions[c.i] - system.positions[c.j]);
                let r2 = r.norm_squared();
                let target2 = c.length * c.length;
                let diff = r2 - target2;
                max_rel = max_rel.max((diff / target2).abs());
                // Old (reference) bond direction.
                let r_old = system.cell.min_image(reference[c.i] - reference[c.j]);
                let denom = 2.0 * (inv_mass_i + inv_mass_j) * r_old.dot(&r);
                if denom.abs() < 1e-18 {
                    continue;
                }
                let g = diff / denom;
                let shift_i = -g * inv_mass_i * r_old;
                let shift_j = g * inv_mass_j * r_old;
                system.positions[c.i] += shift_i;
                system.positions[c.j] += shift_j;
                if let Some(dt) = velocities {
                    if dt > 0.0 {
                        system.velocities[c.i] += shift_i / dt;
                        system.velocities[c.j] += shift_j / dt;
                    }
                }
            }
            if max_rel < self.tolerance {
                return Ok(iteration + 1);
            }
        }
        Err(MdError::not_converged("shake", self.max_iterations))
    }

    /// **RATTLE velocity pass** — removes the component of each
    /// constrained pair's relative velocity that lies along the bond,
    /// so `(vᵢ − vⱼ)·(rᵢ − rⱼ) = 0`.
    ///
    /// Run this *after* [`shake`](Self::shake) (or after the
    /// velocity-Verlet velocity update) using the constrained
    /// positions.
    ///
    /// # Errors
    /// [`MdError::NotConverged`] if the iteration cap is hit.
    pub fn rattle_velocities(&self, system: &mut System) -> Result<usize> {
        let n = system.len();
        for c in &self.constraints {
            if c.i >= n || c.j >= n {
                return Err(MdError::invalid("constraint", "atom index out of range"));
            }
        }
        for iteration in 0..self.max_iterations {
            let mut max_proj: f64 = 0.0;
            for c in &self.constraints {
                let inv_mass_i = 1.0 / system.topology.atoms[c.i].mass;
                let inv_mass_j = 1.0 / system.topology.atoms[c.j].mass;
                let r = system
                    .cell
                    .min_image(system.positions[c.i] - system.positions[c.j]);
                let v = system.velocities[c.i] - system.velocities[c.j];
                let rv = r.dot(&v);
                let r2 = r.norm_squared().max(1e-18);
                max_proj = max_proj.max((rv / r2).abs());
                let k = rv / ((inv_mass_i + inv_mass_j) * r2);
                system.velocities[c.i] -= k * inv_mass_i * r;
                system.velocities[c.j] += k * inv_mass_j * r;
            }
            if max_proj < self.tolerance {
                return Ok(iteration + 1);
            }
        }
        Err(MdError::not_converged("rattle", self.max_iterations))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};

    /// A three-atom chain; the two bonds are constrained.
    fn chain() -> System {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 12.0, 0.0).unwrap());
        top.push_atom(Atom::new("B", 1.0, 0.0).unwrap());
        top.push_atom(Atom::new("C", 1.0, 0.0).unwrap());
        System::new(
            top,
            vec![
                Vector3::zeros(),
                Vector3::new(0.1, 0.0, 0.0),
                Vector3::new(0.2, 0.0, 0.0),
            ],
        )
        .unwrap()
    }

    #[test]
    fn distance_constraint_validates() {
        assert!(DistanceConstraint::new(0, 0, 0.1).is_err());
        assert!(DistanceConstraint::new(0, 1, 0.0).is_err());
        assert!(DistanceConstraint::new(0, 1, 0.1).is_ok());
    }

    #[test]
    fn shake_restores_constrained_lengths() {
        let reference = chain();
        let mut sys = reference.clone();
        // Perturb the atoms off their constrained lengths.
        sys.positions[1].x += 0.03;
        sys.positions[2].x -= 0.02;
        sys.positions[1].y += 0.01;

        let cons = Constraints::new(vec![
            DistanceConstraint::new(0, 1, 0.1).unwrap(),
            DistanceConstraint::new(1, 2, 0.1).unwrap(),
        ]);
        cons.shake(&mut sys, &reference.positions, None).unwrap();
        assert!(cons.max_violation(&sys) < 1e-6, "violation = {}", cons.max_violation(&sys));
    }

    #[test]
    fn shake_conserves_centre_of_mass() {
        let reference = chain();
        let mut sys = reference.clone();
        sys.positions[1] += Vector3::new(0.02, 0.02, 0.0);
        let com_before = sys.center_of_mass();
        let cons = Constraints::new(vec![
            DistanceConstraint::new(0, 1, 0.1).unwrap(),
            DistanceConstraint::new(1, 2, 0.1).unwrap(),
        ]);
        cons.shake(&mut sys, &reference.positions, None).unwrap();
        // SHAKE corrections are mass-weighted along the bond -> COM
        // is preserved.
        assert!((sys.center_of_mass() - com_before).norm() < 1e-6);
    }

    #[test]
    fn rattle_zeroes_bond_relative_velocity() {
        let mut sys = chain();
        sys.set_velocities(vec![
            Vector3::new(0.5, 0.0, 0.0),
            Vector3::new(-0.3, 0.2, 0.0),
            Vector3::new(0.1, -0.1, 0.4),
        ])
        .unwrap();
        let cons = Constraints::new(vec![
            DistanceConstraint::new(0, 1, 0.1).unwrap(),
            DistanceConstraint::new(1, 2, 0.1).unwrap(),
        ]);
        cons.rattle_velocities(&mut sys).unwrap();
        // For each constrained bond, the relative velocity along the
        // bond must vanish.
        for c in [(0usize, 1usize), (1, 2)] {
            let r = sys.positions[c.0] - sys.positions[c.1];
            let v = sys.velocities[c.0] - sys.velocities[c.1];
            assert!(r.dot(&v).abs() < 1e-7, "bond {c:?} rv = {}", r.dot(&v));
        }
    }

    #[test]
    fn shake_reports_iteration_count() {
        let reference = chain();
        let mut sys = reference.clone();
        sys.positions[1].x += 0.01;
        let cons = Constraints::new(vec![DistanceConstraint::new(0, 1, 0.1).unwrap()]);
        let iters = cons.shake(&mut sys, &reference.positions, None).unwrap();
        assert!(iters >= 1);
    }
}
