//! Bonded interactions — the intramolecular force terms.
//!
//! **Roadmap features 12–15.** Four functional forms, each in its own
//! submodule, each returning a potential energy *and* the analytic
//! force on every participating atom:
//!
//! - [`bond`] — harmonic two-body bond stretching (feature 12).
//! - [`angle`] — harmonic three-body angle bending (feature 13).
//! - [`dihedral`] — proper four-body torsions, periodic and
//!   Ryckaert-Bellemans (feature 14).
//! - [`improper`] — harmonic improper dihedrals (feature 15).
//!
//! Every bonded term is evaluated with the analytic gradient, so the
//! forces are exact to machine precision (the unit tests cross-check
//! them against a central finite difference of the energy).
//!
//! ## Shared types
//!
//! [`EnergyForce`] accumulates a scalar energy and a per-atom force
//! vector. [`ForceTerm`] is the trait every interaction term —
//! bonded *and* nonbonded — implements so the [`crate::sim`] driver
//! can sum them uniformly.

pub mod angle;
pub mod bond;
pub mod dihedral;
pub mod improper;

use nalgebra::Vector3;

use crate::error::Result;
use crate::system::System;

/// An energy + per-atom force accumulator.
///
/// `forces[i]` is the total force on atom `i` (kJ/(mol·nm)); `energy`
/// is the total potential (kJ/mol). The virial tensor is accumulated
/// alongside so the pressure estimator in [`crate::analysis`] has it
/// for free.
#[derive(Clone, Debug, PartialEq)]
pub struct EnergyForce {
    /// Total potential energy (kJ/mol).
    pub energy: f64,
    /// Per-atom force (kJ/(mol·nm)).
    pub forces: Vec<Vector3<f64>>,
    /// Accumulated scalar virial `Σ rᵢⱼ · fᵢⱼ` (kJ/mol). Used by the
    /// pressure estimator.
    pub virial: f64,
}

impl EnergyForce {
    /// A zeroed accumulator sized for `n` atoms.
    pub fn zeros(n: usize) -> Self {
        EnergyForce {
            energy: 0.0,
            forces: vec![Vector3::zeros(); n],
            virial: 0.0,
        }
    }

    /// Adds `other` into `self` (energies, forces and virial sum).
    ///
    /// Panics if the force arrays differ in length — both are sized by
    /// the same system, so a mismatch is a programming error.
    pub fn add(&mut self, other: &EnergyForce) {
        assert_eq!(
            self.forces.len(),
            other.forces.len(),
            "EnergyForce::add length mismatch"
        );
        self.energy += other.energy;
        self.virial += other.virial;
        for (a, b) in self.forces.iter_mut().zip(&other.forces) {
            *a += b;
        }
    }

    /// The largest force magnitude over all atoms (kJ/(mol·nm)). Zero
    /// for an empty system. Used as a minimiser convergence gauge.
    pub fn max_force(&self) -> f64 {
        self.forces.iter().map(|f| f.norm()).fold(0.0, f64::max)
    }

    /// The root-mean-square force over all atoms.
    pub fn rms_force(&self) -> f64 {
        if self.forces.is_empty() {
            return 0.0;
        }
        let sumsq: f64 = self.forces.iter().map(|f| f.norm_squared()).sum();
        (sumsq / self.forces.len() as f64).sqrt()
    }
}

/// A force-field term that contributes energy + forces for a system.
///
/// Implemented by every bonded interaction and every nonbonded
/// interaction so the [`crate::sim::Simulation`] driver can hold a
/// `Vec<Box<dyn ForceTerm>>` and sum them.
pub trait ForceTerm {
    /// A short human-readable name (`"harmonic-bonds"`, `"lj"`, …)
    /// used in energy reports.
    fn name(&self) -> &str;

    /// Evaluates this term for `system`, accumulating into `out`.
    ///
    /// # Errors
    /// Implementation-specific — typically a dimension mismatch
    /// between the term's parameter list and the system.
    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accumulator_adds_and_reduces() {
        let mut a = EnergyForce::zeros(2);
        a.energy = 1.0;
        a.forces[0] = Vector3::new(3.0, 0.0, 0.0);
        let mut b = EnergyForce::zeros(2);
        b.energy = 2.0;
        b.forces[1] = Vector3::new(0.0, 4.0, 0.0);
        a.add(&b);
        assert_eq!(a.energy, 3.0);
        assert!((a.max_force() - 4.0).abs() < 1e-12);
        assert!(a.rms_force() > 0.0);
    }
}
