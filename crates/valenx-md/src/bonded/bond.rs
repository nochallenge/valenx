//! Harmonic bond stretching — **roadmap feature 12**.
//!
//! Each bond `i`-`j` carries a potential
//!
//! ```text
//! V(r) = ½ · k · (r − r₀)²
//! ```
//!
//! with `r = |rᵢ − rⱼ|`, equilibrium length `r₀` and force constant
//! `k`. The force on atom `i` is `−∂V/∂rᵢ = −k·(r − r₀)·r̂ᵢⱼ`, equal
//! and opposite on `j` — Newton's third law holds exactly.
//!
//! The bond vector is taken through the minimum-image convention so a
//! bond that straddles a periodic boundary still gives the right
//! length.

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::forcefield::BondParam;
use crate::system::{Bond, System};

/// The harmonic-bond force term.
///
/// Holds a `(Bond, BondParam)` per bond — kept independent of the
/// system's own bond list so the term can be reused / cloned.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicBonds {
    bonds: Vec<(Bond, BondParam)>,
}

impl HarmonicBonds {
    /// Builds the term from explicit `(bond, parameter)` pairs.
    pub fn new(bonds: Vec<(Bond, BondParam)>) -> Self {
        HarmonicBonds { bonds }
    }

    /// Builds the term by zipping a system's `topology.bonds` with a
    /// parallel parameter slice.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the lengths differ.
    pub fn from_system(system: &System, params: &[BondParam]) -> Result<Self> {
        let bonds = &system.topology.bonds;
        if bonds.len() != params.len() {
            return Err(MdError::dimension(format!(
                "{} bonds but {} bond parameters",
                bonds.len(),
                params.len()
            )));
        }
        Ok(HarmonicBonds::new(
            bonds.iter().copied().zip(params.iter().copied()).collect(),
        ))
    }

    /// Number of bonds in the term.
    pub fn len(&self) -> usize {
        self.bonds.len()
    }

    /// Whether the term has no bonds.
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }
}

impl ForceTerm for HarmonicBonds {
    fn name(&self) -> &str {
        "harmonic-bonds"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        if out.forces.len() != n {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        for (bond, param) in &self.bonds {
            if bond.i >= n || bond.j >= n {
                return Err(MdError::invalid("bond", "atom index out of range"));
            }
            let d = system
                .cell
                .min_image(system.positions[bond.i] - system.positions[bond.j]);
            let r = d.norm();
            if r < 1e-12 {
                // Coincident atoms: no defined direction, skip the
                // force but still flag the (large) energy.
                out.energy += 0.5 * param.k * param.r0 * param.r0;
                continue;
            }
            let dr = r - param.r0;
            out.energy += 0.5 * param.k * dr * dr;
            // f on i: -k*dr along the unit vector i->j inverse.
            let force_mag = -param.k * dr;
            let dir = d / r;
            let fi = force_mag * dir;
            out.forces[bond.i] += fi;
            out.forces[bond.j] -= fi;
            // Virial contribution r·f.
            out.virial += d.dot(&fi);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    fn two_atom_system(separation: f64) -> System {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        top.push_atom(Atom::new("B", 1.0, 0.0).unwrap());
        top.add_bond(0, 1).unwrap();
        System::new(
            top,
            vec![Vector3::zeros(), Vector3::new(separation, 0.0, 0.0)],
        )
        .unwrap()
    }

    #[test]
    fn zero_energy_at_equilibrium() {
        let sys = two_atom_system(0.15);
        let term =
            HarmonicBonds::from_system(&sys, &[BondParam::new(0.15, 1000.0).unwrap()]).unwrap();
        let mut ef = EnergyForce::zeros(2);
        term.accumulate(&sys, &mut ef).unwrap();
        assert!(ef.energy.abs() < 1e-9);
        assert!(ef.max_force() < 1e-9);
    }

    #[test]
    fn stretched_bond_pulls_atoms_together() {
        // r = 0.2, r0 = 0.15 -> stretched -> restoring force inward.
        let sys = two_atom_system(0.2);
        let term =
            HarmonicBonds::from_system(&sys, &[BondParam::new(0.15, 1000.0).unwrap()]).unwrap();
        let mut ef = EnergyForce::zeros(2);
        term.accumulate(&sys, &mut ef).unwrap();
        // V = 0.5*1000*0.05^2 = 1.25 kJ/mol.
        assert!((ef.energy - 1.25).abs() < 1e-9);
        // Atom 0 pulled toward +x (toward atom 1).
        assert!(ef.forces[0].x > 0.0);
        assert!(ef.forces[1].x < 0.0);
        // Newton's third law.
        assert!((ef.forces[0] + ef.forces[1]).norm() < 1e-9);
    }

    #[test]
    fn analytic_force_matches_finite_difference() {
        let param = BondParam::new(0.12, 2500.0).unwrap();
        let base = two_atom_system(0.17);
        let term = HarmonicBonds::from_system(&base, &[param]).unwrap();

        let mut ef = EnergyForce::zeros(2);
        term.accumulate(&base, &mut ef).unwrap();

        // Central difference of the energy in atom 0's x coordinate.
        let h = 1e-6;
        let energy_at = |dx: f64| {
            let mut s = base.clone();
            s.positions[0].x += dx;
            let mut e = EnergyForce::zeros(2);
            term.accumulate(&s, &mut e).unwrap();
            e.energy
        };
        let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!((ef.forces[0].x - fd).abs() < 1e-4, "{} vs {}", ef.forces[0].x, fd);
    }

    #[test]
    fn rejects_param_count_mismatch() {
        let sys = two_atom_system(0.15);
        assert!(HarmonicBonds::from_system(&sys, &[]).is_err());
    }
}
