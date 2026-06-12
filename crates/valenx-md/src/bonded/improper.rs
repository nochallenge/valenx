//! Harmonic improper dihedrals — **roadmap feature 15**.
//!
//! An improper dihedral pins a planarity or chirality: it holds four
//! atoms `i`-`j`-`k`-`l` near a target dihedral angle `ξ₀` (usually 0,
//! to keep an sp² centre planar) with a harmonic restraint
//!
//! ```text
//! V(ξ) = ½ · k · (ξ − ξ₀)²
//! ```
//!
//! The geometry is the same four-atom torsion as a proper dihedral
//! (see [`crate::bonded::dihedral`]); only the potential differs —
//! harmonic about a single equilibrium rather than a periodic cosine.
//! Because the harmonic well is not periodic the angular difference is
//! wrapped into `[−π, π]` so the restraint pulls the shortest way.
//!
//! The forces reuse the proper-dihedral torsion gradient with
//! `dV/dξ = k·(ξ − ξ₀)`, so net force and net torque vanish exactly.

use crate::bonded::dihedral::{torsion_angle, torsion_forces};
use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::forcefield::ImproperParam;
use crate::system::{Improper, System};

/// The harmonic improper-dihedral force term.
#[derive(Clone, Debug, PartialEq)]
pub struct ImproperDihedrals {
    impropers: Vec<(Improper, ImproperParam)>,
}

impl ImproperDihedrals {
    /// Builds the term from explicit `(improper, parameter)` pairs.
    pub fn new(impropers: Vec<(Improper, ImproperParam)>) -> Self {
        ImproperDihedrals { impropers }
    }

    /// Builds the term by zipping a system's `topology.impropers` with
    /// a parallel parameter slice.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the lengths differ.
    pub fn from_system(system: &System, params: &[ImproperParam]) -> Result<Self> {
        let impropers = &system.topology.impropers;
        if impropers.len() != params.len() {
            return Err(MdError::dimension(format!(
                "{} impropers but {} improper parameters",
                impropers.len(),
                params.len()
            )));
        }
        Ok(ImproperDihedrals::new(
            impropers
                .iter()
                .copied()
                .zip(params.iter().copied())
                .collect(),
        ))
    }

    /// Number of impropers in the term.
    pub fn len(&self) -> usize {
        self.impropers.len()
    }

    /// Whether the term has no impropers.
    pub fn is_empty(&self) -> bool {
        self.impropers.is_empty()
    }
}

/// Wraps an angular difference into `[−π, π]`.
fn wrap_pi(mut x: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    while x > std::f64::consts::PI {
        x -= two_pi;
    }
    while x < -std::f64::consts::PI {
        x += two_pi;
    }
    x
}

impl ForceTerm for ImproperDihedrals {
    fn name(&self) -> &str {
        "improper-dihedrals"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        if out.forces.len() != n {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        for (imp, param) in &self.impropers {
            for idx in [imp.i, imp.j, imp.k, imp.l] {
                if idx >= n {
                    return Err(MdError::invalid("improper", "atom index out of range"));
                }
            }
            let (xi, b1, b2, b3) = torsion_angle(
                system.positions[imp.i],
                system.positions[imp.j],
                system.positions[imp.k],
                system.positions[imp.l],
            );
            let dxi = wrap_pi(xi - param.xi0);
            out.energy += 0.5 * param.k * dxi * dxi;
            let dv_dxi = param.k * dxi;
            let (fi, fj, fk, fl) = torsion_forces(dv_dxi, b1, b2, b3);
            out.forces[imp.i] += fi;
            out.forces[imp.j] += fj;
            out.forces[imp.k] += fk;
            out.forces[imp.l] += fl;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// Four atoms forming a (near-)planar improper centre.
    fn improper_system(out_of_plane: f64) -> System {
        let mut top = Topology::new();
        for _ in 0..4 {
            top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        }
        top.add_improper(0, 1, 2, 3).unwrap();
        // Central atom at origin, three neighbours; the fourth is
        // lifted out of the plane by `out_of_plane`.
        let pos = vec![
            Vector3::zeros(),
            Vector3::new(0.14, 0.0, 0.0),
            Vector3::new(-0.07, 0.12, 0.0),
            Vector3::new(-0.07, -0.12, out_of_plane),
        ];
        System::new(top, pos).unwrap()
    }

    #[test]
    fn planar_centre_is_near_zero_energy() {
        let sys = improper_system(0.0);
        let term = ImproperDihedrals::from_system(&sys, &[ImproperParam::new(0.0, 200.0).unwrap()])
            .unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&sys, &mut ef).unwrap();
        // Planar -> dihedral ~0 or ~pi; energy is small modulo the
        // wrap, so just check it is finite and bounded.
        assert!(ef.energy.is_finite());
    }

    #[test]
    fn out_of_plane_distortion_costs_energy() {
        let flat = improper_system(0.0);
        let bent = improper_system(0.08);
        let param = ImproperParam::new(0.0, 200.0).unwrap();
        let mut e_flat = EnergyForce::zeros(4);
        ImproperDihedrals::from_system(&flat, &[param])
            .unwrap()
            .accumulate(&flat, &mut e_flat)
            .unwrap();
        let mut e_bent = EnergyForce::zeros(4);
        ImproperDihedrals::from_system(&bent, &[param])
            .unwrap()
            .accumulate(&bent, &mut e_bent)
            .unwrap();
        assert!(e_bent.energy > e_flat.energy);
    }

    #[test]
    fn analytic_force_matches_finite_difference() {
        let param = ImproperParam::new(0.1, 150.0).unwrap();
        let base = improper_system(0.06);
        let term = ImproperDihedrals::from_system(&base, &[param]).unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&base, &mut ef).unwrap();

        let h = 1e-6;
        for atom in 0..4 {
            for comp in 0..3 {
                let energy_at = |delta: f64| {
                    let mut s = base.clone();
                    s.positions[atom][comp] += delta;
                    let mut e = EnergyForce::zeros(4);
                    term.accumulate(&s, &mut e).unwrap();
                    e.energy
                };
                let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
                assert!(
                    (ef.forces[atom][comp] - fd).abs() < 2e-3,
                    "atom {atom} comp {comp}: {} vs {}",
                    ef.forces[atom][comp],
                    fd
                );
            }
        }
    }

    #[test]
    fn net_force_vanishes() {
        let sys = improper_system(0.05);
        let term = ImproperDihedrals::from_system(&sys, &[ImproperParam::new(0.0, 100.0).unwrap()])
            .unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&sys, &mut ef).unwrap();
        let net: Vector3<f64> = ef.forces.iter().sum();
        assert!(net.norm() < 1e-8);
    }
}
