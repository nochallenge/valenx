//! Time integrators — advancing the system one step.
//!
//! **Roadmap features 16–18.** An integrator turns the forces at the
//! current configuration into the configuration one time step `dt`
//! later. Three classic schemes, each in its own submodule:
//!
//! - [`velocity_verlet`] — velocity-Verlet, the MD workhorse:
//!   time-reversible, symplectic, second-order, and it produces
//!   positions and velocities at the *same* instant (feature 16).
//! - [`leapfrog`] — the leapfrog scheme: positions and velocities
//!   interleaved a half-step apart. Algebraically identical to
//!   velocity-Verlet, the form GROMACS uses internally (feature 17).
//! - [`langevin`] — Langevin / Brownian stochastic dynamics: adds a
//!   friction drag and a random kick so the integrator *itself*
//!   samples the canonical ensemble at a target temperature
//!   (feature 18).
//!
//! ## The force callback
//!
//! An integrator does not know how to compute forces — that is the
//! force field's job. Each `step` takes a `force_fn` closure that maps
//! a [`System`] to an [`EnergyForce`]. The [`crate::sim`] driver
//! supplies a closure that sums every [`crate::bonded::ForceTerm`].
//!
//! ## The [`Integrator`] trait
//!
//! All three implement [`Integrator`] so the driver can hold a
//! `Box<dyn Integrator>` and swap schemes without touching the loop.

pub mod langevin;
pub mod leapfrog;
pub mod velocity_verlet;

use crate::bonded::EnergyForce;
use crate::error::Result;
use crate::system::System;

/// A time integrator: advances a [`System`] by one step `dt`.
pub trait Integrator {
    /// A short human-readable name for reports.
    fn name(&self) -> &str;

    /// The integration time step (ps).
    fn dt(&self) -> f64;

    /// Advances `system` in place by one step.
    ///
    /// `force_fn` evaluates the total force on the system; an
    /// integrator calls it once or twice per step. Returns the
    /// [`EnergyForce`] at the *new* configuration so the caller can
    /// report energies without recomputing.
    ///
    /// # Errors
    /// Propagates whatever `force_fn` returns.
    fn step(
        &mut self,
        system: &mut System,
        force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
    ) -> Result<EnergyForce>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrate::velocity_verlet::VelocityVerlet;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// A free particle under zero force keeps a constant velocity —
    /// the most basic integrator sanity check, shared across schemes.
    #[test]
    fn free_particle_moves_ballistically() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        let mut sys = System::new(top, vec![Vector3::zeros()]).unwrap();
        sys.set_velocities(vec![Vector3::new(1.0, 0.0, 0.0)]).unwrap();

        let mut integ = VelocityVerlet::new(0.001).unwrap();
        let mut zero_force = |s: &System| Ok(EnergyForce::zeros(s.len()));
        for _ in 0..100 {
            integ.step(&mut sys, &mut zero_force).unwrap();
        }
        // After 100 steps of 0.001 ps at v = 1 nm/ps: x ~= 0.1 nm.
        assert!((sys.positions[0].x - 0.1).abs() < 1e-9);
    }
}
