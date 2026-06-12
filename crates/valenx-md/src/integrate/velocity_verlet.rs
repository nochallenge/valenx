//! Velocity-Verlet integrator — **roadmap feature 16**.
//!
//! The standard MD integrator. One step of size `dt`:
//!
//! ```text
//! v(t+½dt) = v(t)    + ½·dt·a(t)
//! r(t+dt)  = r(t)    + dt·v(t+½dt)
//! a(t+dt)  = F(r(t+dt)) / m            ← one force evaluation
//! v(t+dt)  = v(t+½dt) + ½·dt·a(t+dt)
//! ```
//!
//! It is **time-reversible** and **symplectic**: the energy of a
//! conservative system does not drift, it oscillates within a bounded
//! band — the property a long MD run depends on. It is second-order
//! accurate in `dt` and needs only *one* force evaluation per step
//! (the `a(t+dt)` it computes is reused as the next step's `a(t)`).
//!
//! Acceleration is `F/m`; with `F` in kJ/(mol·nm) and `m` in u, the
//! result is in nm/ps² — consistent with the crate's unit system, no
//! conversion factor.

use nalgebra::Vector3;

use crate::bonded::EnergyForce;
use crate::error::{MdError, Result};
use crate::integrate::Integrator;
use crate::system::System;

/// The velocity-Verlet integrator.
#[derive(Clone, Debug, PartialEq)]
pub struct VelocityVerlet {
    /// Time step (ps).
    dt: f64,
    /// Cached forces from the previous step's end, reused as this
    /// step's `a(t)`. `None` before the first step.
    last_forces: Option<Vec<Vector3<f64>>>,
}

impl VelocityVerlet {
    /// Builds a velocity-Verlet integrator with the given time step
    /// (ps).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `dt` is not finite and positive.
    pub fn new(dt: f64) -> Result<Self> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        Ok(VelocityVerlet {
            dt,
            last_forces: None,
        })
    }

    /// Discards the cached forces, so the next [`step`](Integrator::step)
    /// recomputes `a(t)`. Call this after externally editing the
    /// system's positions.
    pub fn reset(&mut self) {
        self.last_forces = None;
    }
}

impl Integrator for VelocityVerlet {
    fn name(&self) -> &str {
        "velocity-verlet"
    }

    fn dt(&self) -> f64 {
        self.dt
    }

    fn step(
        &mut self,
        system: &mut System,
        force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
    ) -> Result<EnergyForce> {
        let n = system.len();
        // a(t): cached, or freshly evaluated on the first step.
        let forces_t = match self.last_forces.take() {
            Some(f) if f.len() == n => f,
            _ => force_fn(system)?.forces,
        };
        let dt = self.dt;
        let half_dt = 0.5 * dt;

        // v(t+½dt) and r(t+dt). An index loop: it touches three
        // separate borrows of `system` (atoms, velocities, positions)
        // plus the `forces_t` array, which an iterator zip cannot
        // express without splitting the struct.
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            let mass = system.topology.atoms[i].mass;
            let a = forces_t[i] / mass;
            system.velocities[i] += half_dt * a;
            let v_half = system.velocities[i];
            system.positions[i] += dt * v_half;
        }

        // a(t+dt).
        let ef_new = force_fn(system)?;
        if ef_new.forces.len() != n {
            return Err(MdError::dimension(
                "force evaluation returned the wrong number of forces",
            ));
        }

        // v(t+dt): a half-kick with the new forces.
        for ((v, atom), force) in system
            .velocities
            .iter_mut()
            .zip(&system.topology.atoms)
            .zip(&ef_new.forces)
        {
            *v += half_dt * force / atom.mass;
        }

        self.last_forces = Some(ef_new.forces.clone());
        Ok(ef_new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bonded::bond::HarmonicBonds;
    use crate::bonded::ForceTerm;
    use crate::forcefield::BondParam;
    use crate::system::{Atom, Topology};

    /// A harmonic oscillator: one mass on a spring to a fixed-ish
    /// partner. Energy should stay bounded over a long run.
    fn oscillator() -> (System, HarmonicBonds) {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 10.0, 0.0).unwrap());
        top.push_atom(Atom::new("B", 10.0, 0.0).unwrap());
        top.add_bond(0, 1).unwrap();
        let sys = System::new(top, vec![Vector3::zeros(), Vector3::new(0.12, 0.0, 0.0)]).unwrap();
        let term =
            HarmonicBonds::from_system(&sys, &[BondParam::new(0.1, 1000.0).unwrap()]).unwrap();
        (sys, term)
    }

    #[test]
    fn rejects_bad_timestep() {
        assert!(VelocityVerlet::new(0.0).is_err());
        assert!(VelocityVerlet::new(-1.0).is_err());
        assert!(VelocityVerlet::new(f64::NAN).is_err());
    }

    #[test]
    fn energy_is_conserved_for_an_oscillator() {
        let (mut sys, term) = oscillator();
        let mut force = |s: &System| {
            let mut ef = EnergyForce::zeros(s.len());
            term.accumulate(s, &mut ef)?;
            Ok(ef)
        };
        let total = |s: &System, pe: f64| s.kinetic_energy() + pe;

        let mut integ = VelocityVerlet::new(0.0005).unwrap();
        let pe0 = {
            let mut ef = EnergyForce::zeros(sys.len());
            term.accumulate(&sys, &mut ef).unwrap();
            ef.energy
        };
        let e0 = total(&sys, pe0);

        let mut max_dev: f64 = 0.0;
        for _ in 0..2000 {
            let ef = integ.step(&mut sys, &mut force).unwrap();
            let e = total(&sys, ef.energy);
            max_dev = max_dev.max((e - e0).abs());
        }
        // Symplectic integrator: bounded energy drift.
        assert!(
            max_dev < 0.05 * e0.abs().max(1.0),
            "drift = {max_dev}, e0 = {e0}"
        );
    }

    #[test]
    fn momentum_is_conserved() {
        let (mut sys, term) = oscillator();
        sys.set_velocities(vec![
            Vector3::new(0.3, 0.0, 0.0),
            Vector3::new(-0.3, 0.1, 0.0),
        ])
        .unwrap();
        let p0 = sys.linear_momentum();
        let mut force = |s: &System| {
            let mut ef = EnergyForce::zeros(s.len());
            term.accumulate(s, &mut ef)?;
            Ok(ef)
        };
        let mut integ = VelocityVerlet::new(0.001).unwrap();
        for _ in 0..500 {
            integ.step(&mut sys, &mut force).unwrap();
        }
        assert!((sys.linear_momentum() - p0).norm() < 1e-6);
    }
}
