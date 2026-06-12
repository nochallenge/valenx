//! Leapfrog integrator — **roadmap feature 17**.
//!
//! The leapfrog scheme stores velocities at *half-integer* times,
//! half a step out of phase with the positions — they "leap over"
//! each other:
//!
//! ```text
//! v(t+½dt) = v(t−½dt) + dt·a(t)
//! r(t+dt)  = r(t)     + dt·v(t+½dt)
//! ```
//!
//! It is algebraically identical to velocity-Verlet (the same
//! trajectory of positions) and is the form GROMACS uses by default.
//! It is time-reversible and symplectic, so it has the same bounded-
//! energy property.
//!
//! ## On-step velocities
//!
//! Because the stored velocity is at `t−½dt`, the *kinetic energy at
//! the integer step* — what a thermostat or a temperature report
//! needs — is taken from the half-step average
//! `v(t) ≈ ½·[v(t−½dt) + v(t+½dt)]`. [`LeapFrog::step`] writes that
//! on-step velocity back into the system after each step so the rest
//! of the engine sees a consistent `r(t)`, `v(t)` pair; the
//! half-step velocity is kept internally for the next leap.

use nalgebra::Vector3;

use crate::bonded::EnergyForce;
use crate::error::{MdError, Result};
use crate::integrate::Integrator;
use crate::system::System;

/// The leapfrog integrator.
#[derive(Clone, Debug, PartialEq)]
pub struct LeapFrog {
    /// Time step (ps).
    dt: f64,
    /// Velocities at `t−½dt`. `None` before the first step, when they
    /// are seeded from the system's on-step velocities.
    v_half: Option<Vec<Vector3<f64>>>,
}

impl LeapFrog {
    /// Builds a leapfrog integrator with the given time step (ps).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `dt` is not finite and positive.
    pub fn new(dt: f64) -> Result<Self> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        Ok(LeapFrog { dt, v_half: None })
    }

    /// Discards the cached half-step velocities. Call after externally
    /// editing the system.
    pub fn reset(&mut self) {
        self.v_half = None;
    }
}

impl Integrator for LeapFrog {
    fn name(&self) -> &str {
        "leapfrog"
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
        // a(t) at the current positions.
        let ef = force_fn(system)?;
        if ef.forces.len() != n {
            return Err(MdError::dimension(
                "force evaluation returned the wrong number of forces",
            ));
        }
        let dt = self.dt;

        // v(t−½dt): cached, or seeded from the on-step velocities for
        // the very first step.
        let mut v_minus = match self.v_half.take() {
            Some(v) if v.len() == n => v,
            _ => system.velocities.clone(),
        };

        // Leap: v(t+½dt) and r(t+dt). Also form the on-step velocity.
        let mut v_onstep = vec![Vector3::zeros(); n];
        for i in 0..n {
            let mass = system.topology.atoms[i].mass;
            let a = ef.forces[i] / mass;
            let v_plus = v_minus[i] + dt * a;
            system.positions[i] += dt * v_plus;
            // On-step velocity is the half-step average.
            v_onstep[i] = 0.5 * (v_minus[i] + v_plus);
            v_minus[i] = v_plus; // becomes v(t−½dt) for the next step
        }
        system.velocities = v_onstep;
        self.v_half = Some(v_minus);
        Ok(ef)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bonded::bond::HarmonicBonds;
    use crate::bonded::ForceTerm;
    use crate::forcefield::BondParam;
    use crate::system::{Atom, Topology};

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
        assert!(LeapFrog::new(-0.001).is_err());
        assert!(LeapFrog::new(0.0).is_err());
    }

    #[test]
    fn free_particle_moves_ballistically() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        let mut sys = System::new(top, vec![Vector3::zeros()]).unwrap();
        sys.set_velocities(vec![Vector3::new(2.0, 0.0, 0.0)])
            .unwrap();
        let mut integ = LeapFrog::new(0.001).unwrap();
        let mut zero = |s: &System| Ok(EnergyForce::zeros(s.len()));
        for _ in 0..50 {
            integ.step(&mut sys, &mut zero).unwrap();
        }
        // 50 steps * 0.001 ps * 2 nm/ps = 0.1 nm.
        assert!((sys.positions[0].x - 0.1).abs() < 1e-9);
    }

    #[test]
    fn energy_is_conserved_for_an_oscillator() {
        let (mut sys, term) = oscillator();
        let mut force = |s: &System| {
            let mut ef = EnergyForce::zeros(s.len());
            term.accumulate(s, &mut ef)?;
            Ok(ef)
        };
        let mut integ = LeapFrog::new(0.0005).unwrap();
        let mut energies = Vec::new();
        for _ in 0..2000 {
            let ef = integ.step(&mut sys, &mut force).unwrap();
            energies.push(sys.kinetic_energy() + ef.energy);
        }
        let emax = energies.iter().cloned().fold(f64::MIN, f64::max);
        let emin = energies.iter().cloned().fold(f64::MAX, f64::min);
        // Bounded oscillation, no secular drift.
        assert!(emax - emin < 0.5, "energy band = {}", emax - emin);
    }
}
