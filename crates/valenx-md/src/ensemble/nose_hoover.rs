//! Nosé-Hoover (chain) thermostat — **roadmap feature 21**.
//!
//! The Nosé-Hoover thermostat couples the system to a heat bath
//! through an *extended-system* variable: a fictitious "thermostat
//! coordinate" `η` with its own momentum `p_η` and mass `Q`. The
//! thermostat momentum obeys
//!
//! ```text
//! dp_η/dt = 2·KE − N_f·k_B·T₀
//! ```
//!
//! — it accelerates when the system is too hot and decelerates when
//! too cold — and feeds back a friction `ξ = p_η/Q` onto every atom:
//!
//! ```text
//! dv/dt = F/m − ξ·v
//! ```
//!
//! Unlike Berendsen it is **deterministic, time-reversible** and
//! generates the *exact* canonical distribution (it has a conserved
//! pseudo-Hamiltonian).
//!
//! ## Nosé-Hoover chains
//!
//! A single Nosé-Hoover variable can fail to be ergodic for small or
//! stiff systems. A **chain** thermostats the thermostat: variable 1
//! is driven by the system, variable 2 thermostats variable 1, and so
//! on. [`NoseHoover::with_chain`] builds a chain of the requested
//! length.
//!
//! ## v1 caveat — operator splitting
//!
//! The thermostat is applied as a velocity-scaling operator each step
//! (the standard "apply ½-step thermostat, integrate, apply ½-step
//! thermostat" Trotter factorisation). The chain is propagated with a
//! **first-order** sweep over its variables. A production code uses
//! the higher-order **Suzuki-Yoshida** multi-time-step factorisation
//! of the chain propagator for better energy conservation at large
//! `dt`; that refinement is the documented future step. The
//! single-variable thermostat and the qualitative chain behaviour are
//! correct.

use crate::ensemble::Thermostat;
use crate::error::{MdError, Result};
use crate::system::System;
use crate::units::BOLTZMANN;

/// The Nosé-Hoover (chain) thermostat.
#[derive(Clone, Debug, PartialEq)]
pub struct NoseHoover {
    /// Target temperature T₀ (K).
    target: f64,
    /// Coupling time constant τ (ps) — sets the thermostat masses.
    tau: f64,
    /// Per-link thermostat "position" η.
    eta: Vec<f64>,
    /// Per-link thermostat momentum p_η.
    p_eta: Vec<f64>,
    /// Per-link thermostat mass Q.
    q: Vec<f64>,
    /// Whether the thermostat masses have been sized to the system
    /// yet (done lazily on the first `apply`, when `N_f` is known).
    sized: bool,
}

impl NoseHoover {
    /// Builds a single-variable Nosé-Hoover thermostat.
    ///
    /// * `target` — target temperature (K)
    /// * `tau` — coupling time constant (ps)
    ///
    /// # Errors
    /// [`MdError::Invalid`] on a non-finite / negative temperature or
    /// a non-positive `tau`.
    pub fn new(target: f64, tau: f64) -> Result<Self> {
        Self::with_chain(target, tau, 1)
    }

    /// Builds a Nosé-Hoover *chain* of `chain_length` thermostat
    /// variables. `chain_length == 1` is the plain Nosé-Hoover
    /// thermostat.
    ///
    /// # Errors
    /// [`MdError::Invalid`] on bad parameters or a zero chain length.
    pub fn with_chain(target: f64, tau: f64, chain_length: usize) -> Result<Self> {
        if !(target.is_finite() && target >= 0.0) {
            return Err(MdError::invalid("target", "must be finite and non-negative"));
        }
        if !(tau.is_finite() && tau > 0.0) {
            return Err(MdError::invalid("tau", "must be finite and positive"));
        }
        if chain_length == 0 {
            return Err(MdError::invalid("chain_length", "must be at least 1"));
        }
        Ok(NoseHoover {
            target,
            tau,
            eta: vec![0.0; chain_length],
            p_eta: vec![0.0; chain_length],
            q: vec![0.0; chain_length],
            sized: false,
        })
    }

    /// The chain length (1 = plain Nosé-Hoover).
    pub fn chain_length(&self) -> usize {
        self.eta.len()
    }

    /// The coupling time constant τ (ps).
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// The current thermostat "energy" — the extended-system kinetic +
    /// potential terms. Added to the physical energy it gives a
    /// conserved quantity, which is the standard correctness check for
    /// a Nosé-Hoover run.
    pub fn conserved_energy_term(&self, dof: usize) -> f64 {
        if !self.sized {
            return 0.0;
        }
        let mut e = 0.0;
        let kt = BOLTZMANN * self.target;
        for link in 0..self.eta.len() {
            // Kinetic part of the thermostat link.
            e += 0.5 * self.p_eta[link] * self.p_eta[link] / self.q[link];
            // Potential part: link 0 sees N_f·kT·η, the rest see kT·η.
            let g = if link == 0 { dof as f64 } else { 1.0 };
            e += g * kt * self.eta[link];
        }
        e
    }

    /// Sizes the thermostat masses `Q` once `N_f` is known.
    /// `Q₀ = N_f·k_B·T·τ²`, `Qᵢ = k_B·T·τ²` for the rest.
    fn size_masses(&mut self, dof: usize) {
        let kt = BOLTZMANN * self.target;
        let tau2 = self.tau * self.tau;
        for (link, q) in self.q.iter_mut().enumerate() {
            let g = if link == 0 { dof as f64 } else { 1.0 };
            *q = (g * kt * tau2).max(1e-12);
        }
        self.sized = true;
    }
}

impl Thermostat for NoseHoover {
    fn name(&self) -> &str {
        if self.chain_length() > 1 {
            "nose-hoover-chain"
        } else {
            "nose-hoover"
        }
    }

    fn target_temperature(&self) -> f64 {
        self.target
    }

    fn apply(&mut self, system: &mut System, dt: f64, constraints: usize) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        let dof = system.degrees_of_freedom(constraints);
        if dof == 0 {
            return Ok(());
        }
        if !self.sized {
            self.size_masses(dof);
        }
        let kt = BOLTZMANN * self.target;
        let n = self.eta.len();

        // --- Propagate the thermostat chain by a half step ----------
        // Twice the kinetic energy drives thermostat link 0.
        let two_ke = 2.0 * system.kinetic_energy();
        // Sweep the chain from the last link inward (standard order),
        // updating each momentum from the "force" of the link below.
        // Forces:  G₀ = 2·KE − N_f·kT ;  Gᵢ = p²_{η,i-1}/Q_{i-1} − kT.
        for link in (0..n).rev() {
            let g = if link == 0 {
                two_ke - dof as f64 * kt
            } else {
                self.p_eta[link - 1] * self.p_eta[link - 1] / self.q[link - 1] - kt
            };
            if link == n - 1 {
                // Last link has nothing above damping it.
                self.p_eta[link] += 0.5 * dt * g;
            } else {
                // Damped by the link above it.
                let damp = (-0.25 * dt * self.p_eta[link + 1] / self.q[link + 1]).exp();
                self.p_eta[link] = self.p_eta[link] * damp + 0.5 * dt * g;
                self.p_eta[link] *= damp;
            }
        }

        // Friction factor on the atomic velocities from link 0.
        let xi = self.p_eta[0] / self.q[0];
        let scale = (-xi * dt).exp();
        for v in &mut system.velocities {
            *v *= scale;
        }

        // Advance the thermostat positions η by a full step.
        for link in 0..n {
            self.eta[link] += dt * self.p_eta[link] / self.q[link];
        }

        // --- Propagate the chain by the second half step ------------
        let two_ke = 2.0 * system.kinetic_energy();
        for link in 0..n {
            let g = if link == 0 {
                two_ke - dof as f64 * kt
            } else {
                self.p_eta[link - 1] * self.p_eta[link - 1] / self.q[link - 1] - kt
            };
            if link == n - 1 {
                self.p_eta[link] += 0.5 * dt * g;
            } else {
                let damp = (-0.25 * dt * self.p_eta[link + 1] / self.q[link + 1]).exp();
                self.p_eta[link] = self.p_eta[link] * damp + 0.5 * dt * g;
                self.p_eta[link] *= damp;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::tests::hot_gas;

    #[test]
    fn rejects_bad_parameters() {
        assert!(NoseHoover::new(-1.0, 0.5).is_err());
        assert!(NoseHoover::new(300.0, 0.0).is_err());
        assert!(NoseHoover::with_chain(300.0, 0.5, 0).is_err());
    }

    #[test]
    fn chain_length_is_reported() {
        assert_eq!(NoseHoover::new(300.0, 0.5).unwrap().chain_length(), 1);
        assert_eq!(
            NoseHoover::with_chain(300.0, 0.5, 4).unwrap().chain_length(),
            4
        );
    }

    #[test]
    fn hot_system_relaxes_toward_target() {
        let mut sys = hot_gas(150, 3.0);
        assert!(sys.temperature(0) > 350.0);
        let mut thermo = NoseHoover::new(300.0, 0.2).unwrap();
        // Nosé-Hoover oscillates around the target; average over the
        // tail of a long run.
        let mut tail = Vec::new();
        for step in 0..20_000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
            if step >= 10_000 {
                tail.push(sys.temperature(0));
            }
        }
        let mean: f64 = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((mean - 300.0).abs() < 40.0, "mean tail T = {mean}");
    }

    #[test]
    fn chain_also_relaxes_toward_target() {
        let mut sys = hot_gas(150, 3.0);
        let mut thermo = NoseHoover::with_chain(300.0, 0.2, 3).unwrap();
        let mut tail = Vec::new();
        for step in 0..20_000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
            if step >= 10_000 {
                tail.push(sys.temperature(0));
            }
        }
        let mean: f64 = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((mean - 300.0).abs() < 40.0, "chain mean tail T = {mean}");
    }

    #[test]
    fn conserved_term_is_finite() {
        let mut sys = hot_gas(50, 2.0);
        let mut thermo = NoseHoover::new(300.0, 0.3).unwrap();
        for _ in 0..1000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
        }
        let term = thermo.conserved_energy_term(sys.degrees_of_freedom(0));
        assert!(term.is_finite());
    }
}
