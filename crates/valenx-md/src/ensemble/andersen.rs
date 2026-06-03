//! Stochastic thermostats — **roadmap feature 20**.
//!
//! Two stochastic thermostats that — unlike Berendsen — sample the
//! *true* canonical ensemble:
//!
//! - [`Andersen`] — the **Andersen** thermostat. Every step each atom
//!   is, with probability `ν·dt`, "collided" with the bath: its
//!   velocity is replaced by a fresh draw from the Maxwell-Boltzmann
//!   distribution at the target temperature. Between collisions the
//!   dynamics is purely Newtonian. The collision frequency `ν`
//!   controls the coupling strength. Andersen samples NVT exactly but
//!   the random reassignments interrupt the dynamics, so transport
//!   properties (diffusion) are perturbed.
//!
//! - [`VelocityRescale`] — the **Bussi-Donadio-Parrinello**
//!   velocity-rescale thermostat. A global rescale like Berendsen, but
//!   the target kinetic energy each step is drawn *stochastically*
//!   from the canonical kinetic-energy distribution. This restores the
//!   correct kinetic-energy fluctuations that Berendsen suppresses, so
//!   it *does* sample a true canonical ensemble while staying as
//!   smooth and non-disruptive as Berendsen. It is the recommended
//!   modern thermostat.
//!
//! Both use the deterministic [`crate::rng::Rng`], so a run is
//! reproducible from its seed.

use nalgebra::Vector3;

use crate::ensemble::Thermostat;
use crate::error::{MdError, Result};
use crate::rng::Rng;
use crate::system::System;
use crate::units::BOLTZMANN;

/// The Andersen stochastic-collision thermostat.
#[derive(Clone, Debug)]
pub struct Andersen {
    /// Target temperature (K).
    target: f64,
    /// Collision frequency ν (1/ps).
    nu: f64,
    /// Deterministic random source.
    rng: Rng,
}

impl Andersen {
    /// Builds an Andersen thermostat.
    ///
    /// * `target` — target temperature (K)
    /// * `nu` — collision frequency (1/ps)
    /// * `seed` — seed for the deterministic collisions
    ///
    /// # Errors
    /// [`MdError::Invalid`] on a non-finite / negative temperature or
    /// a non-positive collision frequency.
    pub fn new(target: f64, nu: f64, seed: u64) -> Result<Self> {
        if !(target.is_finite() && target >= 0.0) {
            return Err(MdError::invalid("target", "must be finite and non-negative"));
        }
        if !(nu.is_finite() && nu > 0.0) {
            return Err(MdError::invalid("nu", "must be finite and positive"));
        }
        Ok(Andersen {
            target,
            nu,
            rng: Rng::new(seed),
        })
    }

    /// The collision frequency ν (1/ps).
    pub fn nu(&self) -> f64 {
        self.nu
    }

    /// Draws a Maxwell-Boltzmann velocity for an atom of the given
    /// mass at the target temperature. Each component is normal with
    /// standard deviation `√(k_B·T/m)`.
    fn maxwell_velocity(&mut self, mass: f64) -> Vector3<f64> {
        let sigma = (BOLTZMANN * self.target / mass).sqrt();
        Vector3::new(
            sigma * self.rng.normal(),
            sigma * self.rng.normal(),
            sigma * self.rng.normal(),
        )
    }
}

impl Thermostat for Andersen {
    fn name(&self) -> &str {
        "andersen"
    }

    fn target_temperature(&self) -> f64 {
        self.target
    }

    fn apply(&mut self, system: &mut System, dt: f64, _constraints: usize) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        let p_collide = (self.nu * dt).min(1.0);
        for i in 0..system.len() {
            if self.rng.bernoulli(p_collide) {
                let mass = system.topology.atoms[i].mass;
                system.velocities[i] = self.maxwell_velocity(mass);
            }
        }
        Ok(())
    }
}

/// The Bussi-Donadio-Parrinello velocity-rescale thermostat.
#[derive(Clone, Debug)]
pub struct VelocityRescale {
    /// Target temperature (K).
    target: f64,
    /// Coupling time constant τ (ps).
    tau: f64,
    /// Deterministic random source.
    rng: Rng,
}

impl VelocityRescale {
    /// Builds a velocity-rescale (Bussi) thermostat.
    ///
    /// * `target` — target temperature (K)
    /// * `tau` — coupling time constant (ps)
    /// * `seed` — seed for the stochastic target kinetic energy
    ///
    /// # Errors
    /// [`MdError::Invalid`] on a non-finite / negative temperature or
    /// a non-positive `tau`.
    pub fn new(target: f64, tau: f64, seed: u64) -> Result<Self> {
        if !(target.is_finite() && target >= 0.0) {
            return Err(MdError::invalid("target", "must be finite and non-negative"));
        }
        if !(tau.is_finite() && tau > 0.0) {
            return Err(MdError::invalid("tau", "must be finite and positive"));
        }
        Ok(VelocityRescale {
            target,
            tau,
            rng: Rng::new(seed),
        })
    }

    /// The coupling time constant τ (ps).
    pub fn tau(&self) -> f64 {
        self.tau
    }
}

impl Thermostat for VelocityRescale {
    fn name(&self) -> &str {
        "velocity-rescale"
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
        let ke = system.kinetic_energy();
        if ke <= 1e-12 {
            return Ok(());
        }
        // Target average kinetic energy at the set temperature.
        let ke_target = 0.5 * dof as f64 * BOLTZMANN * self.target;
        // Bussi-Donadio-Parrinello stochastic-rescale update of the
        // kinetic energy. With c = exp(-dt/τ):
        //   KE' = KE + (1-c)·(K̄·R²/Nf − KE)
        //             + 2·√(c·(1-c)·KE·K̄/Nf)·R1
        // where R1 is one standard normal and R² is the sum of Nf−1
        // squared normals (a χ² draw).
        let c = (-dt / self.tau).exp();
        let r1 = self.rng.normal();
        let sum_noises_sq = self.rng.chi_squared(dof.saturating_sub(1));
        let ke_per_dof = ke_target / dof as f64;
        let ke_new = ke
            + (1.0 - c) * (ke_per_dof * (r1 * r1 + sum_noises_sq) - ke)
            + 2.0 * (c * (1.0 - c) * ke * ke_per_dof).max(0.0).sqrt() * r1;
        let ke_new = ke_new.max(1e-12);
        let lambda = (ke_new / ke).sqrt();
        for v in &mut system.velocities {
            *v *= lambda;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::tests::hot_gas;

    #[test]
    fn andersen_rejects_bad_parameters() {
        assert!(Andersen::new(-1.0, 1.0, 1).is_err());
        assert!(Andersen::new(300.0, 0.0, 1).is_err());
    }

    #[test]
    fn andersen_drives_to_target_temperature() {
        let mut sys = hot_gas(200, 3.0);
        let mut thermo = Andersen::new(250.0, 10.0, 42).unwrap();
        for _ in 0..3000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
        }
        let t = sys.temperature(0);
        assert!((t - 250.0).abs() < 30.0, "T = {t}");
    }

    #[test]
    fn velocity_rescale_rejects_bad_parameters() {
        assert!(VelocityRescale::new(-1.0, 0.1, 1).is_err());
        assert!(VelocityRescale::new(300.0, -0.1, 1).is_err());
    }

    #[test]
    fn velocity_rescale_drives_to_target_temperature() {
        let mut sys = hot_gas(200, 3.0);
        let mut thermo = VelocityRescale::new(300.0, 0.1, 7).unwrap();
        for _ in 0..3000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
        }
        let t = sys.temperature(0);
        assert!((t - 300.0).abs() < 25.0, "T = {t}");
    }

    #[test]
    fn velocity_rescale_is_deterministic() {
        let mut a = hot_gas(50, 2.0);
        let mut b = hot_gas(50, 2.0);
        let mut ta = VelocityRescale::new(300.0, 0.2, 99).unwrap();
        let mut tb = VelocityRescale::new(300.0, 0.2, 99).unwrap();
        for _ in 0..200 {
            ta.apply(&mut a, 0.002, 0).unwrap();
            tb.apply(&mut b, 0.002, 0).unwrap();
        }
        for (va, vb) in a.velocities.iter().zip(&b.velocities) {
            assert!((va - vb).norm() < 1e-12);
        }
    }
}
