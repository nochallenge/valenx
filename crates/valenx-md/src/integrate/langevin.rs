//! Langevin / Brownian dynamics — **roadmap feature 18**.
//!
//! Where velocity-Verlet integrates pure Newtonian mechanics, the
//! **Langevin** equation adds two extra forces that together model an
//! implicit heat bath:
//!
//! ```text
//! m·a = F_conservative − γ·m·v + R(t)
//! ```
//!
//! - `−γ·m·v` — a viscous **friction** drag, `γ` the collision
//!   frequency (1/ps).
//! - `R(t)` — a **random** force, Gaussian white noise whose variance
//!   is fixed by the fluctuation-dissipation theorem
//!   `⟨R²⟩ = 2·γ·m·k_B·T/dt` so the integrator samples the canonical
//!   (NVT) ensemble at temperature `T` *without a separate
//!   thermostat*.
//!
//! ## Two regimes, one type
//!
//! - **Langevin dynamics** (the default) keeps the inertial `m·a`
//!   term; it uses the BAOAB-style velocity update — half-kick,
//!   drift, Ornstein-Uhlenbeck friction+noise, drift, half-kick —
//!   which is the accurate modern Langevin splitting.
//! - **Brownian (overdamped) dynamics** is the `γ → ∞` limit where
//!   inertia is negligible and the position update is
//!   `r(t+dt) = r(t) + dt·F/(γ·m) + √(2·k_B·T·dt/(γ·m))·ξ`.
//!   Select it with [`Langevin::brownian`].
//!
//! The generator is the crate's deterministic [`crate::rng::Rng`], so
//! a Langevin run is fully reproducible from its seed.

use nalgebra::Vector3;

use crate::bonded::EnergyForce;
use crate::error::{MdError, Result};
use crate::integrate::Integrator;
use crate::rng::Rng;
use crate::system::System;
use crate::units::BOLTZMANN;

/// Whether the integrator runs full Langevin or overdamped Brownian
/// dynamics.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Regime {
    /// Inertial Langevin dynamics (BAOAB splitting).
    Langevin,
    /// Overdamped Brownian dynamics.
    Brownian,
}

/// The Langevin / Brownian-dynamics integrator.
#[derive(Clone, Debug)]
pub struct Langevin {
    /// Time step (ps).
    dt: f64,
    /// Collision / friction frequency γ (1/ps).
    gamma: f64,
    /// Target temperature (K).
    temperature: f64,
    /// Integration regime.
    regime: Regime,
    /// The deterministic noise source.
    rng: Rng,
    /// Cached conservative forces (Langevin regime only).
    last_forces: Option<Vec<Vector3<f64>>>,
}

impl Langevin {
    /// Builds a full (inertial) Langevin integrator.
    ///
    /// * `dt` — time step (ps)
    /// * `gamma` — friction frequency (1/ps); a few 1/ps is typical
    /// * `temperature` — target temperature (K)
    /// * `seed` — seed for the deterministic random force
    ///
    /// # Errors
    /// [`MdError::Invalid`] on any non-positive / non-finite argument.
    pub fn new(dt: f64, gamma: f64, temperature: f64, seed: u64) -> Result<Self> {
        Self::build(dt, gamma, temperature, seed, Regime::Langevin)
    }

    /// Builds an overdamped Brownian-dynamics integrator. Same
    /// arguments as [`new`](Self::new); `gamma` here is the inverse
    /// mobility (drag), and inertia is dropped.
    ///
    /// # Errors
    /// [`MdError::Invalid`] on any non-positive / non-finite argument.
    pub fn brownian(dt: f64, gamma: f64, temperature: f64, seed: u64) -> Result<Self> {
        Self::build(dt, gamma, temperature, seed, Regime::Brownian)
    }

    fn build(dt: f64, gamma: f64, temperature: f64, seed: u64, regime: Regime) -> Result<Self> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        if !(gamma.is_finite() && gamma > 0.0) {
            return Err(MdError::invalid("gamma", "must be finite and positive"));
        }
        if !(temperature.is_finite() && temperature >= 0.0) {
            return Err(MdError::invalid(
                "temperature",
                "must be finite and non-negative",
            ));
        }
        Ok(Langevin {
            dt,
            gamma,
            temperature,
            regime,
            rng: Rng::new(seed),
            last_forces: None,
        })
    }

    /// The target temperature (K).
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// The friction frequency γ (1/ps).
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// Discards cached conservative forces.
    pub fn reset(&mut self) {
        self.last_forces = None;
    }

    /// A vector of three independent standard-normal draws.
    fn normal3(&mut self) -> Vector3<f64> {
        Vector3::new(self.rng.normal(), self.rng.normal(), self.rng.normal())
    }

    /// One BAOAB Langevin step.
    fn step_langevin(
        &mut self,
        system: &mut System,
        force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
    ) -> Result<EnergyForce> {
        let n = system.len();
        let dt = self.dt;
        let half_dt = 0.5 * dt;
        // Ornstein-Uhlenbeck factors for the O (friction+noise) step.
        let c1 = (-self.gamma * dt).exp();
        let c2 = (1.0 - c1 * c1).sqrt();

        let forces_t = match self.last_forces.take() {
            Some(f) if f.len() == n => f,
            _ => force_fn(system)?.forces,
        };

        // B: half kick.
        for ((v, atom), force) in system
            .velocities
            .iter_mut()
            .zip(&system.topology.atoms)
            .zip(&forces_t)
        {
            *v += half_dt * force / atom.mass;
        }
        // A: half drift.
        for (p, v) in system.positions.iter_mut().zip(&system.velocities) {
            *p += half_dt * v;
        }
        // O: friction + noise (Ornstein-Uhlenbeck on the velocity).
        // An index loop: each iteration draws fresh noise, which needs
        // `&mut self`, so a `system` iterator and `self.normal3()`
        // cannot be held at once.
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            let mass = system.topology.atoms[i].mass;
            let sigma_v = (BOLTZMANN * self.temperature / mass).sqrt();
            let xi = self.normal3();
            system.velocities[i] = c1 * system.velocities[i] + c2 * sigma_v * xi;
        }
        // A: half drift.
        for (p, v) in system.positions.iter_mut().zip(&system.velocities) {
            *p += half_dt * v;
        }
        // B: half kick with the new forces.
        let ef_new = force_fn(system)?;
        if ef_new.forces.len() != n {
            return Err(MdError::dimension(
                "force evaluation returned the wrong number of forces",
            ));
        }
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

    /// One overdamped Brownian step.
    fn step_brownian(
        &mut self,
        system: &mut System,
        force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
    ) -> Result<EnergyForce> {
        let n = system.len();
        let dt = self.dt;
        let ef = force_fn(system)?;
        if ef.forces.len() != n {
            return Err(MdError::dimension(
                "force evaluation returned the wrong number of forces",
            ));
        }
        // An index loop: each iteration draws fresh noise via
        // `&mut self`, which precludes holding a `system` iterator.
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            let mass = system.topology.atoms[i].mass;
            // Mobility μ = 1/(γ·m); drift = μ·F·dt.
            let mobility = 1.0 / (self.gamma * mass);
            let drift = mobility * ef.forces[i] * dt;
            // Diffusive step: √(2·D·dt)·ξ with D = μ·k_B·T.
            let d_coeff = mobility * BOLTZMANN * self.temperature;
            let noise = (2.0 * d_coeff * dt).sqrt() * self.normal3();
            system.positions[i] += drift + noise;
            // Report an instantaneous velocity estimate (drift/dt).
            system.velocities[i] = (drift + noise) / dt;
        }
        Ok(ef)
    }
}

impl Integrator for Langevin {
    fn name(&self) -> &str {
        match self.regime {
            Regime::Langevin => "langevin",
            Regime::Brownian => "brownian",
        }
    }

    fn dt(&self) -> f64 {
        self.dt
    }

    fn step(
        &mut self,
        system: &mut System,
        force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
    ) -> Result<EnergyForce> {
        match self.regime {
            Regime::Langevin => self.step_langevin(system, force_fn),
            Regime::Brownian => self.step_brownian(system, force_fn),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};

    /// Many non-interacting particles — an ideal gas. A Langevin
    /// integrator should drive them to the target temperature.
    fn ideal_gas(n: usize) -> System {
        let mut top = Topology::new();
        for _ in 0..n {
            top.push_atom(Atom::new("A", 16.0, 0.0).unwrap());
        }
        let pos = (0..n)
            .map(|i| Vector3::new(i as f64 * 0.4, 0.0, 0.0))
            .collect();
        System::new(top, pos).unwrap()
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(Langevin::new(0.0, 1.0, 300.0, 1).is_err());
        assert!(Langevin::new(0.001, -1.0, 300.0, 1).is_err());
        assert!(Langevin::new(0.001, 1.0, -5.0, 1).is_err());
    }

    #[test]
    fn langevin_thermostats_to_target_temperature() {
        let mut sys = ideal_gas(200);
        // Start cold (zero velocity); Langevin should heat to 300 K.
        let target = 300.0;
        let mut integ = Langevin::new(0.002, 5.0, target, 12345).unwrap();
        let mut zero = |s: &System| Ok(EnergyForce::zeros(s.len()));
        for _ in 0..4000 {
            integ.step(&mut sys, &mut zero).unwrap();
        }
        let t = sys.temperature(0);
        // Should equilibrate near the target (no conservative forces,
        // so the only steady state is the bath temperature).
        assert!((t - target).abs() < 30.0, "equilibrated T = {t}");
    }

    #[test]
    fn brownian_diffuses_with_einstein_coefficient() {
        // Overdamped: mean-squared displacement grows as 6·D·t.
        let mut sys = ideal_gas(300);
        let gamma = 50.0;
        let temp = 300.0;
        let dt = 0.001;
        let mut integ = Langevin::brownian(dt, gamma, temp, 999).unwrap();
        let start: Vec<Vector3<f64>> = sys.positions.clone();
        let mut zero = |s: &System| Ok(EnergyForce::zeros(s.len()));
        let steps = 2000;
        for _ in 0..steps {
            integ.step(&mut sys, &mut zero).unwrap();
        }
        let msd: f64 = sys
            .positions
            .iter()
            .zip(&start)
            .map(|(p, s)| (p - s).norm_squared())
            .sum::<f64>()
            / sys.len() as f64;
        // Expected D = k_B·T/(γ·m); MSD = 6·D·t.
        let mass = 16.0;
        let d_expected = BOLTZMANN * temp / (gamma * mass);
        let t_total = dt * steps as f64;
        let msd_expected = 6.0 * d_expected * t_total;
        assert!(
            (msd - msd_expected).abs() < 0.5 * msd_expected,
            "MSD {msd} vs expected {msd_expected}"
        );
    }

    #[test]
    fn is_deterministic_for_a_seed() {
        let mut a = ideal_gas(10);
        let mut b = ideal_gas(10);
        let mut ia = Langevin::new(0.001, 2.0, 250.0, 7).unwrap();
        let mut ib = Langevin::new(0.001, 2.0, 250.0, 7).unwrap();
        let mut zero = |s: &System| Ok(EnergyForce::zeros(s.len()));
        for _ in 0..100 {
            ia.step(&mut a, &mut zero).unwrap();
            ib.step(&mut b, &mut zero).unwrap();
        }
        for (pa, pb) in a.positions.iter().zip(&b.positions) {
            assert!((pa - pb).norm() < 1e-12);
        }
    }
}
