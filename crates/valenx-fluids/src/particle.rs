//! The particle and the particle system.
//!
//! Each [`Particle`] is a Lagrangian sample of the fluid carrying a position and
//! velocity (the state advanced by integration) plus the per-step scalars the
//! SPH solver fills in: the smoothed **density** estimate and the **pressure**
//! derived from it via an equation of state. Force / acceleration is
//! accumulated each step and consumed by the integrator.
//!
//! A [`ParticleSystem`] is just a `Vec<Particle>` with finite-value-checked
//! mutators, so a `NaN` position can never silently enter the simulation.

use nalgebra::Vector3;

use crate::error::FluidError;

/// One SPH fluid particle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
    /// World-space position (m).
    pub position: Vector3<f64>,
    /// Velocity (m/s).
    pub velocity: Vector3<f64>,
    /// Acceleration accumulated for the current step (m/s²); reset each step.
    pub acceleration: Vector3<f64>,
    /// Smoothed density estimate at this particle (kg/m³); filled by the solver.
    pub density: f64,
    /// Pressure from the equation of state (Pa); filled by the solver.
    pub pressure: f64,
}

impl Particle {
    /// A particle at `position` with `velocity`, zero accumulated acceleration
    /// and zeroed density/pressure.
    #[must_use]
    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>) -> Self {
        Self {
            position,
            velocity,
            acceleration: Vector3::zeros(),
            density: 0.0,
            pressure: 0.0,
        }
    }

    /// A particle at rest (zero velocity) at `position`.
    #[must_use]
    pub fn at(position: Vector3<f64>) -> Self {
        Self::new(position, Vector3::zeros())
    }
}

/// A collection of fluid particles.
///
/// Particles are stored in a flat `Vec`; the solver indexes them by their
/// position in this vector, which is also what the spatial hash in
/// [`crate::grid`] returns. Mutators validate that no non-finite value enters.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ParticleSystem {
    particles: Vec<Particle>,
}

impl ParticleSystem {
    /// An empty particle system.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a system from an iterator of positions, all at rest.
    ///
    /// # Errors
    /// [`FluidError::NonFinite`] if any position has a non-finite component.
    pub fn from_positions<I>(positions: I) -> Result<Self, FluidError>
    where
        I: IntoIterator<Item = Vector3<f64>>,
    {
        let mut sys = Self::new();
        for p in positions {
            sys.push(Particle::at(p))?;
        }
        Ok(sys)
    }

    /// Add a particle, rejecting non-finite position or velocity.
    ///
    /// # Errors
    /// [`FluidError::NonFinite`] if the particle's position or velocity has a
    /// non-finite component.
    pub fn push(&mut self, p: Particle) -> Result<(), FluidError> {
        if !p.position.iter().all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("particle position".into()));
        }
        if !p.velocity.iter().all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("particle velocity".into()));
        }
        self.particles.push(p);
        Ok(())
    }

    /// The number of particles.
    #[must_use]
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// Whether the system has no particles.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    /// A read-only view of the particles.
    #[must_use]
    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }

    /// A mutable view of the particles (used internally by the solver).
    #[must_use]
    pub fn particles_mut(&mut self) -> &mut [Particle] {
        &mut self.particles
    }

    /// The total linear momentum `Σ m·v` for a uniform particle mass `m`.
    ///
    /// (All particles in a single fluid share one mass in this solver, so the
    /// caller passes it in rather than storing it per particle.)
    #[must_use]
    pub fn total_momentum(&self, mass: f64) -> Vector3<f64> {
        let v: Vector3<f64> = self.particles.iter().map(|p| p.velocity).sum();
        v * mass
    }

    /// The centre of mass (the mean position; mass cancels for a uniform mass).
    ///
    /// Returns the zero vector for an empty system.
    #[must_use]
    pub fn center_of_mass(&self) -> Vector3<f64> {
        if self.particles.is_empty() {
            return Vector3::zeros();
        }
        let sum: Vector3<f64> = self.particles.iter().map(|p| p.position).sum();
        sum / self.particles.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_non_finite_position() {
        let mut sys = ParticleSystem::new();
        assert!(sys
            .push(Particle::at(Vector3::new(f64::NAN, 0.0, 0.0)))
            .is_err());
        assert!(sys
            .push(Particle::at(Vector3::new(0.0, f64::INFINITY, 0.0)))
            .is_err());
        assert!(sys.push(Particle::at(Vector3::new(0.0, 0.0, 0.0))).is_ok());
        assert_eq!(sys.len(), 1);
    }

    #[test]
    fn rejects_non_finite_velocity() {
        let mut sys = ParticleSystem::new();
        let p = Particle::new(Vector3::zeros(), Vector3::new(f64::NAN, 0.0, 0.0));
        assert!(sys.push(p).is_err());
    }

    #[test]
    fn from_positions_builds_resting_particles() {
        let sys = ParticleSystem::from_positions([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ])
        .unwrap();
        assert_eq!(sys.len(), 2);
        assert!(sys.particles().iter().all(|p| p.velocity.norm() == 0.0));
    }

    #[test]
    fn momentum_and_center_of_mass() {
        let mut sys = ParticleSystem::new();
        sys.push(Particle::new(
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
        ))
        .unwrap();
        sys.push(Particle::new(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(-2.0, 0.0, 0.0),
        ))
        .unwrap();
        // Equal and opposite velocities ⇒ zero net momentum.
        assert!(sys.total_momentum(1.5).norm() < 1e-15);
        // Symmetric positions ⇒ centre of mass at the origin.
        assert!(sys.center_of_mass().norm() < 1e-15);
    }
}
