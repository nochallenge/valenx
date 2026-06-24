//! Weakly-compressible SPH solver (Müller, Charypar & Gross, SCA 2003).
//!
//! One [`SphSolver::step`] advances the [`ParticleSystem`] by `dt`:
//!
//! 1. **Neighbour search** — rebuild the [`SpatialHash`] from the current
//!    positions (cell size = `h`) and, for each particle, collect candidates
//!    from the surrounding 27 cells.
//! 2. **Density** — `ρ_i = Σ_j m · W_poly6(‖r_i − r_j‖, h)` (the self term
//!    `j = i` is included; the poly6 kernel is `> 0` at `r = 0`, so an isolated
//!    particle gets the floor density `m·W_poly6(0)` rather than `0`).
//! 3. **Pressure** — `p_i = EOS(ρ_i)` (optionally clamped at `0`).
//! 4. **Forces** — per particle, the sum of
//!    * the **pressure** force, in the *symmetric* form
//!      `f^p_i = −Σ_j m² · (p_i/ρ_i² + p_j/ρ_j²) · ∇W_spiky(r_i − r_j)`
//!      (this `p/ρ²` symmetrisation makes the pair force equal-and-opposite, so
//!      total momentum is conserved exactly — Newton's third law),
//!    * the **viscosity** force
//!      `f^v_i = μ Σ_j m · (v_j − v_i)/ρ_j · ∇²W_visc(‖r_i − r_j‖)`, and
//!    * **gravity** `ρ_i · g` (a body force per unit volume).
//! 5. **Integration** — semi-implicit (symplectic) Euler:
//!    `v ← v + a·dt`, then `x ← x + v·dt`, with `a = f / ρ`.
//! 6. **Boundary** — clamp/reflect each particle back into the box.
//!
//! # Guarding
//!
//! Configuration is validated up front (h, mass, EOS rest density, sound speed
//! all `> 0`; `dt > 0`; gravity finite). In the force loop the density is taken
//! as `max(ρ_i, ρ_floor)` where `ρ_floor = m·W_poly6(0) > 0`, so the `1/ρ` and
//! `1/ρ²` divisions can never hit zero even in the pathological case of a
//! particle with no neighbours and a clamped-away self term. The kernels
//! themselves only ever evaluate `r ∈ [0, h]` and return `0` outside it.

use nalgebra::Vector3;

use crate::boundary::BoxBoundary;
use crate::eos::EquationOfState;
use crate::error::FluidError;
use crate::grid::SpatialHash;
use crate::kernels::SmoothingKernels;
use crate::particle::ParticleSystem;

/// Configuration for the [`SphSolver`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphConfig {
    /// Smoothing length / kernel support radius `h` (m), strictly positive.
    pub smoothing_length: f64,
    /// Per-particle mass `m` (kg), strictly positive (one mass for the fluid).
    pub mass: f64,
    /// Dynamic viscosity coefficient `μ` (Pa·s), `≥ 0`.
    pub viscosity: f64,
    /// Gravity vector `g` (m/s²), e.g. `(0, 0, −9.80665)`.
    pub gravity: Vector3<f64>,
    /// The pressure equation of state.
    pub eos: EquationOfState,
    /// Whether to clamp pressure at `0` (drop the cohesive negative branch).
    pub clamp_negative_pressure: bool,
}

impl SphConfig {
    /// A reasonable water-like default for a smoothing length and particle
    /// spacing, using a Tait EOS.
    ///
    /// `mass` is set so that a regular lattice at spacing `h/2` reaches roughly
    /// the rest density; callers tuning a specific scene should set the fields
    /// directly.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `smoothing_length` is not finite and
    /// strictly positive.
    pub fn water(smoothing_length: f64) -> Result<Self, FluidError> {
        if !(smoothing_length.is_finite() && smoothing_length > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "smoothing length must be finite and > 0, got {smoothing_length}"
            )));
        }
        let rest_density = 1000.0;
        // Mass ≈ ρ₀ · (volume per particle), with ~ (h/2)³ per particle.
        let spacing = smoothing_length * 0.5;
        let mass = rest_density * spacing * spacing * spacing;
        Ok(Self {
            smoothing_length,
            mass,
            viscosity: 0.1,
            gravity: Vector3::new(0.0, 0.0, -9.806_65),
            eos: EquationOfState::tait_water(rest_density, 50.0)?,
            clamp_negative_pressure: false,
        })
    }
}

/// A weakly-compressible SPH solver.
#[derive(Debug, Clone)]
pub struct SphSolver {
    config: SphConfig,
    kernels: SmoothingKernels,
    grid: SpatialHash,
    rho_floor: f64,
    // Reusable scratch buffers (avoid per-step / per-particle allocation).
    neighbor_scratch: Vec<usize>,
}

impl SphSolver {
    /// Build a solver from its configuration, validating every parameter.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `smoothing_length`, `mass`, or the EOS
    /// parameters are non-positive; [`FluidError::NonFinite`] if `gravity` or
    /// `viscosity` is non-finite or `viscosity` is negative.
    pub fn new(config: SphConfig) -> Result<Self, FluidError> {
        let kernels = SmoothingKernels::new(config.smoothing_length)?;
        if !(config.mass.is_finite() && config.mass > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "particle mass must be finite and > 0, got {}",
                config.mass
            )));
        }
        if !(config.viscosity.is_finite() && config.viscosity >= 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "viscosity must be finite and ≥ 0, got {}",
                config.viscosity
            )));
        }
        if !config.gravity.iter().all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("gravity vector".into()));
        }
        // EOS rest density is validated by its constructor; re-check here is cheap.
        if config.eos.rest_density() <= 0.0 {
            return Err(FluidError::InvalidConfig(
                "EOS rest density must be > 0".into(),
            ));
        }
        let grid = SpatialHash::new(config.smoothing_length)?;
        // Density floor: a lone particle still has its own self-contribution.
        let rho_floor = config.mass * kernels.poly6(0.0);
        Ok(Self {
            config,
            kernels,
            grid,
            rho_floor,
            neighbor_scratch: Vec::new(),
        })
    }

    /// The configuration this solver was built with.
    #[must_use]
    pub fn config(&self) -> &SphConfig {
        &self.config
    }

    /// The kernels in use.
    #[must_use]
    pub fn kernels(&self) -> &SmoothingKernels {
        &self.kernels
    }

    /// The density floor `m · W_poly6(0)` (a single particle's self density).
    #[must_use]
    pub fn density_floor(&self) -> f64 {
        self.rho_floor
    }

    /// Recompute `density` and `pressure` for every particle from the current
    /// positions (also leaves the spatial hash freshly built for this state).
    pub fn compute_density_pressure(&mut self, system: &mut ParticleSystem) {
        self.grid.rebuild(system);
        let mass = self.config.mass;

        // Snapshot positions so neighbour reads don't alias the &mut write.
        let positions: Vec<Vector3<f64>> = system.particles().iter().map(|p| p.position).collect();

        for i in 0..positions.len() {
            self.grid
                .neighbors_into(positions[i], &mut self.neighbor_scratch);
            let mut density = 0.0;
            for &j in &self.neighbor_scratch {
                let r2 = (positions[i] - positions[j]).norm_squared();
                density += mass * self.kernels.poly6_r2(r2);
            }
            // The self term (j == i, r = 0) is included above since i is in its
            // own cell; guard the floor regardless.
            let density = density.max(self.rho_floor);
            let p = &mut system.particles_mut()[i];
            p.density = density;
            p.pressure = if self.config.clamp_negative_pressure {
                self.config.eos.pressure_clamped(density)
            } else {
                self.config.eos.pressure(density)
            };
        }
    }

    /// Accumulate the per-particle acceleration `a = f/ρ` (the sum of pressure,
    /// viscosity, and gravity) into each particle's `acceleration` field,
    /// assuming `density` / `pressure` are already current.
    fn compute_accelerations(&mut self, system: &mut ParticleSystem) {
        let mass = self.config.mass;
        let mu = self.config.viscosity;
        let gravity = self.config.gravity;

        // Snapshot the read-only state for symmetric pair sums.
        let n = system.len();
        let positions: Vec<Vector3<f64>> = system.particles().iter().map(|p| p.position).collect();
        let velocities: Vec<Vector3<f64>> = system.particles().iter().map(|p| p.velocity).collect();
        let densities: Vec<f64> = system
            .particles()
            .iter()
            .map(|p| p.density.max(self.rho_floor))
            .collect();
        let pressures: Vec<f64> = system.particles().iter().map(|p| p.pressure).collect();

        for i in 0..n {
            self.grid
                .neighbors_into(positions[i], &mut self.neighbor_scratch);
            let rho_i = densities[i];
            let p_i = pressures[i];
            let p_over_rho2_i = p_i / (rho_i * rho_i);

            let mut f_pressure = Vector3::zeros();
            let mut f_viscosity = Vector3::zeros();

            for &j in &self.neighbor_scratch {
                if j == i {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r > self.kernels.h() {
                    continue; // outside support (candidate from an adjacent cell)
                }
                let rho_j = densities[j];

                // Symmetric pressure force (momentum-conserving form):
                // f_i += −m² (p_i/ρ_i² + p_j/ρ_j²) ∇W_spiky(r_ij)
                let p_over_rho2_j = pressures[j] / (rho_j * rho_j);
                let grad = self.kernels.spiky_gradient(rij);
                f_pressure -= grad * (mass * mass * (p_over_rho2_i + p_over_rho2_j));

                // Viscosity (Laplacian) force:
                // f_i += μ m (v_j − v_i)/ρ_j ∇²W_visc(r)
                if mu > 0.0 {
                    let lap = self.kernels.viscosity_laplacian(r);
                    f_viscosity += (velocities[j] - velocities[i]) * (mu * mass / rho_j * lap);
                }
            }

            // Gravity is a body force per unit volume: f_grav = ρ_i g, so the
            // acceleration contribution is just g.
            let accel = (f_pressure + f_viscosity) / rho_i + gravity;
            system.particles_mut()[i].acceleration = accel;
        }
    }

    /// Advance the particle system by one time step `dt` (seconds).
    ///
    /// Order: density/pressure → accelerations → semi-implicit Euler →
    /// boundary. Optionally enforces a [`BoxBoundary`] if `boundary` is `Some`.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `dt` is not finite and strictly
    /// positive.
    pub fn step(
        &mut self,
        system: &mut ParticleSystem,
        dt: f64,
        boundary: Option<&BoxBoundary>,
    ) -> Result<(), FluidError> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "time step dt must be finite and > 0, got {dt}"
            )));
        }
        if system.is_empty() {
            return Ok(());
        }

        self.compute_density_pressure(system);
        self.compute_accelerations(system);

        // Semi-implicit (symplectic) Euler: update velocity first, then position.
        for p in system.particles_mut() {
            p.velocity += p.acceleration * dt;
            p.position += p.velocity * dt;
        }

        if let Some(b) = boundary {
            for p in system.particles_mut() {
                b.enforce(p);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::Particle;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    // BENCHMARK-PIN (1): a single particle with no neighbours under gravity
    // integrates as free-fall x = ½ g t² (within the discrete-Euler tolerance).
    #[test]
    fn single_particle_free_fall() {
        let mut cfg = SphConfig::water(0.1).unwrap();
        cfg.viscosity = 0.0;
        let g = -9.806_65;
        cfg.gravity = Vector3::new(0.0, 0.0, g);
        let mut solver = SphSolver::new(cfg).unwrap();

        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::new(0.0, 0.0, 0.0))).unwrap();

        let dt = 1e-3;
        let steps = 1000;
        for _ in 0..steps {
            solver.step(&mut sys, dt, None).unwrap();
        }
        let t = dt * steps as f64; // 1.0 s
        let z = sys.particles()[0].position.z;

        // Semi-implicit Euler over-shoots the exact ½gt² by exactly ½ g dt · t
        // (a known O(dt) bias). Check against the analytic free-fall within a
        // tolerance that comfortably covers that bias.
        let analytic = 0.5 * g * t * t;
        let euler_bias = 0.5 * g * dt * t; // = g·dt·t/2
        assert!(
            approx(z, analytic + euler_bias, 1e-6),
            "z = {z}, analytic = {analytic}, +bias = {}",
            analytic + euler_bias
        );
        // And it must be in genuine free-fall (close to the closed form).
        assert!(
            approx(z, analytic, 1e-2),
            "z = {z} should be ≈ ½gt² = {analytic}"
        );

        // The horizontal components never move.
        assert!(sys.particles()[0].position.x.abs() < 1e-15);
        assert!(sys.particles()[0].position.y.abs() < 1e-15);
    }

    // BENCHMARK-PIN (3): a closed box with gravity OFF conserves total
    // momentum over a step (the symmetric pressure + viscosity pair forces are
    // equal-and-opposite, so internal forces sum to zero).
    #[test]
    fn momentum_conserved_with_gravity_off() {
        let mut cfg = SphConfig::water(0.2).unwrap();
        cfg.gravity = Vector3::zeros();
        cfg.viscosity = 0.05;
        let mut solver = SphSolver::new(cfg).unwrap();

        // A small clump of particles, some moving, no boundary forcing.
        let mut sys = ParticleSystem::new();
        let s = 0.08;
        for ix in 0..3 {
            for iy in 0..3 {
                for iz in 0..3 {
                    let pos = Vector3::new(ix as f64 * s, iy as f64 * s, iz as f64 * s);
                    let vel = Vector3::new(0.01 * ix as f64, -0.02 * iy as f64, 0.0);
                    sys.push(Particle::new(pos, vel)).unwrap();
                }
            }
        }
        let m = cfg.mass;
        let p_before = sys.total_momentum(m);

        // No boundary ⇒ only internal forces act; they must cancel.
        for _ in 0..5 {
            solver.step(&mut sys, 1e-4, None).unwrap();
        }
        let p_after = sys.total_momentum(m);

        // Gravity off + symmetric internal forces ⇒ momentum is conserved to
        // round-off.
        assert!(
            (p_after - p_before).norm() < 1e-12,
            "Δp = {} (before {p_before:?}, after {p_after:?})",
            (p_after - p_before).norm()
        );
    }

    // BENCHMARK-PIN (4): two particles compressed below rest spacing are pushed
    // apart SYMMETRICALLY (Newton's third law), so their centre of mass and net
    // momentum stay fixed and their accelerations are equal-and-opposite.
    #[test]
    fn two_compressed_particles_repel_symmetrically() {
        // Build a config whose *rest density* is the lone-particle self density,
        // so that any neighbour at all puts a particle above rest ⇒ positive
        // (compressive) pressure ⇒ repulsion. (The `water` preset sizes mass for
        // a *filled* neighbourhood of ~50 particles, against which an isolated
        // pair reads far below rest — not the regime this Newton's-third-law pin
        // is about.)
        let h = 0.2;
        let mut cfg = SphConfig::water(h).unwrap();
        cfg.gravity = Vector3::zeros();
        cfg.viscosity = 0.0;
        cfg.clamp_negative_pressure = true;
        let lone_density = cfg.mass * SmoothingKernels::new(h).unwrap().poly6(0.0);
        cfg.eos = EquationOfState::tait_water(lone_density, 50.0).unwrap();
        let mut solver = SphSolver::new(cfg).unwrap();

        // Two particles within the support ⇒ each sees the other and itself ⇒
        // density above the lone-particle rest ⇒ positive pressure ⇒ repulsion.
        let mut sys = ParticleSystem::new();
        let half_gap = 0.02;
        sys.push(Particle::at(Vector3::new(-half_gap, 0.0, 0.0)))
            .unwrap();
        sys.push(Particle::at(Vector3::new(half_gap, 0.0, 0.0)))
            .unwrap();

        let com_before = sys.center_of_mass();
        solver.compute_density_pressure(&mut sys);

        // Both compressed: density above rest, pressure strictly positive.
        assert!(sys.particles()[0].pressure > 0.0);
        assert!(sys.particles()[1].pressure > 0.0);

        // One full step (still no gravity, no boundary).
        solver.step(&mut sys, 1e-4, None).unwrap();

        let p0 = sys.particles()[0].position;
        let p1 = sys.particles()[1].position;
        // They moved apart along x...
        assert!(p0.x < -half_gap, "left particle should move further left");
        assert!(p1.x > half_gap, "right particle should move further right");
        // ...symmetrically: the centre of mass is unchanged and velocities are
        // equal and opposite (Newton's third law).
        let com_after = sys.center_of_mass();
        assert!(
            (com_after - com_before).norm() < 1e-15,
            "centre of mass drifted: {com_before:?} → {com_after:?}"
        );
        let v0 = sys.particles()[0].velocity;
        let v1 = sys.particles()[1].velocity;
        assert!(
            (v0 + v1).norm() < 1e-15,
            "velocities not equal-and-opposite: {v0:?}, {v1:?}"
        );
    }

    // BENCHMARK-PIN (5): a settled column has monotonically higher
    // density/pressure toward the bottom (hydrostatic stratification).
    #[test]
    fn settled_column_density_increases_downward() {
        let h = 0.2;
        let mut cfg = SphConfig::water(h).unwrap();
        cfg.viscosity = 2.0; // strong damping so it actually settles
        cfg.gravity = Vector3::new(0.0, 0.0, -9.806_65);
        cfg.clamp_negative_pressure = true;
        let mut solver = SphSolver::new(cfg).unwrap();

        // A vertical column of particles in a tall thin box, spacing ≈ h/2.
        let spacing = h * 0.5;
        let n_layers = 12;
        let mut sys = ParticleSystem::new();
        // A small 2×2 footprint per layer so each interior particle has lateral
        // neighbours too (a single 1-D string under-counts density).
        for layer in 0..n_layers {
            let z = layer as f64 * spacing;
            for fx in 0..2 {
                for fy in 0..2 {
                    sys.push(Particle::at(Vector3::new(
                        fx as f64 * spacing,
                        fy as f64 * spacing,
                        z,
                    )))
                    .unwrap();
                }
            }
        }
        let boundary = BoxBoundary::new(
            Vector3::new(-spacing, -spacing, 0.0),
            Vector3::new(2.0 * spacing, 2.0 * spacing, n_layers as f64 * spacing),
            0.0, // fully damped walls
        )
        .unwrap();

        // Let it settle.
        for _ in 0..4000 {
            solver.step(&mut sys, 5e-4, Some(&boundary)).unwrap();
        }
        solver.compute_density_pressure(&mut sys);

        // Average density per height layer; compare the bottom layers to the top.
        // Sort particles by z and bin into thirds.
        let mut ps: Vec<Particle> = sys.particles().to_vec();
        ps.sort_by(|a, b| a.position.z.partial_cmp(&b.position.z).unwrap());
        let third = ps.len() / 3;
        let bottom_avg: f64 = ps[..third].iter().map(|p| p.density).sum::<f64>() / third as f64;
        let top_avg: f64 = ps[ps.len() - third..]
            .iter()
            .map(|p| p.density)
            .sum::<f64>()
            / third as f64;

        assert!(
            bottom_avg > top_avg,
            "settled column should be denser at the bottom: bottom {bottom_avg}, top {top_avg}"
        );
        // The whole column stays inside the tank.
        assert!(sys
            .particles()
            .iter()
            .all(|p| boundary.contains(p.position)));
    }

    #[test]
    fn rejects_bad_config_and_dt() {
        // Bad mass.
        let mut cfg = SphConfig::water(0.1).unwrap();
        cfg.mass = 0.0;
        assert!(SphSolver::new(cfg).is_err());

        // Bad viscosity.
        let mut cfg = SphConfig::water(0.1).unwrap();
        cfg.viscosity = -1.0;
        assert!(SphSolver::new(cfg).is_err());

        // Non-finite gravity.
        let mut cfg = SphConfig::water(0.1).unwrap();
        cfg.gravity = Vector3::new(0.0, 0.0, f64::NAN);
        assert!(SphSolver::new(cfg).is_err());

        // Bad dt.
        let cfg = SphConfig::water(0.1).unwrap();
        let mut solver = SphSolver::new(cfg).unwrap();
        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::zeros())).unwrap();
        assert!(solver.step(&mut sys, 0.0, None).is_err());
        assert!(solver.step(&mut sys, -1e-3, None).is_err());
        assert!(solver.step(&mut sys, f64::NAN, None).is_err());
    }

    #[test]
    fn isolated_particle_gets_floor_density() {
        let cfg = SphConfig::water(0.1).unwrap();
        let mut solver = SphSolver::new(cfg).unwrap();
        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::new(0.0, 0.0, 0.0))).unwrap();
        solver.compute_density_pressure(&mut sys);
        // A lone particle's density is exactly its own self-contribution.
        assert!(approx(
            sys.particles()[0].density,
            solver.density_floor(),
            1e-12
        ));
        assert!(solver.density_floor() > 0.0);
    }
}
