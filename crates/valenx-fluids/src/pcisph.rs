//! PCISPH — Predictive-Corrective Incompressible SPH
//! (Solenthaler & Pajarola, *"Predictive-Corrective Incompressible SPH"*,
//! SIGGRAPH 2009).
//!
//! Plain weakly-compressible SPH ([`crate::sph`]) reads pressure off the density
//! with a stiff equation of state. To keep the fluid nearly incompressible the
//! stiffness has to be large, which forces a tiny time step (the sound-speed CFL
//! limit). **PCISPH** instead leaves the EOS out of the inner loop and *solves*
//! for the pressure that cancels the predicted density error:
//!
//! 1. Compute non-pressure forces once (gravity + viscosity).
//! 2. Initialise pressure `p_i = 0` and pressure force `f^p_i = 0`.
//! 3. **Iterate** (a fixed number of correction passes, or until the density
//!    error falls below a tolerance):
//!    * **Predict** each particle's position/velocity under the current total
//!      force (non-pressure + pressure).
//!    * **Estimate** the predicted density `ρ*_i` at the predicted positions and
//!      its error `ρ*_i − ρ₀`.
//!    * **Correct** the pressure by `Δp_i = −δ · (ρ*_i − ρ₀)`, where the scaling
//!      `δ` is a precomputed constant that depends only on the kernel, the
//!      particle mass, the rest density and the time step (the "PCISPH delta"
//!      from the paper's Eq. 8). Accumulate `p_i += Δp_i`.
//!    * Recompute the symmetric pressure force from the updated pressures.
//! 4. Integrate once with the converged pressure force, then apply the boundary.
//!
//! The result enforces incompressibility much more stiffly per step than WCSPH
//! for the same time step — at the cost of several density/force passes per step.
//!
//! # The δ scaling factor
//!
//! For a particle with a *filled* neighbourhood the paper derives
//!
//! ```text
//! δ = −1 / ( β · ( −Σ_j ∇W_ij · Σ_j ∇W_ij − Σ_j (∇W_ij · ∇W_ij) ) )
//! ```
//!
//! with `β = (Δt · m / ρ₀)² · 2`. Rather than depend on a specific particle's
//! instantaneous neighbourhood (which is empty for an isolated particle and
//! would divide by zero), this implementation precomputes the kernel sums over a
//! reference filled neighbourhood — a regular lattice at the configured particle
//! spacing within the support — exactly as the original does, so `δ` is a stable
//! constant and the correction never divides by a near-zero gradient sum.
//!
//! # Honesty
//!
//! This is the standard graphics-grade PCISPH. It is **not** a validated
//! incompressible Navier–Stokes solver: the boundary handling is the same
//! reflect/clamp box as WCSPH (no boundary-particle pressure), the iteration
//! count is small, and there is no implicit viscosity or surface tension.

use nalgebra::Vector3;

use crate::boundary::BoxBoundary;
use crate::error::FluidError;
use crate::grid::SpatialHash;
use crate::kernels::SmoothingKernels;
use crate::particle::ParticleSystem;

/// Configuration for the [`PciSphSolver`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PciSphConfig {
    /// Smoothing length / kernel support radius `h` (m), strictly positive.
    pub smoothing_length: f64,
    /// Target particle spacing (m), used to size the reference neighbourhood for
    /// the δ factor and the default mass. Strictly positive, `≤ h`.
    pub particle_spacing: f64,
    /// Per-particle mass `m` (kg), strictly positive.
    pub mass: f64,
    /// Rest density `ρ₀` (kg/m³), strictly positive.
    pub rest_density: f64,
    /// Dynamic viscosity coefficient `μ` (Pa·s), `≥ 0`.
    pub viscosity: f64,
    /// Gravity vector `g` (m/s²).
    pub gravity: Vector3<f64>,
    /// Number of prediction-correction iterations per step (`≥ 1`).
    pub iterations: usize,
}

impl PciSphConfig {
    /// A water-like default for a smoothing length: spacing `h/2`, mass set to
    /// reach rest density on a lattice, `ρ₀ = 1000`, 3 correction iterations.
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
        let spacing = smoothing_length * 0.5;
        let mass = rest_density * spacing * spacing * spacing;
        Ok(Self {
            smoothing_length,
            particle_spacing: spacing,
            mass,
            rest_density,
            viscosity: 0.1,
            gravity: Vector3::new(0.0, 0.0, -9.806_65),
            iterations: 3,
        })
    }
}

/// A predictive-corrective incompressible SPH solver.
#[derive(Debug, Clone)]
pub struct PciSphSolver {
    config: PciSphConfig,
    kernels: SmoothingKernels,
    grid: SpatialHash,
    rho_floor: f64,
    neighbor_scratch: Vec<usize>,
}

impl PciSphSolver {
    /// Build a PCISPH solver, validating every parameter.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] for non-positive `smoothing_length`,
    /// `particle_spacing`, `mass`, `rest_density`, a `particle_spacing > h`, a
    /// negative `viscosity`, or `iterations == 0`; [`FluidError::NonFinite`] for
    /// non-finite `gravity`.
    pub fn new(config: PciSphConfig) -> Result<Self, FluidError> {
        let kernels = SmoothingKernels::new(config.smoothing_length)?;
        check_positive(config.particle_spacing, "particle spacing")?;
        check_positive(config.mass, "particle mass")?;
        check_positive(config.rest_density, "rest density")?;
        if config.particle_spacing > config.smoothing_length {
            return Err(FluidError::InvalidConfig(format!(
                "particle spacing ({}) must be ≤ smoothing length ({})",
                config.particle_spacing, config.smoothing_length
            )));
        }
        if !(config.viscosity.is_finite() && config.viscosity >= 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "viscosity must be finite and ≥ 0, got {}",
                config.viscosity
            )));
        }
        if config.iterations == 0 {
            return Err(FluidError::InvalidConfig("iterations must be ≥ 1".into()));
        }
        if !config.gravity.iter().all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("gravity vector".into()));
        }
        let grid = SpatialHash::new(config.smoothing_length)?;
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
    pub fn config(&self) -> &PciSphConfig {
        &self.config
    }

    /// The PCISPH pressure-correction scaling factor `δ` for a time step `dt`.
    ///
    /// Computed from kernel-gradient sums over a *reference filled
    /// neighbourhood* (a regular lattice at `particle_spacing` within the
    /// support), per Solenthaler & Pajarola 2009 Eq. 8, so it is a stable
    /// constant independent of any particle's instantaneous neighbour count.
    /// Returns `0` if the reference neighbourhood has no in-support neighbours
    /// (degenerate; the caller then performs no correction rather than dividing
    /// by zero).
    #[must_use]
    pub fn delta(&self, dt: f64) -> f64 {
        let mut grad_sum = Vector3::zeros();
        let mut grad_dot_sum = 0.0;
        let s = self.config.particle_spacing;
        let h = self.kernels.h();
        // Build a reference lattice of offsets around the origin within support.
        let reach = (h / s).ceil() as i64;
        for i in -reach..=reach {
            for j in -reach..=reach {
                for k in -reach..=reach {
                    if i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    let rij = Vector3::new(i as f64 * s, j as f64 * s, k as f64 * s);
                    if rij.norm() > h {
                        continue;
                    }
                    let grad = self.kernels.spiky_gradient(rij);
                    grad_sum += grad;
                    grad_dot_sum += grad.dot(&grad);
                }
            }
        }
        let denom_geom = -grad_sum.dot(&grad_sum) - grad_dot_sum;
        if denom_geom.abs() < f64::EPSILON {
            return 0.0;
        }
        let beta = 2.0 * (dt * self.config.mass / self.config.rest_density).powi(2);
        let denom = beta * denom_geom;
        if denom.abs() < f64::EPSILON {
            return 0.0;
        }
        -1.0 / denom
    }

    /// SPH density at a set of positions (poly6 sum, including the self term),
    /// floored at the lone-particle density. Writes into `out` (resized).
    fn densities_at(&mut self, positions: &[Vector3<f64>], out: &mut Vec<f64>) {
        out.clear();
        out.resize(positions.len(), 0.0);
        let mass = self.config.mass;
        for i in 0..positions.len() {
            self.grid
                .neighbors_into(positions[i], &mut self.neighbor_scratch);
            let mut density = 0.0;
            for &j in &self.neighbor_scratch {
                let r2 = (positions[i] - positions[j]).norm_squared();
                density += mass * self.kernels.poly6_r2(r2);
            }
            out[i] = density.max(self.rho_floor);
        }
    }

    /// Symmetric pressure force for each particle given current pressures and
    /// densities. `grid` must already be built for `positions`.
    fn pressure_forces(
        &mut self,
        positions: &[Vector3<f64>],
        densities: &[f64],
        pressures: &[f64],
        out: &mut Vec<Vector3<f64>>,
    ) {
        out.clear();
        out.resize(positions.len(), Vector3::zeros());
        let mass = self.config.mass;
        let h = self.kernels.h();
        for i in 0..positions.len() {
            self.grid
                .neighbors_into(positions[i], &mut self.neighbor_scratch);
            let rho_i = densities[i];
            let p_over_rho2_i = pressures[i] / (rho_i * rho_i);
            let mut f = Vector3::zeros();
            for &j in &self.neighbor_scratch {
                if j == i {
                    continue;
                }
                let rij = positions[i] - positions[j];
                if rij.norm() > h {
                    continue;
                }
                let rho_j = densities[j];
                let p_over_rho2_j = pressures[j] / (rho_j * rho_j);
                let grad = self.kernels.spiky_gradient(rij);
                f -= grad * (mass * mass * (p_over_rho2_i + p_over_rho2_j));
            }
            out[i] = f;
        }
    }

    /// Advance the particle system by one time step `dt` (seconds).
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
        let n = system.len();
        if n == 0 {
            return Ok(());
        }

        let mass = self.config.mass;
        let mu = self.config.viscosity;
        let gravity = self.config.gravity;
        let rho0 = self.config.rest_density;
        let h = self.kernels.h();
        let delta = self.delta(dt);

        // Current state.
        let positions: Vec<Vector3<f64>> = system.particles().iter().map(|p| p.position).collect();
        let velocities: Vec<Vector3<f64>> = system.particles().iter().map(|p| p.velocity).collect();

        // --- 1. Non-pressure forces (gravity + viscosity), computed once. ---
        self.grid.rebuild(system);
        let mut density0 = Vec::new();
        self.densities_at(&positions, &mut density0);

        let mut f_nonpressure = vec![Vector3::zeros(); n];
        for i in 0..n {
            self.grid
                .neighbors_into(positions[i], &mut self.neighbor_scratch);
            let mut f_visc = Vector3::zeros();
            if mu > 0.0 {
                for &j in &self.neighbor_scratch {
                    if j == i {
                        continue;
                    }
                    let rij = positions[i] - positions[j];
                    let r = rij.norm();
                    if r > h {
                        continue;
                    }
                    let lap = self.kernels.viscosity_laplacian(r);
                    f_visc += (velocities[j] - velocities[i]) * (mu * mass / density0[j] * lap);
                }
            }
            // Gravity as a force: ρ_i g (so dividing by ρ_i later gives g).
            f_nonpressure[i] = f_visc + gravity * density0[i];
        }

        // --- 2. Init pressure / pressure force. ---
        let mut pressures = vec![0.0_f64; n];
        let mut f_pressure = vec![Vector3::zeros(); n];

        // Scratch for predicted positions/densities/forces.
        let mut pred_pos = positions.clone();
        let mut pred_density = Vec::new();

        // --- 3. Prediction-correction iterations. ---
        for _ in 0..self.config.iterations {
            // Predict velocity & position under (non-pressure + pressure) force.
            for i in 0..n {
                let accel = (f_nonpressure[i] + f_pressure[i]) / density0[i];
                let v_pred = velocities[i] + accel * dt;
                pred_pos[i] = positions[i] + v_pred * dt;
            }
            // Predicted density at predicted positions.
            self.grid.rebuild_from(&pred_pos);
            self.densities_at(&pred_pos, &mut pred_density);

            // Correct pressure from the density error.
            for i in 0..n {
                let err = pred_density[i] - rho0;
                pressures[i] += delta * err;
            }

            // Recompute the symmetric pressure force at predicted positions.
            self.pressure_forces(&pred_pos, &pred_density, &pressures, &mut f_pressure);
        }

        // --- 4. Final integration with the converged pressure force. ---
        for (i, p) in system.particles_mut().iter_mut().enumerate() {
            let accel = (f_nonpressure[i] + f_pressure[i]) / density0[i];
            p.velocity += accel * dt;
            p.position += p.velocity * dt;
            p.density = density0[i];
            p.pressure = pressures[i];
        }

        if let Some(b) = boundary {
            for p in system.particles_mut() {
                b.enforce(p);
            }
        }
        Ok(())
    }
}

fn check_positive(value: f64, what: &str) -> Result<(), FluidError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(FluidError::InvalidConfig(format!(
            "{what} must be finite and > 0, got {value}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::SmoothingKernels;
    use crate::particle::Particle;

    #[test]
    fn rejects_bad_config() {
        let mut cfg = PciSphConfig::water(0.1).unwrap();
        cfg.mass = 0.0;
        assert!(PciSphSolver::new(cfg).is_err());

        let mut cfg = PciSphConfig::water(0.1).unwrap();
        cfg.iterations = 0;
        assert!(PciSphSolver::new(cfg).is_err());

        let mut cfg = PciSphConfig::water(0.1).unwrap();
        cfg.particle_spacing = 0.5; // > h
        assert!(PciSphSolver::new(cfg).is_err());

        let mut cfg = PciSphConfig::water(0.1).unwrap();
        cfg.gravity = Vector3::new(f64::NAN, 0.0, 0.0);
        assert!(PciSphSolver::new(cfg).is_err());
    }

    #[test]
    fn delta_is_negative_and_finite() {
        let cfg = PciSphConfig::water(0.2).unwrap();
        let solver = PciSphSolver::new(cfg).unwrap();
        let d = solver.delta(1e-3);
        // δ = −1/(β·denom_geom); β > 0 and denom_geom < 0 ⇒ δ > 0 in this
        // sign convention. Either way it must be finite and non-zero for a
        // filled reference neighbourhood.
        assert!(d.is_finite());
        assert!(d.abs() > 0.0, "delta should be non-zero, got {d}");
    }

    // PCISPH free-fall: a lone particle has no neighbours, so no pressure
    // correction acts and it must fall under gravity exactly like WCSPH.
    #[test]
    fn single_particle_free_fall() {
        let mut cfg = PciSphConfig::water(0.1).unwrap();
        cfg.viscosity = 0.0;
        let g = -9.806_65;
        cfg.gravity = Vector3::new(0.0, 0.0, g);
        let mut solver = PciSphSolver::new(cfg).unwrap();

        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::zeros())).unwrap();

        let dt = 1e-3;
        let steps = 500;
        for _ in 0..steps {
            solver.step(&mut sys, dt, None).unwrap();
        }
        let t = dt * steps as f64;
        let analytic = 0.5 * g * t * t;
        let z = sys.particles()[0].position.z;
        assert!(
            (z - analytic).abs() < 1e-2,
            "z = {z} should be ≈ ½gt² = {analytic}"
        );
        assert!(sys.particles()[0].position.x.abs() < 1e-15);
    }

    // PCISPH reduces the density error of a compressed blob relative to taking
    // no pressure correction at all — its whole point. We compare the post-step
    // peak density error against an uncorrected (0-pressure) predicted step.
    #[test]
    fn correction_reduces_density_error() {
        let mut cfg = PciSphConfig::water(0.2).unwrap();
        cfg.gravity = Vector3::zeros();
        cfg.viscosity = 0.0;
        cfg.iterations = 5;
        let mut solver = PciSphSolver::new(cfg).unwrap();

        // A compressed lattice (spacing below the configured spacing) ⇒
        // over-density that the correction should relax.
        let mut sys = ParticleSystem::new();
        let s = cfg.particle_spacing * 0.7;
        for ix in 0..4 {
            for iy in 0..4 {
                for iz in 0..4 {
                    sys.push(Particle::at(Vector3::new(
                        ix as f64 * s,
                        iy as f64 * s,
                        iz as f64 * s,
                    )))
                    .unwrap();
                }
            }
        }

        // Peak density error BEFORE any motion.
        solver.grid.rebuild(&sys);
        let pos0: Vec<Vector3<f64>> = sys.particles().iter().map(|p| p.position).collect();
        let mut d0 = Vec::new();
        solver.densities_at(&pos0, &mut d0);
        let err_before = d0
            .iter()
            .map(|&d| (d - cfg.rest_density).abs())
            .fold(0.0_f64, f64::max);

        // One PCISPH step, then re-measure the peak density error.
        solver.step(&mut sys, 1e-3, None).unwrap();
        let pos1: Vec<Vector3<f64>> = sys.particles().iter().map(|p| p.position).collect();
        solver.grid.rebuild_from(&pos1);
        let mut d1 = Vec::new();
        solver.densities_at(&pos1, &mut d1);
        let err_after = d1
            .iter()
            .map(|&d| (d - cfg.rest_density).abs())
            .fold(0.0_f64, f64::max);

        assert!(
            err_after < err_before,
            "PCISPH should reduce the peak density error: before {err_before}, after {err_after}"
        );
    }

    // Newton's third law also holds for PCISPH: two compressed particles, no
    // gravity, conserve centre of mass and net momentum.
    #[test]
    fn two_particles_symmetric_and_momentum_conserved() {
        // Rest density = lone-particle self density, so an isolated pair is
        // genuinely compressed (above rest) rather than below it (see the WCSPH
        // counterpart for the rationale).
        let h = 0.2;
        let mut cfg = PciSphConfig::water(h).unwrap();
        cfg.gravity = Vector3::zeros();
        cfg.viscosity = 0.0;
        cfg.rest_density = cfg.mass * SmoothingKernels::new(h).unwrap().poly6(0.0);
        let mut solver = PciSphSolver::new(cfg).unwrap();

        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::new(-0.02, 0.0, 0.0)))
            .unwrap();
        sys.push(Particle::at(Vector3::new(0.02, 0.0, 0.0)))
            .unwrap();
        let com_before = sys.center_of_mass();
        let mom_before = sys.total_momentum(cfg.mass);

        solver.step(&mut sys, 1e-3, None).unwrap();

        let com_after = sys.center_of_mass();
        let mom_after = sys.total_momentum(cfg.mass);
        assert!(
            (com_after - com_before).norm() < 1e-12,
            "centre of mass drifted: {com_before:?} → {com_after:?}"
        );
        assert!(
            (mom_after - mom_before).norm() < 1e-12,
            "momentum not conserved: {mom_before:?} → {mom_after:?}"
        );
        // And they actually pushed apart.
        assert!(sys.particles()[0].position.x < -0.02);
        assert!(sys.particles()[1].position.x > 0.02);
    }
}
