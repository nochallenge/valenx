//! `valenx-mpm` — Material Point Method (MLS-MPM): large-deformation
//! continuum mechanics for engineering materials, structural failure,
//! and geotechnics.
//!
//! # What this crate is for
//!
//! Mesh-based finite-element solvers tangle and crash when elements invert —
//! exactly the regime of metal forming, soil/granular collapse, foam crush,
//! and ductile structural failure. The Material Point Method sidesteps element
//! inversion by carrying the material state on Lagrangian *particles* and
//! solving momentum on a disposable Eulerian *grid* each step. This is the V1
//! 2D solver: **MLS-MPM** (Hu et al. 2018) with the **APIC** transfer (Jiang
//! et al. 2015), quadratic B-spline interpolation, and a choice of neo-Hookean
//! or fixed-corotated elasticity.
//!
//! Scope is strictly engineering continuum mechanics — structural and material
//! response, geotechnics, and forming processes. It is not a model of
//! weapons-penetration or lethality.
//!
//! # The step
//!
//! Each [`MpmSolver::step`] runs the standard MLS-MPM pipeline:
//! 1. **P2G** — scatter particle mass and (affine + internal-stress) momentum
//!    to the grid via quadratic B-spline weights.
//! 2. **Grid update** — convert momentum to velocity, apply gravity, and
//!    enforce wall/floor boundary conditions.
//! 3. **G2P** — gather grid velocity back to particles, reconstruct the affine
//!    matrix `C`, advect positions, and update each deformation gradient `F`.
//!
//! The step is deterministic (fixed iteration order, no randomness) and fails
//! loud on empty/non-finite/out-of-bounds input rather than emitting `NaN`.
//!
//! # Example
//!
//! ```
//! use valenx_mpm::{
//!     ConstitutiveModel, ElasticParams, MpmConfig, MpmSolver, Particle, Vec2,
//! };
//!
//! // One particle in free fall under gravity (no other forces).
//! let particle = Particle::new(Vec2::new(0.5, 0.8), Vec2::zero(), 1.0);
//! let config = MpmConfig {
//!     dx: 1.0 / 64.0,
//!     grid_size: 64,
//!     dt: 1.0e-4,
//!     gravity: Vec2::new(0.0, -9.81),
//!     material: ElasticParams::from_youngs(1.0e5, 0.3),
//!     model: ConstitutiveModel::FixedCorotated,
//!     floor: 0.05,
//! };
//! let mut solver = MpmSolver::new(vec![particle], config);
//! solver.step().expect("valid step");
//! assert!(solver.center_of_mass().is_finite());
//! ```

#![forbid(unsafe_code)]

pub mod math;
pub mod model;
pub mod solver;

pub use math::{Mat2, Vec2};
pub use model::{ConstitutiveModel, ElasticParams};
pub use solver::{Grid, MpmConfig, MpmError, MpmSolver, Particle};

/// Builds a rectangular block of particles laid out on a regular lattice.
///
/// Convenience helper for tests and demos: fills `[lo, hi]` (world metres) with
/// `nx × ny` particles, each given `velocity`, and assigns mass so the block's
/// total mass equals `total_mass`. Deterministic row-major ordering.
///
/// # Panics
/// Panics if `nx == 0` or `ny == 0`, if `total_mass <= 0`, or if any bound is
/// non-finite or `hi <= lo` on either axis.
#[must_use]
pub fn block_of_particles(
    lo: Vec2,
    hi: Vec2,
    nx: usize,
    ny: usize,
    velocity: Vec2,
    total_mass: f64,
) -> Vec<Particle> {
    assert!(
        nx > 0 && ny > 0,
        "block needs at least one particle per axis"
    );
    assert!(
        total_mass.is_finite() && total_mass > 0.0,
        "total mass must be positive"
    );
    assert!(
        lo.is_finite() && hi.is_finite(),
        "block bounds must be finite"
    );
    assert!(
        hi.x > lo.x && hi.y > lo.y,
        "block hi must exceed lo on both axes"
    );

    let count = nx * ny;
    let mass_each = total_mass / count as f64;
    // Reference volume (area, m²) each particle represents = block area / count.
    let volume_each = ((hi.x - lo.x) * (hi.y - lo.y)) / count as f64;
    let mut out = Vec::with_capacity(count);
    for j in 0..ny {
        for i in 0..nx {
            // Cell-centred placement within the block extent.
            let fx = (i as f64 + 0.5) / nx as f64;
            let fy = (j as f64 + 0.5) / ny as f64;
            let pos = Vec2::new(lo.x + fx * (hi.x - lo.x), lo.y + fy * (hi.y - lo.y));
            out.push(Particle::with_volume(pos, velocity, mass_each, volume_each));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fall_config(gravity: f64) -> MpmConfig {
        MpmConfig {
            dx: 1.0 / 64.0,
            grid_size: 64,
            dt: 1.0e-4,
            gravity: Vec2::new(0.0, gravity),
            material: ElasticParams::from_youngs(1.0e5, 0.3),
            model: ConstitutiveModel::FixedCorotated,
            floor: 0.0,
        }
    }

    /// (a) FREE-FALL: with only gravity acting, the centre of mass must follow
    /// the analytic ballistic drop `Δy = ½ g t²` (here g < 0). MPM integrates
    /// velocity then position, so we compare against the matching discrete sum
    /// to isolate physics correctness from integrator order.
    #[test]
    fn free_fall_matches_ballistic_drop() {
        let g = -9.81;
        let cfg = fall_config(g);
        let dt = cfg.dt;
        // A small block away from any boundary so no BC interferes.
        let particles = block_of_particles(
            Vec2::new(0.40, 0.60),
            Vec2::new(0.50, 0.70),
            4,
            4,
            Vec2::zero(),
            1.0,
        );
        let y0 = {
            let s = MpmSolver::new(particles.clone(), cfg.clone());
            s.center_of_mass().y
        };
        let mut solver = MpmSolver::new(particles, cfg);

        let steps = 200;
        for _ in 0..steps {
            solver.step().expect("free-fall step");
        }
        let t = steps as f64 * dt;
        let dy = solver.center_of_mass().y - y0;

        // Closed-form continuous drop.
        let analytic = 0.5 * g * t * t;
        // Symplectic-Euler discrete drop: Σ_{k=1..N} g·dt² · k.
        let n = steps as f64;
        let discrete = g * dt * dt * (n * (n + 1.0) / 2.0);

        // Must match the discrete integral very tightly, and the closed form
        // to within one integration step's worth of drift.
        assert!(
            (dy - discrete).abs() < 1e-9,
            "drop {dy} vs discrete {discrete}"
        );
        assert!(
            (dy - analytic).abs() < (g.abs() * dt * t),
            "drop {dy} vs analytic {analytic}"
        );
        // Velocity must equal g·t (within a step).
        let vy = solver.center_of_mass_velocity().y;
        assert!(
            (vy - g * t).abs() < g.abs() * dt * 1.5,
            "vy {vy} vs {}",
            g * t
        );
    }

    /// (b) MASS CONSERVATION: P2G scatters with partition-of-unity weights, so
    /// the total grid mass after P2G must equal the total particle mass.
    #[test]
    fn p2g_conserves_mass() {
        let cfg = fall_config(-9.81);
        let particles = block_of_particles(
            Vec2::new(0.30, 0.30),
            Vec2::new(0.70, 0.70),
            8,
            8,
            Vec2::new(0.5, -0.5),
            3.7,
        );
        let mut solver = MpmSolver::new(particles, cfg);
        let particle_mass = solver.total_particle_mass();
        solver.step().expect("step");
        let grid_mass = solver.grid.total_mass();
        assert!(
            (grid_mass - particle_mass).abs() < 1e-9 * particle_mass.max(1.0),
            "grid mass {grid_mass} != particle mass {particle_mass}"
        );
    }

    /// (c) MOMENTUM CONSERVATION: in a gravity-free, internal-stress-free step
    /// (F = I ⇒ zero stress) the round trip P2G→grid→G2P must preserve total
    /// linear momentum to round-off.
    #[test]
    fn momentum_conserved_without_gravity() {
        let mut cfg = fall_config(0.0); // no gravity
        cfg.floor = -1.0; // floor below domain so it never triggers
                          // Keep the block well inside the domain and away from walls.
        let particles = block_of_particles(
            Vec2::new(0.40, 0.40),
            Vec2::new(0.60, 0.60),
            6,
            6,
            Vec2::new(0.8, 0.3),
            2.0,
        );
        let mut solver = MpmSolver::new(particles, cfg);
        let p0 = solver.total_momentum();
        // Single step keeps every particle in the interior (no wall BC hit).
        solver.step().expect("step");
        let p1 = solver.total_momentum();
        assert!(
            (p1.x - p0.x).abs() < 1e-9 && (p1.y - p0.y).abs() < 1e-9,
            "momentum {p0:?} -> {p1:?}"
        );
    }

    /// (d) ELASTIC BLOCK ON A FLOOR: an elastic block dropped onto the floor
    /// must settle/bounce with BOUNDED energy over many steps — no explosion,
    /// no `NaN`. We assert the run stays finite and the kinetic energy never
    /// exceeds a generous multiple of the initial gravitational budget.
    #[test]
    fn elastic_block_settles_without_exploding() {
        let cfg = MpmConfig {
            dx: 1.0 / 64.0,
            grid_size: 64,
            dt: 2.0e-4,
            gravity: Vec2::new(0.0, -9.81),
            material: ElasticParams::from_youngs(5.0e4, 0.2),
            model: ConstitutiveModel::FixedCorotated,
            floor: 0.10,
        };
        let particles = block_of_particles(
            Vec2::new(0.40, 0.30),
            Vec2::new(0.60, 0.50),
            10,
            10,
            Vec2::zero(),
            1.0,
        );
        let mut solver = MpmSolver::new(particles, cfg);

        // Energy ceiling: total mass · |g| · initial drop height, times a
        // comfortable factor. The block can fall at most ~0.5 m in the domain.
        let m = solver.total_particle_mass();
        let ceiling = m * 9.81 * 0.5 * 20.0;

        for k in 0..600 {
            solver
                .step()
                .unwrap_or_else(|e| panic!("step {k} failed: {e}"));
            let ke = solver.kinetic_energy();
            assert!(ke.is_finite(), "KE non-finite at step {k}");
            assert!(
                ke < ceiling,
                "energy blew up at step {k}: {ke} >= {ceiling}"
            );
            for p in &solver.particles {
                assert!(p.is_valid(), "invalid particle at step {k}");
            }
        }
        // The block should have moved downward toward/onto the floor.
        assert!(solver.center_of_mass().y < 0.40, "block did not fall");
    }

    /// (e) BAD INPUT fails loud — empty set, non-finite, out-of-bounds, bad cfg.
    #[test]
    fn bad_input_fails_loud() {
        let cfg = fall_config(-9.81);

        // Empty particle set.
        let mut empty = MpmSolver::new(Vec::new(), cfg.clone());
        assert_eq!(empty.step(), Err(MpmError::EmptyParticleSet));

        // Non-finite velocity.
        let mut bad = Particle::new(Vec2::new(0.5, 0.5), Vec2::zero(), 1.0);
        bad.velocity = Vec2::new(f64::NAN, 0.0);
        let mut s = MpmSolver::new(vec![bad], cfg.clone());
        assert_eq!(s.step(), Err(MpmError::NonFiniteParticle { index: 0 }));

        // Out-of-bounds position (outside the padded domain).
        let oob = Particle::new(Vec2::new(100.0, 0.5), Vec2::zero(), 1.0);
        let mut s2 = MpmSolver::new(vec![oob], cfg.clone());
        assert_eq!(s2.step(), Err(MpmError::ParticleOutOfBounds { index: 0 }));

        // Invalid configuration (dt = 0).
        let mut bad_cfg = cfg;
        bad_cfg.dt = 0.0;
        let good = Particle::new(Vec2::new(0.5, 0.5), Vec2::zero(), 1.0);
        let mut s3 = MpmSolver::new(vec![good], bad_cfg);
        assert_eq!(s3.step(), Err(MpmError::InvalidConfiguration));
    }

    /// Determinism: two identical solvers produce bit-identical state.
    #[test]
    fn step_is_deterministic() {
        let cfg = fall_config(-9.81);
        let particles = block_of_particles(
            Vec2::new(0.40, 0.40),
            Vec2::new(0.60, 0.60),
            5,
            5,
            Vec2::new(0.1, 0.0),
            1.0,
        );
        let mut a = MpmSolver::new(particles.clone(), cfg.clone());
        let mut b = MpmSolver::new(particles, cfg);
        for _ in 0..50 {
            a.step().unwrap();
            b.step().unwrap();
        }
        for (pa, pb) in a.particles.iter().zip(b.particles.iter()) {
            assert_eq!(pa, pb);
        }
    }

    /// The neo-Hookean model also runs stably end-to-end.
    #[test]
    fn neo_hookean_runs_stable() {
        let cfg = MpmConfig {
            dx: 1.0 / 64.0,
            grid_size: 64,
            dt: 1.0e-4,
            gravity: Vec2::new(0.0, -9.81),
            material: ElasticParams::from_youngs(1.0e5, 0.3),
            model: ConstitutiveModel::NeoHookean,
            floor: 0.10,
        };
        let particles = block_of_particles(
            Vec2::new(0.42, 0.40),
            Vec2::new(0.58, 0.56),
            8,
            8,
            Vec2::zero(),
            1.0,
        );
        let mut solver = MpmSolver::new(particles, cfg);
        for k in 0..300 {
            solver.step().unwrap_or_else(|e| panic!("nh step {k}: {e}"));
            assert!(solver.kinetic_energy().is_finite());
        }
    }
}
