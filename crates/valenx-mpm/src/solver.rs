//! The 2D MLS-MPM solver: particles, background grid, and one explicit step.
//!
//! Transfer scheme is APIC (Affine Particle-In-Cell, Jiang et al. 2015) with
//! the MLS-MPM force formulation (Hu et al. 2018), using quadratic B-spline
//! interpolation on a uniform Cartesian grid. The step is fully deterministic:
//! a fixed left-to-right particle iteration order and no randomness.
//!
//! Engineering framing: this models large-deformation continuum solids —
//! metal forming, soil/granular flow, structural collapse — where Lagrangian
//! mesh elements would invert and crash a classical FEM solver.

use crate::math::{Mat2, Vec2};
use crate::model::{ConstitutiveModel, ElasticParams};

/// Error conditions raised by the solver when given invalid state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MpmError {
    /// `step` was called with no particles.
    EmptyParticleSet,
    /// A particle held a non-finite position, velocity, `C`, or `F`.
    NonFiniteParticle {
        /// Index of the offending particle.
        index: usize,
    },
    /// A particle left the padded grid domain (its support touched the border).
    ParticleOutOfBounds {
        /// Index of the offending particle.
        index: usize,
    },
    /// The constitutive model returned no stress (singular/`NaN` `F`).
    StressEvaluationFailed {
        /// Index of the offending particle.
        index: usize,
    },
    /// A configuration value (`dx`, `dt`, grid size) was non-positive/non-finite.
    InvalidConfiguration,
}

impl core::fmt::Display for MpmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MpmError::EmptyParticleSet => write!(f, "MPM step requires at least one particle"),
            MpmError::NonFiniteParticle { index } => {
                write!(f, "particle {index} has a non-finite field")
            }
            MpmError::ParticleOutOfBounds { index } => {
                write!(f, "particle {index} left the grid domain")
            }
            MpmError::StressEvaluationFailed { index } => {
                write!(
                    f,
                    "stress evaluation failed for particle {index} (singular F)"
                )
            }
            MpmError::InvalidConfiguration => write!(f, "invalid MPM configuration"),
        }
    }
}

impl std::error::Error for MpmError {}

/// A single material point (particle) carrying the full continuum state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Particle {
    /// World-space position, in metres.
    pub position: Vec2,
    /// Velocity, in m/s.
    pub velocity: Vec2,
    /// Mass, in kilograms (must be positive).
    pub mass: f64,
    /// Initial (reference) volume `V₀`, in m² for this 2D solver (area per
    /// particle, must be positive). Scales the internal-stress force in the
    /// MLS-MPM transfer; getting it wrong rescales the effective stiffness.
    pub volume: f64,
    /// Affine velocity matrix `C` (the APIC local velocity gradient).
    pub affine_c: Mat2,
    /// Deformation gradient `F` (initially identity, evolves each step).
    pub deformation_f: Mat2,
}

impl Particle {
    /// Creates a particle at rest with identity deformation and zero affine `C`.
    ///
    /// The reference volume defaults to `mass` (i.e. unit reference density of
    /// `1 kg/m²`). For quantitatively meaningful stress, set the true per-
    /// particle area with [`Particle::with_volume`] or via
    /// [`crate::block_of_particles`], which derives it from the block geometry.
    #[must_use]
    pub fn new(position: Vec2, velocity: Vec2, mass: f64) -> Self {
        Self {
            position,
            velocity,
            mass,
            volume: mass,
            affine_c: Mat2::zero(),
            deformation_f: Mat2::identity(),
        }
    }

    /// Creates a particle with an explicit reference volume `V₀` (area, m²).
    #[must_use]
    pub fn with_volume(position: Vec2, velocity: Vec2, mass: f64, volume: f64) -> Self {
        Self {
            position,
            velocity,
            mass,
            volume,
            affine_c: Mat2::zero(),
            deformation_f: Mat2::identity(),
        }
    }

    /// Returns `true` iff every field is finite and `mass > 0`.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.position.is_finite()
            && self.velocity.is_finite()
            && self.mass.is_finite()
            && self.mass > 0.0
            && self.volume.is_finite()
            && self.volume > 0.0
            && self.affine_c.is_finite()
            && self.deformation_f.is_finite()
    }
}

/// One background-grid node: accumulated mass and momentum/velocity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct GridNode {
    /// Accumulated nodal mass (kg).
    mass: f64,
    /// During P2G: accumulated momentum (kg·m/s). After the grid update it is
    /// reinterpreted in place as velocity (m/s).
    vel: Vec2,
}

/// A uniform Cartesian background grid covering `[0, size·dx]²`.
#[derive(Clone, Debug)]
pub struct Grid {
    size: usize,
    dx: f64,
    inv_dx: f64,
    nodes: Vec<GridNode>,
}

impl Grid {
    /// Builds a `size × size` node grid with cell spacing `dx` (metres).
    ///
    /// # Panics
    /// Panics if `size == 0` or `dx` is non-finite/non-positive.
    #[must_use]
    pub fn new(size: usize, dx: f64) -> Self {
        assert!(size > 0, "grid size must be positive");
        assert!(
            dx.is_finite() && dx > 0.0,
            "grid spacing must be finite and positive"
        );
        Self {
            size,
            dx,
            inv_dx: 1.0 / dx,
            nodes: vec![GridNode::default(); size * size],
        }
    }

    /// Number of nodes along one axis.
    #[must_use]
    pub fn size(&self) -> usize {
        self.size
    }

    /// Cell spacing `dx` (metres).
    #[must_use]
    pub fn dx(&self) -> f64 {
        self.dx
    }

    /// Total mass currently accumulated on the grid (kg).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.nodes.iter().map(|n| n.mass).sum()
    }

    #[inline]
    fn clear(&mut self) {
        for n in &mut self.nodes {
            n.mass = 0.0;
            n.vel = Vec2::zero();
        }
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.size + i
    }
}

/// Quadratic B-spline weights and the base node for one particle.
///
/// Returns the integer base node `(bi, bj)` (the lower-left of the 3x3 stencil)
/// plus the per-axis weights `wx[0..3]`, `wy[0..3]`, and the fractional offset
/// of the particle from the base+1 node in grid units (used for `C` / APIC).
struct Stencil {
    base_i: i64,
    base_j: i64,
    wx: [f64; 3],
    wy: [f64; 3],
}

impl Stencil {
    /// Computes the quadratic-spline stencil for a particle at grid-space
    /// coordinate `gp = position / dx`.
    #[inline]
    fn new(gp: Vec2) -> Self {
        // Base node: round to nearest, then step back one (quadratic support).
        let base_i = (gp.x - 0.5).floor() as i64;
        let base_j = (gp.y - 0.5).floor() as i64;
        // fx is the offset of the particle from the base node (in [0.5, 1.5)).
        let fx = gp.x - base_i as f64;
        let fy = gp.y - base_j as f64;
        Self {
            base_i,
            base_j,
            wx: Self::weights(fx),
            wy: Self::weights(fy),
        }
    }

    /// Quadratic B-spline weights for the three stencil nodes given the
    /// fractional coordinate `f ∈ [0.5, 1.5)` of the particle from the base.
    #[inline]
    fn weights(f: f64) -> [f64; 3] {
        let a = 1.5 - f; // distance to node 0
        let b = f - 1.0; // signed distance to node 1
        let c = f - 0.5; // distance to node 2 (mirror)
        [0.5 * a * a, 0.75 - b * b, 0.5 * c * c]
    }
}

/// Configuration for an MLS-MPM simulation.
#[derive(Clone, Debug)]
pub struct MpmConfig {
    /// Background grid spacing `dx` (metres).
    pub dx: f64,
    /// Number of grid nodes per axis.
    pub grid_size: usize,
    /// Time step `dt` (seconds).
    pub dt: f64,
    /// Gravitational acceleration (m/s²), applied as a body force.
    pub gravity: Vec2,
    /// Lamé elastic parameters.
    pub material: ElasticParams,
    /// Elastic constitutive model.
    pub model: ConstitutiveModel,
    /// Floor height (metres). Nodes at or below this `y` get a no-penetration
    /// boundary condition; walls are applied at the domain edges.
    pub floor: f64,
}

impl MpmConfig {
    /// Validates the configuration, returning [`MpmError::InvalidConfiguration`]
    /// when any field is non-finite or out of range.
    fn validate(&self) -> Result<(), MpmError> {
        if self.dx.is_finite()
            && self.dx > 0.0
            && self.grid_size > 0
            && self.dt.is_finite()
            && self.dt > 0.0
            && self.gravity.is_finite()
            && self.material.is_finite()
            && self.floor.is_finite()
        {
            Ok(())
        } else {
            Err(MpmError::InvalidConfiguration)
        }
    }
}

/// The MLS-MPM solver: owns the particles, the grid, and the configuration.
#[derive(Clone, Debug)]
pub struct MpmSolver {
    /// The material points.
    pub particles: Vec<Particle>,
    /// The background grid (reused across steps).
    pub grid: Grid,
    /// Simulation configuration.
    pub config: MpmConfig,
}

impl MpmSolver {
    /// Constructs a solver from particles and configuration.
    ///
    /// # Panics
    /// Panics (via [`Grid::new`]) if the configured grid size/spacing is
    /// invalid. Per-particle and per-step validation is performed in [`Self::step`].
    #[must_use]
    pub fn new(particles: Vec<Particle>, config: MpmConfig) -> Self {
        let grid = Grid::new(config.grid_size, config.dx);
        Self {
            particles,
            grid,
            config,
        }
    }

    /// Total particle mass (kg).
    #[must_use]
    pub fn total_particle_mass(&self) -> f64 {
        self.particles.iter().map(|p| p.mass).sum()
    }

    /// Total particle momentum (kg·m/s) = Σ mᵢ vᵢ.
    #[must_use]
    pub fn total_momentum(&self) -> Vec2 {
        self.particles
            .iter()
            .fold(Vec2::zero(), |acc, p| acc + p.velocity * p.mass)
    }

    /// Centre of mass (metres) = (Σ mᵢ xᵢ) / (Σ mᵢ).
    #[must_use]
    pub fn center_of_mass(&self) -> Vec2 {
        let m = self.total_particle_mass();
        if m <= 0.0 {
            return Vec2::zero();
        }
        let weighted = self
            .particles
            .iter()
            .fold(Vec2::zero(), |acc, p| acc + p.position * p.mass);
        weighted * (1.0 / m)
    }

    /// Centre-of-mass velocity (m/s) = (Σ mᵢ vᵢ) / (Σ mᵢ).
    #[must_use]
    pub fn center_of_mass_velocity(&self) -> Vec2 {
        let m = self.total_particle_mass();
        if m <= 0.0 {
            return Vec2::zero();
        }
        self.total_momentum() * (1.0 / m)
    }

    /// Total kinetic energy (joules) = ½ Σ mᵢ ‖vᵢ‖².
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self
            .particles
            .iter()
            .map(|p| p.mass * p.velocity.length_squared())
            .sum::<f64>()
    }

    /// Advances the simulation by one explicit MLS-MPM step.
    ///
    /// Pipeline: validate → particle-to-grid (P2G) → grid update (body force +
    /// boundary conditions) → grid-to-particle (G2P, updating velocity, the
    /// affine `C`, position, and the deformation gradient `F`).
    ///
    /// # Errors
    /// Returns an [`MpmError`] if the particle set is empty, any particle is
    /// non-finite, a particle leaves the domain, a stress evaluation fails, or
    /// the configuration is invalid — the solver never silently produces
    /// `NaN`/`inf` state.
    pub fn step(&mut self) -> Result<(), MpmError> {
        self.config.validate()?;
        if self.particles.is_empty() {
            return Err(MpmError::EmptyParticleSet);
        }
        self.validate_particles()?;
        self.p2g()?;
        self.grid_update();
        self.g2p()?;
        Ok(())
    }

    fn validate_particles(&self) -> Result<(), MpmError> {
        let max = self.grid.size as i64;
        for (index, p) in self.particles.iter().enumerate() {
            if !p.is_valid() {
                return Err(MpmError::NonFiniteParticle { index });
            }
            // The quadratic stencil touches base..=base+2; keep a one-node
            // margin so writes stay inside the array.
            let gp = p.position * self.grid.inv_dx;
            let st = Stencil::new(gp);
            if st.base_i < 1 || st.base_j < 1 || st.base_i + 2 >= max || st.base_j + 2 >= max {
                return Err(MpmError::ParticleOutOfBounds { index });
            }
        }
        Ok(())
    }

    /// Particle-to-grid: scatter mass and (APIC-affine) momentum.
    fn p2g(&mut self) -> Result<(), MpmError> {
        self.grid.clear();
        let dx = self.grid.dx;
        let inv_dx = self.grid.inv_dx;
        // MLS-MPM stress→force scaling: −dt · V₀ · (4/dx²) · P·Fᵀ, where the
        // per-particle reference volume V₀ is applied below. `D⁻¹ = 4/dx²` is
        // the (inverse) inertia weighting for quadratic B-splines. Omitting V₀
        // overscales the internal force by ~1/V₀ and makes the scheme blow up.
        let d_inv = 4.0 * inv_dx * inv_dx;
        let dt = self.config.dt;

        for (index, p) in self.particles.iter().enumerate() {
            let gp = p.position * inv_dx;
            let st = Stencil::new(gp);

            // Affine momentum term: APIC uses C; MLS-MPM folds the internal
            // stress into the same affine channel as an extra contribution.
            let stress = self
                .config
                .model
                .pk1_stress(p.deformation_f, self.config.material)
                .ok_or(MpmError::StressEvaluationFailed { index })?;
            // affine = mass·C − dt·V₀·D⁻¹ · P·Fᵀ
            let p_ft = stress.mul_mat(p.deformation_f.transpose());
            let stress_coeff = -dt * p.volume * d_inv;
            let affine = p.affine_c.scale(p.mass).plus(p_ft.scale(stress_coeff));

            for gj in 0..3usize {
                for gi in 0..3usize {
                    let w = st.wx[gi] * st.wy[gj];
                    let ni = (st.base_i + gi as i64) as usize;
                    let nj = (st.base_j + gj as i64) as usize;
                    // dpos: node position minus particle position (world units).
                    let node_x = (st.base_i + gi as i64) as f64 * dx;
                    let node_y = (st.base_j + gj as i64) as f64 * dx;
                    let dpos = Vec2::new(node_x - p.position.x, node_y - p.position.y);

                    let idx = self.grid.idx(ni, nj);
                    let node = &mut self.grid.nodes[idx];
                    node.mass += w * p.mass;
                    // momentum += w·(mass·v + affine·dpos)
                    let momentum = p.velocity * p.mass + affine.mul_vec(dpos);
                    node.vel += momentum * w;
                }
            }
        }
        Ok(())
    }

    /// Grid update: convert momentum→velocity, apply body force and BCs.
    fn grid_update(&mut self) {
        let dt = self.config.dt;
        let gravity = self.config.gravity;
        let dx = self.grid.dx;
        let floor = self.config.floor;
        let size = self.grid.size;

        for j in 0..size {
            for i in 0..size {
                let idx = self.grid.idx(i, j);
                let node = &mut self.grid.nodes[idx];
                if node.mass <= 0.0 {
                    node.vel = Vec2::zero();
                    continue;
                }
                // momentum → velocity
                node.vel = node.vel * (1.0 / node.mass);
                // body force (gravity)
                node.vel += gravity * dt;

                // Boundary conditions: sticky walls at domain edges (2-node
                // pad) and a no-penetration floor.
                let y = j as f64 * dx;
                if i < 2 || i + 2 >= size {
                    node.vel.x = 0.0;
                }
                if j < 2 || j + 2 >= size {
                    node.vel.y = 0.0;
                }
                // Floor: block downward motion at/below the floor height.
                if y <= floor && node.vel.y < 0.0 {
                    node.vel.y = 0.0;
                }
            }
        }
    }

    /// Grid-to-particle: gather velocity + affine `C`, advect, update `F`.
    fn g2p(&mut self) -> Result<(), MpmError> {
        let dx = self.grid.dx;
        let inv_dx = self.grid.inv_dx;
        let dt = self.config.dt;
        // APIC reconstruction gain for quadratic splines: 4/dx².
        let c_gain = 4.0 * inv_dx * inv_dx;

        for (index, p) in self.particles.iter_mut().enumerate() {
            let gp = p.position * inv_dx;
            let st = Stencil::new(gp);

            let mut new_v = Vec2::zero();
            let mut new_c = Mat2::zero();

            for gj in 0..3usize {
                for gi in 0..3usize {
                    let w = st.wx[gi] * st.wy[gj];
                    let ni = (st.base_i + gi as i64) as usize;
                    let nj = (st.base_j + gj as i64) as usize;
                    let idx = self.grid.idx(ni, nj);
                    let g_vel = self.grid.nodes[idx].vel;

                    let node_x = (st.base_i + gi as i64) as f64 * dx;
                    let node_y = (st.base_j + gj as i64) as f64 * dx;
                    let dpos = Vec2::new(node_x - p.position.x, node_y - p.position.y);

                    new_v += g_vel * w;
                    // C += gain · w · (g_vel ⊗ dpos)
                    new_c = new_c.plus(g_vel.outer(dpos).scale(c_gain * w));
                }
            }

            p.velocity = new_v;
            p.affine_c = new_c;
            p.position += new_v * dt;

            // Deformation gradient update: F ← (I + dt·C) · F.
            let grad = Mat2::identity().plus(new_c.scale(dt));
            let new_f = grad.mul_mat(p.deformation_f);
            if !new_f.is_finite() {
                return Err(MpmError::NonFiniteParticle { index });
            }
            p.deformation_f = new_f;
        }
        Ok(())
    }
}
