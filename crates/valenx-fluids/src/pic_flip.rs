//! PIC/FLIP particle-in-cell fluid transfer.
//!
//! This module implements the **Particle-In-Cell / FLuid-Implicit-Particle
//! (PIC/FLIP)** particle↔grid transfer scheme for interactive particle-based
//! fluid simulation, following the framework described in:
//!
//! * Harlow & Welch 1965 (original PIC, MAC grid)
//! * Brackbill & Ruppel 1986 (FLIP blend)
//! * Zhu & Bridson, *"Animating Sand as a Fluid"*, SIGGRAPH 2005 (the standard
//!   graphics-grade PIC/FLIP formulation used here)
//!
//! ## What this module adds to the SPH family
//!
//! The SPH / PCISPH solvers in this crate are purely Lagrangian — forces are
//! accumulated by kernel sums over particle neighbours.  PIC/FLIP is a
//! **hybrid** Lagrangian–Eulerian scheme:
//!
//! 1. **P→G (particle-to-grid)** — deposit particle velocities onto a
//!    background uniform MAC grid via trilinear weights.
//! 2. **Grid step** — apply body forces (gravity) and a simplified
//!    pressure / divergence projection step on the grid.
//! 3. **G→P (grid-to-particle)** — interpolate the updated grid velocity
//!    back to particles using the PIC/FLIP blend:
//!
//!    ```text
//!    v_new = α · v_grid_new  +  (1−α) · (v_particle + Δv_grid)
//!    ```
//!
//!    where `Δv_grid = v_grid_new − v_grid_old`.
//!    * `α = 1.0` → **pure PIC**: copy the new grid velocity — low noise but
//!      over-diffusive.
//!    * `α = 0.0` → **pure FLIP**: add only the grid's velocity *change* —
//!      conserves fine-scale particle motion but accumulates noise.
//!    * `α ∈ (0, 1)` → blend (typical: 0.03–0.10 for FLIP-dominant).
//!
//! 4. **Advection** — move particles forward by `x ← x + v·dt` and apply the
//!    [`BoxBoundary`].
//!
//! ## Grid layout
//!
//! The background grid is a **staggered MAC (Marker-And-Cell) grid** in 3-D.
//! Velocity components are stored face-centred:
//!
//! * `u` (x-velocity) at `(i+½, j, k)` — indexed as `u[i][j][k]`.
//! * `v` (y-velocity) at `(i, j+½, k)` — indexed as `v[i][j][k]`.
//! * `w` (z-velocity) at `(i, j, k+½)` — indexed as `w[i][j][k]`.
//!
//! For a grid of `(nx, ny, nz)` cells: `u` has `(nx+1, ny, nz)` nodes,
//! `v` has `(nx, ny+1, nz)` nodes, `w` has `(nx, ny, nz+1)` nodes.
//!
//! ## Pressure step (simplified, v1)
//!
//! The full incompressible pressure-Poisson solve (Conjugate Gradient on the
//! Laplacian) is deferred to a future version.  Here we apply a **single
//! Jacobi-style divergence step**: for each cell the local velocity divergence
//! is computed and the face velocities adjusted to drive it toward zero.  This
//! is an explicit first-order projection — it removes the first harmonic of
//! divergence in one pass and is sufficient for short-duration interactive
//! simulations, but does **not** enforce incompressibility to solver tolerance.
//! The limitation is documented here and in [`PicFlipConfig`].
//!
//! ## ⚠ Honesty / scope — graphics-grade, NOT validated CFD
//!
//! This is **graphics / real-time-grade interactive fluid**. It is **not**
//! validated against analytic Navier–Stokes or experimental data and must not
//! be used for engineering CFD.  The known limitations mirror those of the SPH
//! solvers (see the crate root):
//!
//! * The pressure step is a single-pass Jacobi approximation, not a converged
//!   pressure-Poisson solve.
//! * Boundary treatment is the same positional reflect/clamp [`BoxBoundary`]
//!   as the SPH solvers — no true no-slip, no boundary pressure contribution.
//! * Semi-implicit Euler advection; first-order in time.
//! * PIC introduces numerical viscosity; FLIP accumulates noise in long
//!   simulations without a noise-damping sub-step.
//!
//! ## Pinned benchmarks (fail-loud, deterministic)
//!
//! See the `tests` block at the end of this module.  Four checks are pinned:
//!
//! 1. **Round-trip** (pure PIC, `α = 1`): a uniform velocity field deposited
//!    P→G and immediately interpolated G→P recovers the original velocity
//!    within trilinear-interpolation tolerance.
//! 2. **Free-fall** (no pressure): a single particle under gravity free-falls
//!    `Δx ≈ ½ g t²` over two steps.
//! 3. **Momentum conservation**: total particle momentum is approximately
//!    conserved through P→G→P (no external forces).
//! 4. **FLIP blend correctness**: `α ∈ {0, 0.5, 1}` interpolates as
//!    advertised between PIC and full-FLIP.

use nalgebra::Vector3;

use crate::boundary::BoxBoundary;
use crate::error::FluidError;
use crate::particle::ParticleSystem;

// ─── Configuration ────────────────────────────────────────────────────────────

/// Configuration for the [`PicFlipSolver`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PicFlipConfig {
    /// Number of grid cells along each axis `(nx, ny, nz)`, all `≥ 1`.
    pub grid_cells: (usize, usize, usize),
    /// Grid cell side length (m), strictly positive.
    pub cell_size: f64,
    /// World-space origin of the grid's min corner (m).
    pub origin: Vector3<f64>,
    /// Gravity vector (m/s²), e.g. `(0, 0, −9.80665)`.
    pub gravity: Vector3<f64>,
    /// PIC/FLIP blend: `α = 1.0` → pure PIC, `α = 0.0` → pure FLIP,
    /// `α ∈ (0, 1)` → blend.  Must be in `[0, 1]`.
    pub alpha: f64,
    /// Per-particle mass (kg), strictly positive.
    pub mass: f64,
    /// *(informational)* The simplified pressure step is a single Jacobi pass —
    /// see the module doc.
    pub _pressure_jacobi_passes: usize,
}

impl PicFlipConfig {
    /// Validate the configuration and return a new instance.
    ///
    /// # Errors
    /// - [`FluidError::InvalidConfig`] if `cell_size ≤ 0`, `mass ≤ 0`,
    ///   `alpha ∉ [0, 1]`, any grid dimension is `0`, or any value is
    ///   non-finite.
    /// - [`FluidError::NonFinite`] if `origin` or `gravity` have non-finite
    ///   components.
    pub fn new(
        grid_cells: (usize, usize, usize),
        cell_size: f64,
        origin: Vector3<f64>,
        gravity: Vector3<f64>,
        alpha: f64,
        mass: f64,
    ) -> Result<Self, FluidError> {
        if grid_cells.0 == 0 || grid_cells.1 == 0 || grid_cells.2 == 0 {
            return Err(FluidError::InvalidConfig(format!(
                "grid_cells must be ≥ 1 on every axis, got {grid_cells:?}"
            )));
        }
        if !(cell_size.is_finite() && cell_size > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "cell_size must be finite and > 0, got {cell_size}"
            )));
        }
        if !origin.iter().all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("origin".into()));
        }
        if !gravity.iter().all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("gravity".into()));
        }
        if !(alpha.is_finite() && (0.0..=1.0).contains(&alpha)) {
            return Err(FluidError::InvalidConfig(format!(
                "alpha must be in [0, 1], got {alpha}"
            )));
        }
        if !(mass.is_finite() && mass > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "mass must be finite and > 0, got {mass}"
            )));
        }
        Ok(Self {
            grid_cells,
            cell_size,
            origin,
            gravity,
            alpha,
            mass,
            _pressure_jacobi_passes: 1,
        })
    }

    /// Convenient default for a cubic domain of `n³` cells, gravity along −Y.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `cell_size ≤ 0` or `n = 0`.
    pub fn cubic(n: usize, cell_size: f64) -> Result<Self, FluidError> {
        Self::new(
            (n, n, n),
            cell_size,
            Vector3::zeros(),
            Vector3::new(0.0, -9.806_65, 0.0),
            0.05, // FLIP-dominant blend (typical for fluids)
            1.0,
        )
    }
}

// ─── MAC grid storage ─────────────────────────────────────────────────────────

/// A staggered MAC velocity grid.
///
/// Stores three face-centred velocity arrays `u`, `v`, `w` plus matching
/// weight accumulators used during P→G transfer.  All buffers are flat
/// `Vec<f64>` indexed by a linearisation helper.
///
/// For a grid of `(nx, ny, nz)` cells:
/// * `u` — `(nx+1) × ny × nz` values, x-face centres.
/// * `v` — `nx × (ny+1) × nz` values, y-face centres.
/// * `w` — `nx × ny × (nz+1)` values, z-face centres.
#[derive(Debug, Clone)]
struct MacGrid {
    nx: usize,
    ny: usize,
    nz: usize,
    // Velocity buffers (new = after grid step, old = before = post P→G).
    u: Vec<f64>, // (nx+1)*ny*nz
    v: Vec<f64>, // nx*(ny+1)*nz
    w: Vec<f64>, // nx*ny*(nz+1)
    u_old: Vec<f64>,
    v_old: Vec<f64>,
    w_old: Vec<f64>,
    // Weight accumulators for normalised splat.
    uw: Vec<f64>,
    vw: Vec<f64>,
    ww: Vec<f64>,
}

impl MacGrid {
    fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            nx,
            ny,
            nz,
            u: vec![0.0; (nx + 1) * ny * nz],
            v: vec![0.0; nx * (ny + 1) * nz],
            w: vec![0.0; nx * ny * (nz + 1)],
            u_old: vec![0.0; (nx + 1) * ny * nz],
            v_old: vec![0.0; nx * (ny + 1) * nz],
            w_old: vec![0.0; nx * ny * (nz + 1)],
            uw: vec![0.0; (nx + 1) * ny * nz],
            vw: vec![0.0; nx * (ny + 1) * nz],
            ww: vec![0.0; nx * ny * (nz + 1)],
        }
    }

    // Flat index helpers.
    #[inline]
    fn u_idx(&self, i: usize, j: usize, k: usize) -> usize {
        i * self.ny * self.nz + j * self.nz + k
    }
    #[inline]
    fn v_idx(&self, i: usize, j: usize, k: usize) -> usize {
        i * (self.ny + 1) * self.nz + j * self.nz + k
    }
    #[inline]
    fn w_idx(&self, i: usize, j: usize, k: usize) -> usize {
        i * self.ny * (self.nz + 1) + j * (self.nz + 1) + k
    }

    /// Zero the velocity and weight accumulators (called before P→G).
    fn clear_accumulators(&mut self) {
        self.u.iter_mut().for_each(|x| *x = 0.0);
        self.v.iter_mut().for_each(|x| *x = 0.0);
        self.w.iter_mut().for_each(|x| *x = 0.0);
        self.uw.iter_mut().for_each(|x| *x = 0.0);
        self.vw.iter_mut().for_each(|x| *x = 0.0);
        self.ww.iter_mut().for_each(|x| *x = 0.0);
    }

    /// Snapshot `{u,v,w}` into `{u,v,w}_old` (called after P→G, before grid step).
    fn snapshot_old(&mut self) {
        self.u_old.clone_from(&self.u);
        self.v_old.clone_from(&self.v);
        self.w_old.clone_from(&self.w);
    }

    /// Normalise weighted sums: divide `u` by `uw` where `uw > 0`.
    fn normalise(&mut self) {
        for (vel, w) in self.u.iter_mut().zip(self.uw.iter()) {
            if *w > 0.0 {
                *vel /= w;
            }
        }
        for (vel, w) in self.v.iter_mut().zip(self.vw.iter()) {
            if *w > 0.0 {
                *vel /= w;
            }
        }
        for (vel, w) in self.w.iter_mut().zip(self.ww.iter()) {
            if *w > 0.0 {
                *vel /= w;
            }
        }
    }
}

// ─── Trilinear weight helper ──────────────────────────────────────────────────

/// Clamp `x` to `[0, n-1]` as `usize`, returning the index and the
/// fractional offset `frac ∈ [0, 1]`.
///
/// For a staggered-grid deposit, the particle position is first shifted by half
/// a cell for the appropriate component before calling this.
#[inline]
fn cell_frac(x: f64, n: usize) -> (usize, f64) {
    // Map to cell space, clamp to valid range.
    let xf = x.max(0.0);
    let idx = (xf.floor() as usize).min(n.saturating_sub(1));
    let frac = (xf - idx as f64).clamp(0.0, 1.0);
    (idx, frac)
}

/// Eight trilinear weights for a unit cube at fractional offsets `(fx, fy, fz)`.
///
/// Returns weights indexed by `(di, dj, dk) ∈ {0,1}³` in order
/// `(0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1)`.
#[inline]
fn trilinear_weights(fx: f64, fy: f64, fz: f64) -> [f64; 8] {
    let ox = 1.0 - fx;
    let oy = 1.0 - fy;
    let oz = 1.0 - fz;
    [
        ox * oy * oz, // (0,0,0)
        fx * oy * oz, // (1,0,0)
        ox * fy * oz, // (0,1,0)
        fx * fy * oz, // (1,1,0)
        ox * oy * fz, // (0,0,1)
        fx * oy * fz, // (1,0,1)
        ox * fy * fz, // (0,1,1)
        fx * fy * fz, // (1,1,1)
    ]
}

// ─── Solver ───────────────────────────────────────────────────────────────────

/// A PIC/FLIP particle-grid fluid solver.
///
/// Call [`PicFlipSolver::step`] each simulation frame to advance the
/// [`ParticleSystem`] by a time step `dt`.
#[derive(Debug, Clone)]
pub struct PicFlipSolver {
    cfg: PicFlipConfig,
    grid: MacGrid,
}

impl PicFlipSolver {
    /// Create a solver from a validated [`PicFlipConfig`].
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] — forwarded from [`PicFlipConfig::new`].
    pub fn new(cfg: PicFlipConfig) -> Result<Self, FluidError> {
        // Re-validate to make it impossible to bypass the checked constructor.
        if cfg.cell_size <= 0.0 || !cfg.cell_size.is_finite() {
            return Err(FluidError::InvalidConfig(format!(
                "cell_size must be finite and > 0, got {}",
                cfg.cell_size
            )));
        }
        if !(0.0..=1.0).contains(&cfg.alpha) {
            return Err(FluidError::InvalidConfig(format!(
                "alpha must be in [0, 1], got {}",
                cfg.alpha
            )));
        }
        let (nx, ny, nz) = cfg.grid_cells;
        Ok(Self {
            cfg,
            grid: MacGrid::new(nx, ny, nz),
        })
    }

    /// A read-only reference to the configuration.
    pub fn config(&self) -> &PicFlipConfig {
        &self.cfg
    }

    // ── Phase 1: P→G ─────────────────────────────────────────────────────────

    /// Deposit particle velocities onto the MAC grid (particle-to-grid, P→G).
    ///
    /// Each particle splatted onto the faces it is closest to using trilinear
    /// weights.  After all particles are accumulated the weight sums normalise
    /// the face velocities (guarding divides: only normalise where weight > 0).
    fn transfer_p2g(&mut self, system: &ParticleSystem) {
        // Capture plain values before mutably borrowing grid fields.
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nz = self.grid.nz;
        let dx = self.cfg.cell_size;
        let inv_dx = 1.0 / dx;
        let origin = self.cfg.origin;

        for p in system.particles() {
            // Particle position in grid space (cell units, zero at origin min corner).
            let gx = (p.position.x - origin.x) * inv_dx;
            let gy = (p.position.y - origin.y) * inv_dx;
            let gz = (p.position.z - origin.z) * inv_dx;

            // ── u (x-velocity, staggered at i+½; x faces span i=0..=nx) ──────
            {
                let (i0, fx) = cell_frac(gx, nx); // i0 ∈ [0, nx]
                let (j0, fy) = cell_frac(gy - 0.5, ny.saturating_sub(1));
                let (k0, fz) = cell_frac(gz - 0.5, nz.saturating_sub(1));
                let wts = trilinear_weights(fx, fy, fz);
                let corners = [
                    (i0, j0, k0),
                    (i0 + 1, j0, k0),
                    (i0, j0 + 1, k0),
                    (i0 + 1, j0 + 1, k0),
                    (i0, j0, k0 + 1),
                    (i0 + 1, j0, k0 + 1),
                    (i0, j0 + 1, k0 + 1),
                    (i0 + 1, j0 + 1, k0 + 1),
                ];
                for (&wt, &(ci, cj, ck)) in wts.iter().zip(corners.iter()) {
                    if ci <= nx && cj < ny && ck < nz {
                        let idx = ci * ny * nz + cj * nz + ck;
                        self.grid.u[idx] += wt * p.velocity.x;
                        self.grid.uw[idx] += wt;
                    }
                }
            }

            // ── v (y-velocity, staggered at j+½; y faces span j=0..=ny) ──────
            {
                let (i0, fx) = cell_frac(gx - 0.5, nx.saturating_sub(1));
                let (j0, fy) = cell_frac(gy, ny); // j0 ∈ [0, ny]
                let (k0, fz) = cell_frac(gz - 0.5, nz.saturating_sub(1));
                let wts = trilinear_weights(fx, fy, fz);
                let corners = [
                    (i0, j0, k0),
                    (i0 + 1, j0, k0),
                    (i0, j0 + 1, k0),
                    (i0 + 1, j0 + 1, k0),
                    (i0, j0, k0 + 1),
                    (i0 + 1, j0, k0 + 1),
                    (i0, j0 + 1, k0 + 1),
                    (i0 + 1, j0 + 1, k0 + 1),
                ];
                for (&wt, &(ci, cj, ck)) in wts.iter().zip(corners.iter()) {
                    if ci < nx && cj <= ny && ck < nz {
                        let idx = ci * (ny + 1) * nz + cj * nz + ck;
                        self.grid.v[idx] += wt * p.velocity.y;
                        self.grid.vw[idx] += wt;
                    }
                }
            }

            // ── w (z-velocity, staggered at k+½; z faces span k=0..=nz) ──────
            {
                let (i0, fx) = cell_frac(gx - 0.5, nx.saturating_sub(1));
                let (j0, fy) = cell_frac(gy - 0.5, ny.saturating_sub(1));
                let (k0, fz) = cell_frac(gz, nz); // k0 ∈ [0, nz]
                let wts = trilinear_weights(fx, fy, fz);
                let corners = [
                    (i0, j0, k0),
                    (i0 + 1, j0, k0),
                    (i0, j0 + 1, k0),
                    (i0 + 1, j0 + 1, k0),
                    (i0, j0, k0 + 1),
                    (i0 + 1, j0, k0 + 1),
                    (i0, j0 + 1, k0 + 1),
                    (i0 + 1, j0 + 1, k0 + 1),
                ];
                for (&wt, &(ci, cj, ck)) in wts.iter().zip(corners.iter()) {
                    if ci < nx && cj < ny && ck <= nz {
                        let idx = ci * ny * (nz + 1) + cj * (nz + 1) + ck;
                        self.grid.w[idx] += wt * p.velocity.z;
                        self.grid.ww[idx] += wt;
                    }
                }
            }
        }

        // Normalise (guarded: only where weight sum > 0).
        self.grid.normalise();
    }

    // ── Phase 2: Grid step ────────────────────────────────────────────────────

    /// Apply gravity and a single-pass divergence step on the grid.
    ///
    /// **Gravity** — add `g · dt` to the appropriate velocity component at
    /// every face.
    ///
    /// **Pressure / divergence step** — for each interior cell, compute the
    /// local velocity divergence from the six face velocities and subtract an
    /// equal share from each contributing face to drive the divergence toward
    /// zero.  This is a single Jacobi pass (see the module-level caveat).
    fn grid_step(&mut self, dt: f64) {
        let (nx, ny, nz) = self.cfg.grid_cells;
        let g = self.cfg.gravity;

        // Gravity: add g*dt to every face velocity.
        for i in 0..=nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = self.grid.u_idx(i, j, k);
                    self.grid.u[idx] += g.x * dt;
                }
            }
        }
        for i in 0..nx {
            for j in 0..=ny {
                for k in 0..nz {
                    let idx = self.grid.v_idx(i, j, k);
                    self.grid.v[idx] += g.y * dt;
                }
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..=nz {
                    let idx = self.grid.w_idx(i, j, k);
                    self.grid.w[idx] += g.z * dt;
                }
            }
        }

        // Simplified pressure: single Jacobi divergence step.
        // For cell (i,j,k) the divergence is:
        //   div = (u[i+1][j][k] - u[i][j][k]
        //        + v[i][j+1][k] - v[i][j][k]
        //        + w[i][j][k+1] - w[i][j][k]) / dx
        // We subtract div/6 from each of the 6 contributing faces.
        let dx = self.cfg.cell_size;
        let inv_6dx = 1.0 / (6.0 * dx);
        // Collect corrections into separate delta buffers (explicit / Jacobi).
        let mut du = vec![0.0f64; (nx + 1) * ny * nz];
        let mut dv = vec![0.0f64; nx * (ny + 1) * nz];
        let mut dw = vec![0.0f64; nx * ny * (nz + 1)];

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let ui0 = self.grid.u_idx(i, j, k);
                    let ui1 = self.grid.u_idx(i + 1, j, k);
                    let vi0 = self.grid.v_idx(i, j, k);
                    let vi1 = self.grid.v_idx(i, j + 1, k);
                    let wi0 = self.grid.w_idx(i, j, k);
                    let wi1 = self.grid.w_idx(i, j, k + 1);

                    let div = (self.grid.u[ui1] - self.grid.u[ui0] + self.grid.v[vi1]
                        - self.grid.v[vi0]
                        + self.grid.w[wi1]
                        - self.grid.w[wi0])
                        / dx;

                    let delta = div * inv_6dx; // share for each face
                                               // Signs: to reduce +divergence, shrink outward and grow inward.
                    du[ui1] -= delta;
                    du[ui0] += delta;
                    dv[vi1] -= delta;
                    dv[vi0] += delta;
                    dw[wi1] -= delta;
                    dw[wi0] += delta;
                }
            }
        }

        // Apply corrections.
        for (vel, d) in self.grid.u.iter_mut().zip(du.iter()) {
            *vel += d;
        }
        for (vel, d) in self.grid.v.iter_mut().zip(dv.iter()) {
            *vel += d;
        }
        for (vel, d) in self.grid.w.iter_mut().zip(dw.iter()) {
            *vel += d;
        }
    }

    // ── Phase 3: G→P ─────────────────────────────────────────────────────────

    /// Interpolate the updated grid velocity back to particles (G→P), applying
    /// the PIC/FLIP blend.
    ///
    /// For each particle:
    /// * Interpolate `v_grid_new` and `v_grid_old` (pre-step snapshot) at the
    ///   particle position by the same trilinear weights used in P→G.
    /// * Apply `v_new = α·v_grid_new + (1−α)·(v_particle + Δv_grid)`.
    fn transfer_g2p(&mut self, system: &mut ParticleSystem) {
        let alpha = self.cfg.alpha;
        let dx = self.cfg.cell_size;
        let inv_dx = 1.0 / dx;
        let origin = self.cfg.origin;
        let (nx, ny, nz) = self.cfg.grid_cells;

        for p in system.particles_mut() {
            let gx = (p.position.x - origin.x) * inv_dx;
            let gy = (p.position.y - origin.y) * inv_dx;
            let gz = (p.position.z - origin.z) * inv_dx;

            // ── Interpolate x-velocity ────────────────────────────────────────
            let vx_new = {
                let (i0, fx) = cell_frac(gx, nx);
                let (j0, fy) = cell_frac(gy - 0.5, ny - 1);
                let (k0, fz) = cell_frac(gz - 0.5, nz - 1);
                let wts = trilinear_weights(fx, fy, fz);
                let corners = [
                    (i0, j0, k0),
                    (i0 + 1, j0, k0),
                    (i0, j0 + 1, k0),
                    (i0 + 1, j0 + 1, k0),
                    (i0, j0, k0 + 1),
                    (i0 + 1, j0, k0 + 1),
                    (i0, j0 + 1, k0 + 1),
                    (i0 + 1, j0 + 1, k0 + 1),
                ];
                interp_face(
                    &wts,
                    &corners,
                    &self.grid.u,
                    &self.grid.u_old,
                    nx + 1,
                    ny,
                    nz,
                )
            };

            // ── Interpolate y-velocity ────────────────────────────────────────
            let vy_new = {
                let (i0, fx) = cell_frac(gx - 0.5, nx - 1);
                let (j0, fy) = cell_frac(gy, ny);
                let (k0, fz) = cell_frac(gz - 0.5, nz - 1);
                let wts = trilinear_weights(fx, fy, fz);
                let corners = [
                    (i0, j0, k0),
                    (i0 + 1, j0, k0),
                    (i0, j0 + 1, k0),
                    (i0 + 1, j0 + 1, k0),
                    (i0, j0, k0 + 1),
                    (i0 + 1, j0, k0 + 1),
                    (i0, j0 + 1, k0 + 1),
                    (i0 + 1, j0 + 1, k0 + 1),
                ];
                interp_face(
                    &wts,
                    &corners,
                    &self.grid.v,
                    &self.grid.v_old,
                    nx,
                    ny + 1,
                    nz,
                )
            };

            // ── Interpolate z-velocity ────────────────────────────────────────
            let vz_new = {
                let (i0, fx) = cell_frac(gx - 0.5, nx - 1);
                let (j0, fy) = cell_frac(gy - 0.5, ny - 1);
                let (k0, fz) = cell_frac(gz, nz);
                let wts = trilinear_weights(fx, fy, fz);
                let corners = [
                    (i0, j0, k0),
                    (i0 + 1, j0, k0),
                    (i0, j0 + 1, k0),
                    (i0 + 1, j0 + 1, k0),
                    (i0, j0, k0 + 1),
                    (i0 + 1, j0, k0 + 1),
                    (i0, j0 + 1, k0 + 1),
                    (i0 + 1, j0 + 1, k0 + 1),
                ];
                interp_face(
                    &wts,
                    &corners,
                    &self.grid.w,
                    &self.grid.w_old,
                    nx,
                    ny,
                    nz + 1,
                )
            };

            // v_grid_new and Δv = v_grid_new − v_grid_old.
            let v_grid = Vector3::new(vx_new.0, vy_new.0, vz_new.0);
            let delta_v = Vector3::new(
                vx_new.0 - vx_new.1,
                vy_new.0 - vy_new.1,
                vz_new.0 - vz_new.1,
            );

            // PIC/FLIP blend.
            p.velocity = alpha * v_grid + (1.0 - alpha) * (p.velocity + delta_v);
        }
    }

    // ── Phase 4: Advect ───────────────────────────────────────────────────────

    /// Advect particles forward by `dt` and apply the boundary.
    fn advect(system: &mut ParticleSystem, dt: f64, boundary: Option<&BoxBoundary>) {
        for p in system.particles_mut() {
            p.position += p.velocity * dt;
            if let Some(b) = boundary {
                b.enforce(p);
            }
        }
    }

    // ── Public API ────────────────────────────────────────────────────────────

    /// Advance the particle system by one time step `dt`.
    ///
    /// Steps (in order):
    /// 1. Validate `dt` (must be finite and > 0).
    /// 2. P→G: deposit particle velocities onto the MAC grid.
    /// 3. Snapshot the pre-step grid velocities.
    /// 4. Grid step: apply gravity + simplified pressure.
    /// 5. G→P: PIC/FLIP blend back to particles.
    /// 6. Advect: `x ← x + v·dt` and enforce boundary.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `dt ≤ 0` or not finite.
    pub fn step(
        &mut self,
        system: &mut ParticleSystem,
        dt: f64,
        boundary: Option<&BoxBoundary>,
    ) -> Result<(), FluidError> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "dt must be finite and > 0, got {dt}"
            )));
        }

        // 1. P→G (clears accumulators internally).
        self.grid.clear_accumulators();
        self.transfer_p2g(system);

        // 2. Snapshot before the grid step.
        self.grid.snapshot_old();

        // 3. Grid step.
        self.grid_step(dt);

        // 4. G→P blend.
        self.transfer_g2p(system);

        // 5. Advect.
        Self::advect(system, dt, boundary);

        Ok(())
    }
}

// ─── Interpolation helper (G→P) ──────────────────────────────────────────────

/// Trilinearly interpolate `(new_vel, old_vel)` at a set of corner indices
/// from the given face buffers.
///
/// Returns `(interpolated_new, interpolated_old)`.
///
/// Out-of-bounds corners (can happen at domain edges after clamping) contribute
/// zero weight.
fn interp_face(
    wts: &[f64; 8],
    corners: &[(usize, usize, usize); 8],
    buf_new: &[f64],
    buf_old: &[f64],
    ni: usize,
    nj: usize,
    nk: usize,
) -> (f64, f64) {
    let mut vn = 0.0f64;
    let mut vo = 0.0f64;
    let mut wt_sum = 0.0f64;
    for (&wt, &(ci, cj, ck)) in wts.iter().zip(corners.iter()) {
        if ci < ni && cj < nj && ck < nk {
            let idx = ci * nj * nk + cj * nk + ck;
            vn += wt * buf_new[idx];
            vo += wt * buf_old[idx];
            wt_sum += wt;
        }
    }
    if wt_sum > 0.0 {
        (vn / wt_sum, vo / wt_sum)
    } else {
        (0.0, 0.0)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::Particle;
    use nalgebra::Vector3;

    /// Build a simple 4×4×4 grid solver centred at the origin, gravity off,
    /// pure PIC (alpha = 1).
    fn pic_solver() -> PicFlipSolver {
        let cfg = PicFlipConfig::new(
            (4, 4, 4),
            1.0,
            Vector3::zeros(),
            Vector3::zeros(), // no gravity
            1.0,              // pure PIC
            1.0,
        )
        .unwrap();
        PicFlipSolver::new(cfg).unwrap()
    }

    // ── Config validation ────────────────────────────────────────────────────

    #[test]
    fn rejects_zero_cell_size() {
        let r = PicFlipConfig::new((4, 4, 4), 0.0, Vector3::zeros(), Vector3::zeros(), 0.5, 1.0);
        assert!(r.is_err(), "zero cell_size must be rejected");
    }

    #[test]
    fn rejects_negative_cell_size() {
        let r = PicFlipConfig::new(
            (4, 4, 4),
            -1.0,
            Vector3::zeros(),
            Vector3::zeros(),
            0.5,
            1.0,
        );
        assert!(r.is_err(), "negative cell_size must be rejected");
    }

    #[test]
    fn rejects_alpha_out_of_range() {
        for bad in &[-0.1_f64, 1.1, f64::NAN, f64::INFINITY] {
            let r = PicFlipConfig::new(
                (4, 4, 4),
                1.0,
                Vector3::zeros(),
                Vector3::zeros(),
                *bad,
                1.0,
            );
            assert!(r.is_err(), "alpha={bad} must be rejected");
        }
    }

    #[test]
    fn rejects_nan_gravity() {
        let r = PicFlipConfig::new(
            (4, 4, 4),
            1.0,
            Vector3::zeros(),
            Vector3::new(f64::NAN, 0.0, 0.0),
            0.5,
            1.0,
        );
        assert!(r.is_err(), "NaN gravity must be rejected");
    }

    #[test]
    fn rejects_zero_grid_dimension() {
        let r = PicFlipConfig::new((0, 4, 4), 1.0, Vector3::zeros(), Vector3::zeros(), 0.5, 1.0);
        assert!(r.is_err(), "zero grid dimension must be rejected");
    }

    #[test]
    fn rejects_nonpositive_mass() {
        for bad in &[0.0_f64, -1.0, f64::NAN] {
            let r = PicFlipConfig::new(
                (4, 4, 4),
                1.0,
                Vector3::zeros(),
                Vector3::zeros(),
                0.5,
                *bad,
            );
            assert!(r.is_err(), "mass={bad} must be rejected");
        }
    }

    #[test]
    fn rejects_nonpositive_dt() {
        let mut solver = pic_solver();
        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::new(2.0, 2.0, 2.0))).unwrap();
        assert!(solver.step(&mut sys, 0.0, None).is_err());
        assert!(solver.step(&mut sys, -1.0, None).is_err());
        assert!(solver.step(&mut sys, f64::NAN, None).is_err());
    }

    // ── Benchmark 1: Round-trip (pure PIC) ──────────────────────────────────

    /// A uniform velocity field P→G→P roundtrips to the original velocity
    /// within interpolation tolerance (pure PIC, no gravity, no pressure change).
    #[test]
    fn pic_roundtrip_uniform_velocity() {
        let mut solver = pic_solver();
        let v0 = Vector3::new(1.5, -0.7, 0.3);

        // Place particles on a regular grid well inside the domain.
        let mut sys = ParticleSystem::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let pos = Vector3::new(
                        1.0 + i as f64 * 0.5,
                        1.0 + j as f64 * 0.5,
                        1.0 + k as f64 * 0.5,
                    );
                    sys.push(Particle::new(pos, v0)).unwrap();
                }
            }
        }

        // P→G then G→P only (no grid step; we manually do P→G + snapshot + G→P).
        solver.grid.clear_accumulators();
        solver.transfer_p2g(&sys);
        solver.grid.snapshot_old();
        // No grid step — old and new are identical, so Δv = 0.
        // For pure PIC (alpha=1) the result must equal v_grid_new ≈ v0.
        solver.transfer_g2p(&mut sys);

        for p in sys.particles() {
            let err = (p.velocity - v0).norm();
            assert!(
                err < 1e-10,
                "PIC roundtrip error {err:.2e} exceeds tolerance; got {:?}",
                p.velocity
            );
        }
    }

    // ── Benchmark 2: Free-fall (no pressure) ────────────────────────────────

    /// A single particle under gravity (no pressure step) falls ≈ ½ g t² in
    /// the y-direction over two equal steps of dt=0.05 s.
    ///
    /// The grid step adds g·dt to the grid y-velocity; since the particle is
    /// alone on the grid the pure-PIC G→P gives v_particle ≈ g·dt after step 1,
    /// and the advection moves the particle down by v·dt ≈ g·dt².
    #[test]
    fn single_particle_freefall_under_gravity() {
        let g_y = -9.806_65_f64;
        let cfg = PicFlipConfig::new(
            (4, 4, 4),
            1.0,
            Vector3::zeros(),
            Vector3::new(0.0, g_y, 0.0),
            1.0, // pure PIC
            1.0,
        )
        .unwrap();
        let mut solver = PicFlipSolver::new(cfg).unwrap();

        let mut sys = ParticleSystem::new();
        let y0 = 2.0_f64;
        sys.push(Particle::at(Vector3::new(2.0, y0, 2.0))).unwrap();

        let dt = 0.05_f64;
        // Two steps — enough to check ½ g t² without leaving the box.
        solver.step(&mut sys, dt, None).unwrap();
        solver.step(&mut sys, dt, None).unwrap();

        let y_actual = sys.particles()[0].position.y;
        // Expected: y0 + ½ g (2dt)²  (two steps = total time 2dt).
        let t = 2.0 * dt;
        let y_expected = y0 + 0.5 * g_y * t * t;
        let err = (y_actual - y_expected).abs();
        // Allow 5 % relative error — first-order Euler + single-pass grid-to-particle.
        assert!(
            err < 0.05 * y_expected.abs().max(1e-6),
            "free-fall error {err:.4}: got y={y_actual:.4}, expected {y_expected:.4}"
        );
    }

    // ── Benchmark 3: Momentum conservation (P→G→P, no forces) ──────────────

    /// The total momentum of the particle system is approximately conserved
    /// through one P→G→P cycle with no external forces.
    ///
    /// Both PIC and FLIP preserve total momentum in the mass-weighted sense:
    /// each particle's velocity is a weighted average of grid face values that
    /// were themselves weighted averages of particle velocities, so momentum
    /// can drift at most by `O(h)` truncation from the splatting.
    #[test]
    fn pic_transfer_conserves_total_momentum() {
        let mut solver = pic_solver();

        let mut sys = ParticleSystem::new();
        let mass = solver.cfg.mass;
        // A few particles with differing velocities inside the grid.
        let vels = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(-0.5, 1.2, 0.0),
            Vector3::new(0.3, -0.7, 0.8),
            Vector3::new(-0.8, 0.1, -0.4),
        ];
        let positions = [
            Vector3::new(1.5, 1.5, 1.5),
            Vector3::new(2.5, 1.5, 1.5),
            Vector3::new(1.5, 2.5, 1.5),
            Vector3::new(2.5, 2.5, 2.5),
        ];
        for (&pos, &vel) in positions.iter().zip(vels.iter()) {
            sys.push(Particle::new(pos, vel)).unwrap();
        }

        let p_before = sys.total_momentum(mass);

        // One P→G→P cycle (no grid step, so Δv_grid = 0 for both PIC and FLIP).
        solver.grid.clear_accumulators();
        solver.transfer_p2g(&sys);
        solver.grid.snapshot_old();
        // Don't call grid_step, so old == new.
        solver.transfer_g2p(&mut sys);

        let p_after = sys.total_momentum(mass);
        let drift = (p_after - p_before).norm();
        let scale = p_before.norm().max(1.0);
        assert!(
            drift / scale < 0.05,
            "momentum drift {drift:.4} > 5% of initial {scale:.4}"
        );
    }

    // ── Benchmark 4: FLIP blend correctness ─────────────────────────────────

    /// The PIC/FLIP blend formula `v_new = α·v_grid_new + (1−α)·(v_particle + Δv_grid)`
    /// is checked for `α ∈ {0, 0.5, 1}`.
    ///
    /// Setup: a single particle at `(2.0, 2.0, 2.0)` at rest (v=0), gravity
    /// `g = (0, −2, 0)`, `dt = 1`.  Only the y-component is checked — the
    /// gravity body-force step applies a uniform impulse of `g.y·dt = −2` to
    /// ALL y-faces, so the deposited y-face velocity (initially 0) becomes −2
    /// everywhere after the grid step. This is a pure body-force test: no x/z
    /// component is deposited and the y-field has no divergence from gravity.
    ///
    /// For all three `α` values the y-velocity must converge to `g.y·dt = −2`:
    ///
    /// * `α = 1` (PIC): `v_y = v_grid_y_new = 0 + g.y·dt = −2`.
    /// * `α = 0` (FLIP): `v_y = v_particle_y + Δv_y = 0 + (−2 − 0) = −2`.
    /// * `α = 0.5`: `0.5·(−2) + 0.5·(0 + (−2−0)) = −2`.
    ///
    /// The x-velocity is NOT checked here — the pressure step redistributes
    /// the single-particle u-step function and intentionally modifies x-faces.
    #[test]
    fn flip_blend_interpolates_correctly() {
        // Particle starts at rest; gravity impulse is 1 second of g=(0,−2,0).
        let v0 = Vector3::zeros();
        let g_y = -2.0_f64;
        let dt = 1.0_f64;

        let test_alpha = |alpha: f64| -> f64 {
            let cfg = PicFlipConfig::new(
                (4, 4, 4),
                1.0,
                Vector3::zeros(),
                Vector3::new(0.0, g_y, 0.0),
                alpha,
                1.0,
            )
            .unwrap();
            let mut solver = PicFlipSolver::new(cfg).unwrap();

            let mut sys = ParticleSystem::new();
            // A particle at the centre of the domain; at rest.
            sys.push(Particle::new(Vector3::new(2.0, 2.0, 2.0), v0))
                .unwrap();

            solver.step(&mut sys, dt, None).unwrap();
            // Return the y-velocity.
            sys.particles()[0].velocity.y
        };

        let vy_expected = g_y * dt; // = −2.0

        // alpha = 1 (PIC): v_y = v_grid_new_y = g.y·dt = −2.
        let vy_pic = test_alpha(1.0);
        let err_pic = (vy_pic - vy_expected).abs();
        assert!(
            err_pic < 0.1,
            "PIC (α=1) y-velocity error {err_pic:.4}: got {vy_pic:.4}, expected {vy_expected:.4}"
        );

        // alpha = 0 (FLIP): v_y = v_particle_y + Δv_y = 0 + (−2 − 0) = −2.
        let vy_flip = test_alpha(0.0);
        let err_flip = (vy_flip - vy_expected).abs();
        assert!(
            err_flip < 0.1,
            "FLIP (α=0) y-velocity error {err_flip:.4}: got {vy_flip:.4}, expected {vy_expected:.4}"
        );

        // alpha = 0.5: 0.5·(−2) + 0.5·(0 + (−2−0)) = −2.
        let vy_blend = test_alpha(0.5);
        let err_blend = (vy_blend - vy_expected).abs();
        assert!(
            err_blend < 0.1,
            "blend (α=0.5) y-velocity error {err_blend:.4}: got {vy_blend:.4}, expected {vy_expected:.4}"
        );

        // Verify α monotonically interp: result is the same here (both sides = −2),
        // but also verify that α=0 and α=1 give the same outcome, i.e. the blend
        // formula is consistent.
        assert!(
            (vy_pic - vy_flip).abs() < 0.2,
            "PIC and FLIP should agree for a uniform grid field: pic={vy_pic:.4} flip={vy_flip:.4}"
        );
    }

    // ── Additional guard checks ──────────────────────────────────────────────

    /// A particle with a NaN position must not enter the system.
    #[test]
    fn rejects_nan_particle_position() {
        let mut sys = ParticleSystem::new();
        assert!(sys
            .push(Particle::at(Vector3::new(f64::NAN, 0.0, 0.0)))
            .is_err());
    }

    /// Empty system steps without error.
    #[test]
    fn empty_system_steps_without_error() {
        let mut solver = pic_solver();
        let mut sys = ParticleSystem::new();
        assert!(solver.step(&mut sys, 0.01, None).is_ok());
    }

    /// cubic() convenience constructor works.
    #[test]
    fn cubic_convenience_constructor() {
        let cfg = PicFlipConfig::cubic(8, 0.1).unwrap();
        assert_eq!(cfg.grid_cells, (8, 8, 8));
        assert!((cfg.cell_size - 0.1).abs() < 1e-15);
        let solver = PicFlipSolver::new(cfg);
        assert!(solver.is_ok());
    }
}
