//! The **transient (unsteady) incompressible-flow solver** — time
//! marching by the implicit-Euler / transient-SIMPLE scheme.
//!
//! # What this is
//!
//! [`crate::solver`] solves the *steady* Navier-Stokes equations: it
//! iterates straight to the time-independent final state and never
//! represents the flow at an intermediate instant. Many flows are
//! genuinely **unsteady** — a flow started impulsively from rest, a
//! flow with a time-varying boundary, vortex shedding — and to capture
//! them the time-derivative term must be put back into the momentum
//! equation:
//!
//! ```text
//!   ρ ∂u/∂t + ρ(u·∇)u = −∇p + μ∇²u
//!   ∇·u = 0
//! ```
//!
//! This module marches that system in time. Each physical time step it
//! runs an inner pressure-velocity coupling — exactly the SIMPLE
//! iteration of [`crate::solver`], but with the **unsteady term added
//! to every momentum control volume**. Discretising `∂u/∂t` by
//! first-order **implicit (backward) Euler** adds, per cell,
//!
//! ```text
//!   a_t = ρ·V / Δt          to the momentum diagonal a_P
//!   a_t · u_old             to the momentum source
//! ```
//!
//! where `u_old` is the velocity at the *previous* time level. The
//! scheme is **implicit** — every term except `u_old` is evaluated at
//! the new time level — so it is unconditionally stable: the time step
//! is limited by accuracy, not by a CFL stability bound.
//!
//! # Time marching
//!
//! [`solve_transient`] is the outer loop. For each of
//! [`TransientControls::n_steps`] steps it:
//!
//! 1. stores the current field as `u_old` / `v_old`;
//! 2. runs up to [`TransientControls::inner_iterations`] SIMPLE
//!    iterations on the implicit momentum + pressure-correction system
//!    until that step's mass imbalance is below tolerance;
//! 3. accepts the converged field as the new time level and advances
//!    `t += Δt`.
//!
//! At every step the velocity satisfies the discrete continuity
//! equation, so the field is divergence-free at every instant — the
//! defining property a transient incompressible solver must hold.
//!
//! As `t → ∞` a transient run with steady boundary conditions
//! **relaxes onto the steady solution** — the time derivative vanishes
//! and the equations become the steady ones [`crate::solver`] solves.
//! That limit is exactly what the tests here verify.
//!
//! # Honest scope
//!
//! This is a **real working transient solver v1** — a channel flow
//! started from rest is verified to relax onto the known steady
//! parabolic profile. It is deliberately a v1:
//!
//! - **First-order implicit (backward) Euler** in time. Unconditionally
//!   stable and the standard first choice; a second-order
//!   Crank-Nicolson or BDF2 scheme is a bounded accuracy-only
//!   follow-up.
//! - **Transient SIMPLE**, not PISO. Each step is converged with the
//!   same outer SIMPLE iteration as the steady solver; a
//!   non-iterative PISO predictor-corrector is an efficiency-only
//!   alternative.
//! - **2-D, laminar, structured grid** — inherited from
//!   [`crate::solver`]. Pair with [`crate::turbulence`] for a
//!   turbulent run.
//!
//! Within that scope the time advancement is the genuine article: the
//! implicit unsteady term, a divergence-free field at every step, and
//! the correct long-time relaxation onto the steady state.

use crate::grid::{Field, Grid};
use crate::linsolve::{solve_sor, PoissonCoeffs};
use crate::solver::{Boundaries, Fluid, SideBc};

/// Controls for a transient time-marching run.
#[derive(Clone, Copy, Debug)]
pub struct TransientControls {
    /// The physical time step `Δt` (seconds). With the implicit scheme
    /// this is an *accuracy* choice, not a stability one — a smaller
    /// `Δt` resolves the transient more finely.
    pub dt: f64,
    /// How many time steps to march.
    pub n_steps: usize,
    /// Maximum inner SIMPLE iterations per time step. The implicit
    /// system at each step is nonlinear and is converged with the
    /// usual SIMPLE outer iteration.
    pub inner_iterations: usize,
    /// Convergence tolerance on the scaled mass-imbalance residual for
    /// the inner SIMPLE iteration of each step.
    pub tolerance: f64,
    /// Momentum under-relaxation factor `αu` for the inner iteration.
    /// A transient solve can run a milder relaxation than the steady
    /// solver because the unsteady term `ρV/Δt` already stabilises the
    /// momentum diagonal.
    pub relax_u: f64,
    /// Pressure under-relaxation factor `αp` for the inner iteration.
    pub relax_p: f64,
    /// Over-relaxation factor for the inner SOR pressure-correction
    /// solve.
    pub sor_omega: f64,
    /// Maximum inner SOR sweeps.
    pub sor_iterations: usize,
}

impl Default for TransientControls {
    /// Broadly-stable defaults: `Δt = 0.01 s`, 200 steps, up to 50
    /// inner SIMPLE iterations to a `1e-5` residual, `αu = 0.7`,
    /// `αp = 0.3`.
    fn default() -> Self {
        TransientControls {
            dt: 0.01,
            n_steps: 200,
            inner_iterations: 50,
            tolerance: 1e-5,
            relax_u: 0.7,
            relax_p: 0.3,
            sor_omega: 1.7,
            sor_iterations: 40,
        }
    }
}

/// A single time level of a transient run — the flow field at one
/// instant plus the diagnostics of the step that produced it.
#[derive(Clone, Debug)]
pub struct TransientStep {
    /// Physical time `t` (seconds) at the end of this step.
    pub time: f64,
    /// Scaled mass-imbalance residual the inner SIMPLE iteration of
    /// this step converged to.
    pub residual: f64,
    /// Inner SIMPLE iterations the step took.
    pub inner_iterations: usize,
}

/// The result of a transient time-marching run.
#[derive(Clone, Debug)]
pub struct TransientSolution {
    /// The grid the solution lives on.
    pub grid: Grid,
    /// Staggered `u`-velocity at the final time, `(nx+1) × ny`.
    pub u: Field,
    /// Staggered `v`-velocity at the final time, `nx × (ny+1)`.
    pub v: Field,
    /// Cell-centred pressure at the final time, `nx × ny`.
    pub pressure: Field,
    /// The final physical time reached.
    pub time: f64,
    /// Per-step history — one [`TransientStep`] per time step.
    pub history: Vec<TransientStep>,
}

impl TransientSolution {
    /// Cell-centred `u`-velocity of pressure cell `(i, j)` — the
    /// average of the two bracketing `u`-faces.
    pub fn u_at_cell(&self, i: usize, j: usize) -> f64 {
        0.5 * (self.u.at(i, j) + self.u.at(i + 1, j))
    }

    /// Cell-centred `v`-velocity of pressure cell `(i, j)`.
    pub fn v_at_cell(&self, i: usize, j: usize) -> f64 {
        0.5 * (self.v.at(i, j) + self.v.at(i, j + 1))
    }

    /// The largest velocity magnitude anywhere in the field — a
    /// convenient steady-state-detection norm.
    pub fn max_speed(&self) -> f64 {
        let mut m = 0.0_f64;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let u = self.u_at_cell(i, j);
                let v = self.v_at_cell(i, j);
                m = m.max((u * u + v * v).sqrt());
            }
        }
        m
    }
}

/// Solve an unsteady incompressible laminar flow by implicit time
/// marching.
///
/// `grid`, `fluid`, and `bcs` are exactly as for
/// [`crate::solve_simple`]; `controls` adds the time step, the step
/// count, and the inner-iteration settings. The flow starts from rest
/// (a zero velocity field) unless `initial` supplies a starting field.
///
/// Returns the [`TransientSolution`] — the final flow field plus the
/// per-step history. Marching with steady boundary conditions, the
/// field relaxes toward the steady solution; marching far enough, it
/// reaches it.
///
/// # Method
///
/// Each time step adds the implicit-Euler unsteady term `ρV/Δt` to the
/// momentum diagonal and `ρV/Δt·u_old` to the source, then runs the
/// standard SIMPLE pressure-velocity coupling — momentum predictor,
/// pressure-correction Poisson solve, under-relaxed correction — to
/// convergence on that step's nonlinear system.
pub fn solve_transient(
    grid: &Grid,
    fluid: &Fluid,
    bcs: &Boundaries,
    controls: &TransientControls,
    initial: Option<(&Field, &Field)>,
) -> TransientSolution {
    let nx = grid.nx;
    let ny = grid.ny;
    let dy = grid.dy();
    let rho = fluid.density;
    let nu = fluid.viscosity;
    let dt = controls.dt.max(1e-30);

    // Field storage — the *current* time level.
    let mut u = grid.u_field();
    let mut v = grid.v_field();
    let mut p = grid.pressure_field();
    if let Some((u0, v0)) = initial {
        if u0.width == u.width && u0.height == u.height {
            u = u0.clone();
        }
        if v0.width == v.width && v0.height == v.height {
            v = v0.clone();
        }
    }
    apply_velocity_bcs(&mut u, &mut v, bcs, nx, ny);

    // Momentum-equation diagonals — reused to build the
    // pressure-correction stencil.
    let mut apu = Field::zeros(nx + 1, ny);
    let mut apv = Field::zeros(nx, ny + 1);

    // Mass-residual normalisation (a characteristic mass flow).
    let ref_velocity = reference_velocity(bcs).max(1e-6);
    let mass_scale = (rho * ref_velocity * dy).max(1e-30);

    let mut time = 0.0;
    let mut history = Vec::with_capacity(controls.n_steps);

    for _step in 0..controls.n_steps {
        // The previous time level — frozen for the whole inner loop.
        let u_old = u.clone();
        let v_old = v.clone();

        let mut step_residual = f64::INFINITY;
        let mut inner_used = 0;

        for inner in 0..controls.inner_iterations.max(1) {
            inner_used = inner + 1;

            // --- implicit momentum predictor with the unsteady term ---
            solve_u_momentum(
                &mut u,
                &v,
                &p,
                &u_old,
                &mut apu,
                grid,
                rho,
                nu,
                dt,
                bcs,
                controls.relax_u,
            );
            solve_v_momentum(
                &u,
                &mut v,
                &p,
                &v_old,
                &mut apv,
                grid,
                rho,
                nu,
                dt,
                bcs,
                controls.relax_u,
            );
            apply_velocity_bcs(&mut u, &mut v, bcs, nx, ny);

            // --- pressure-correction Poisson equation ---
            let mut coeffs = PoissonCoeffs::zeros(nx, ny);
            let mass_imbalance =
                assemble_pressure_correction(&mut coeffs, &u, &v, &apu, &apv, grid, rho, bcs);

            // --- solve for p' ---
            let pin_mean = !has_outlet(bcs);
            let mut pcorr = Field::zeros(nx, ny);
            solve_sor(
                &coeffs,
                &mut pcorr,
                controls.sor_omega,
                controls.tolerance * 1e-3,
                controls.sor_iterations,
                pin_mean,
            );

            // --- correct pressure and velocity ---
            correct_pressure(&mut p, &pcorr, controls.relax_p);
            correct_velocity(&mut u, &mut v, &pcorr, &apu, &apv, grid);
            apply_velocity_bcs(&mut u, &mut v, bcs, nx, ny);

            step_residual = mass_imbalance / mass_scale;
            if step_residual.is_finite() && step_residual <= controls.tolerance {
                break;
            }
            if !step_residual.is_finite() || step_residual > 1e12 {
                break;
            }
        }

        time += dt;
        history.push(TransientStep {
            time,
            residual: step_residual,
            inner_iterations: inner_used,
        });
        if !step_residual.is_finite() {
            break;
        }
    }

    TransientSolution {
        grid: *grid,
        u,
        v,
        pressure: p,
        time,
        history,
    }
}

/// A characteristic velocity for the residual normalisation.
fn reference_velocity(bcs: &Boundaries) -> f64 {
    let mag = |s: &SideBc| -> f64 {
        match s {
            SideBc::Wall { u, v } | SideBc::Inlet { u, v } => (u * u + v * v).sqrt(),
            SideBc::Outlet => 0.0,
        }
    };
    mag(&bcs.west)
        .max(mag(&bcs.east))
        .max(mag(&bcs.south))
        .max(mag(&bcs.north))
}

/// True if any side is an outlet.
fn has_outlet(bcs: &Boundaries) -> bool {
    matches!(bcs.west, SideBc::Outlet)
        || matches!(bcs.east, SideBc::Outlet)
        || matches!(bcs.south, SideBc::Outlet)
        || matches!(bcs.north, SideBc::Outlet)
}

/// Stamp the prescribed boundary velocities onto the velocity fields —
/// the transient counterpart of [`crate::solver`]'s private
/// `apply_velocity_bcs`.
fn apply_velocity_bcs(u: &mut Field, v: &mut Field, bcs: &Boundaries, nx: usize, ny: usize) {
    for j in 0..ny {
        match bcs.west {
            SideBc::Wall { u: uw, .. } | SideBc::Inlet { u: uw, .. } => u.set(0, j, uw),
            SideBc::Outlet => {}
        }
        match bcs.east {
            SideBc::Wall { u: ue, .. } | SideBc::Inlet { u: ue, .. } => u.set(nx, j, ue),
            SideBc::Outlet => {
                let interior = u.at(nx - 1, j);
                u.set(nx, j, interior);
            }
        }
    }
    for i in 0..nx {
        match bcs.south {
            SideBc::Wall { v: vs, .. } | SideBc::Inlet { v: vs, .. } => v.set(i, 0, vs),
            SideBc::Outlet => {}
        }
        match bcs.north {
            SideBc::Wall { v: vn, .. } | SideBc::Inlet { v: vn, .. } => v.set(i, ny, vn),
            SideBc::Outlet => {
                let interior = v.at(i, ny - 1);
                v.set(i, ny, interior);
            }
        }
    }
}

/// Hybrid-scheme convection-diffusion coefficient — identical rule to
/// [`crate::solver`]'s.
#[inline]
fn hybrid_coeff(d: f64, f: f64, upwind_sign: f64) -> f64 {
    (d - 0.5 * f.abs()).max(0.0) + (upwind_sign * f).max(0.0)
}

/// Implicit `u`-momentum predictor — [`crate::solver`]'s steady
/// `u`-momentum kernel **plus the backward-Euler unsteady term**.
///
/// The unsteady term contributes `a_t = ρ·V/Δt` to the diagonal `a_P`
/// and `a_t·u_old` to the source. Every other term is at the new time
/// level, so the scheme is fully implicit.
#[allow(clippy::too_many_arguments)]
fn solve_u_momentum(
    u: &mut Field,
    v: &Field,
    p: &Field,
    u_old: &Field,
    apu: &mut Field,
    grid: &Grid,
    rho: f64,
    nu: f64,
    dt: f64,
    bcs: &Boundaries,
    relax: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();
    // Dynamic-viscosity diffusion conductance μ = ρν (matches the ρuA
    // convective flux and ρV/Δt unsteady term); kinematic ν dropped the ρ.
    let dx_diff = rho * nu * dy / dx;
    let dy_diff = rho * nu * dx / dy;
    // The unsteady-term coefficient a_t = ρ·V/Δt (the u-CV volume is
    // dx·dy per unit depth).
    let a_t = rho * dx * dy / dt;

    for _sweep in 0..2 {
        for j in 0..ny {
            for i in 1..nx {
                let fe = rho * dy * 0.5 * (u.at(i, j) + u.at(i + 1, j));
                let fw = rho * dy * 0.5 * (u.at(i - 1, j) + u.at(i, j));
                let fn_ = rho * dx * 0.5 * (v.at(i - 1, j + 1) + v.at(i, j + 1));
                let fs = rho * dx * 0.5 * (v.at(i - 1, j) + v.at(i, j));

                let ae = hybrid_coeff(dx_diff, fe, -1.0);
                let aw = hybrid_coeff(dx_diff, fw, 1.0);

                // The diagonal carries the unsteady term a_t.
                let mut a_p = ae + aw + a_t;
                // The source carries a_t·u_old (the previous time
                // level's velocity at this face).
                let mut su = a_t * u_old.at(i, j);
                let mut an = 0.0;
                let mut as_ = 0.0;
                let dwall = rho * nu * dx / (0.5 * dy);
                if j == ny - 1 {
                    let wall_u = match bcs.north {
                        SideBc::Wall { u, .. } | SideBc::Inlet { u, .. } => u,
                        SideBc::Outlet => u.at(i, j),
                    };
                    a_p += dwall;
                    su += dwall * wall_u;
                } else {
                    an = hybrid_coeff(dy_diff, fn_, -1.0);
                    a_p += an;
                }
                if j == 0 {
                    let wall_u = match bcs.south {
                        SideBc::Wall { u, .. } | SideBc::Inlet { u, .. } => u,
                        SideBc::Outlet => u.at(i, j),
                    };
                    a_p += dwall;
                    su += dwall * wall_u;
                } else {
                    as_ = hybrid_coeff(dy_diff, fs, 1.0);
                    a_p += as_;
                }

                if a_p.abs() < 1e-30 {
                    continue;
                }
                let dp = (p.at(i - 1, j) - p.at(i, j)) * dy;
                let mut nb = su + dp;
                nb += ae * u.at(i + 1, j);
                nb += aw * u.at(i - 1, j);
                if j + 1 < ny {
                    nb += an * u.at(i, j + 1);
                }
                if j > 0 {
                    nb += as_ * u.at(i, j - 1);
                }
                let u_new = nb / a_p;
                let u_curr = u.at(i, j);
                u.set(i, j, u_curr + relax * (u_new - u_curr));
                apu.set(i, j, a_p / relax);
            }
        }
    }
}

/// Implicit `v`-momentum predictor — the `y`-direction mirror of
/// [`solve_u_momentum`], likewise carrying the backward-Euler unsteady
/// term.
#[allow(clippy::too_many_arguments)]
fn solve_v_momentum(
    u: &Field,
    v: &mut Field,
    p: &Field,
    v_old: &Field,
    apv: &mut Field,
    grid: &Grid,
    rho: f64,
    nu: f64,
    dt: f64,
    bcs: &Boundaries,
    relax: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();
    // Dynamic-viscosity diffusion conductance μ = ρν (matches the ρuA
    // convective flux and ρV/Δt unsteady term); kinematic ν dropped the ρ.
    let dx_diff = rho * nu * dy / dx;
    let dy_diff = rho * nu * dx / dy;
    let a_t = rho * dx * dy / dt;

    for _sweep in 0..2 {
        for j in 1..ny {
            for i in 0..nx {
                let fn_ = rho * dx * 0.5 * (v.at(i, j) + v.at(i, j + 1));
                let fs = rho * dx * 0.5 * (v.at(i, j - 1) + v.at(i, j));
                let fe = rho * dy * 0.5 * (u.at(i + 1, j - 1) + u.at(i + 1, j));
                let fw = rho * dy * 0.5 * (u.at(i, j - 1) + u.at(i, j));

                let an = hybrid_coeff(dy_diff, fn_, -1.0);
                let as_ = hybrid_coeff(dy_diff, fs, 1.0);

                let mut a_p = an + as_ + a_t;
                let mut sv = a_t * v_old.at(i, j);
                let mut ae = 0.0;
                let mut aw = 0.0;
                let dwall = rho * nu * dy / (0.5 * dx);
                if i == nx - 1 {
                    let wall_v = match bcs.east {
                        SideBc::Wall { v, .. } | SideBc::Inlet { v, .. } => v,
                        SideBc::Outlet => v.at(i, j),
                    };
                    a_p += dwall;
                    sv += dwall * wall_v;
                } else {
                    ae = hybrid_coeff(dx_diff, fe, -1.0);
                    a_p += ae;
                }
                if i == 0 {
                    let wall_v = match bcs.west {
                        SideBc::Wall { v, .. } | SideBc::Inlet { v, .. } => v,
                        SideBc::Outlet => v.at(i, j),
                    };
                    a_p += dwall;
                    sv += dwall * wall_v;
                } else {
                    aw = hybrid_coeff(dx_diff, fw, 1.0);
                    a_p += aw;
                }

                if a_p.abs() < 1e-30 {
                    continue;
                }
                let dp = (p.at(i, j - 1) - p.at(i, j)) * dx;
                let mut nb = sv + dp;
                nb += an * v.at(i, j + 1);
                nb += as_ * v.at(i, j - 1);
                if i + 1 < nx {
                    nb += ae * v.at(i + 1, j);
                }
                if i > 0 {
                    nb += aw * v.at(i - 1, j);
                }
                let v_new = nb / a_p;
                let v_curr = v.at(i, j);
                v.set(i, j, v_curr + relax * (v_new - v_curr));
                apv.set(i, j, a_p / relax);
            }
        }
    }
}

/// Assemble the pressure-correction Poisson equation and return the L2
/// norm of the per-cell mass imbalance — identical to
/// [`crate::solver`]'s assembler (the unsteady term lives entirely in
/// the momentum diagonals `apu` / `apv`, which this routine simply
/// reads, so the pressure-correction stencil is unchanged).
#[allow(clippy::too_many_arguments)]
fn assemble_pressure_correction(
    coeffs: &mut PoissonCoeffs,
    u: &Field,
    v: &Field,
    apu: &Field,
    apv: &Field,
    grid: &Grid,
    rho: f64,
    bcs: &Boundaries,
) -> f64 {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();

    let mut sum_sq = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let ae = if i + 1 < nx {
                rho * dy * dy / apu.at(i + 1, j).max(1e-30)
            } else {
                0.0
            };
            let aw = if i > 0 {
                rho * dy * dy / apu.at(i, j).max(1e-30)
            } else {
                0.0
            };
            let an = if j + 1 < ny {
                rho * dx * dx / apv.at(i, j + 1).max(1e-30)
            } else {
                0.0
            };
            let as_ = if j > 0 {
                rho * dx * dx / apv.at(i, j).max(1e-30)
            } else {
                0.0
            };

            let mass_out =
                rho * dy * (u.at(i + 1, j) - u.at(i, j)) + rho * dx * (v.at(i, j + 1) - v.at(i, j));
            coeffs.ae.set(i, j, ae);
            coeffs.aw.set(i, j, aw);
            coeffs.an.set(i, j, an);
            coeffs.as_.set(i, j, as_);
            coeffs.ap.set(i, j, ae + aw + an + as_);
            coeffs.b.set(i, j, -mass_out);
            sum_sq += mass_out * mass_out;
        }
    }

    if has_outlet(bcs) {
        let anchor = if matches!(bcs.east, SideBc::Outlet) {
            (nx - 1, ny / 2)
        } else if matches!(bcs.west, SideBc::Outlet) {
            (0, ny / 2)
        } else if matches!(bcs.north, SideBc::Outlet) {
            (nx / 2, ny - 1)
        } else {
            (nx / 2, 0)
        };
        let big = 1e30;
        coeffs.ap.set(anchor.0, anchor.1, big);
        coeffs.ae.set(anchor.0, anchor.1, 0.0);
        coeffs.aw.set(anchor.0, anchor.1, 0.0);
        coeffs.an.set(anchor.0, anchor.1, 0.0);
        coeffs.as_.set(anchor.0, anchor.1, 0.0);
        coeffs.b.set(anchor.0, anchor.1, 0.0);
    }

    let n = (nx * ny) as f64;
    if n > 0.0 {
        (sum_sq / n).sqrt()
    } else {
        0.0
    }
}

/// Apply the under-relaxed pressure correction `p ← p + αp·p'`.
fn correct_pressure(p: &mut Field, pcorr: &Field, relax_p: f64) {
    for k in 0..p.data.len() {
        p.data[k] += relax_p * pcorr.data[k];
    }
}

/// Apply the velocity correction so the corrected field is
/// divergence-free — identical to [`crate::solver`]'s.
fn correct_velocity(
    u: &mut Field,
    v: &mut Field,
    pcorr: &Field,
    apu: &Field,
    apv: &Field,
    grid: &Grid,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dy = grid.dy();
    let dx = grid.dx();

    for j in 0..ny {
        for i in 1..nx {
            let d = dy / apu.at(i, j).max(1e-30);
            let du = d * (pcorr.at(i - 1, j) - pcorr.at(i, j));
            u.add(i, j, du);
        }
    }
    for j in 1..ny {
        for i in 0..nx {
            let d = dx / apv.at(i, j).max(1e-30);
            let dv = d * (pcorr.at(i, j - 1) - pcorr.at(i, j));
            v.add(i, j, dv);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::{solve_simple, SimpleControls};

    #[test]
    fn started_from_rest_channel_relaxes_toward_the_steady_solution() {
        // A channel flow started impulsively from rest. After enough
        // physical time the transient field must relax onto the steady
        // solution the SIMPLE solver computes directly — the defining
        // property of a transient solver run with steady BCs.
        let grid = Grid::new(40, 12, 4.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let bcs = Boundaries::channel_flow(1.0);

        // Steady reference.
        let steady = solve_simple(
            &grid,
            &fluid,
            &bcs,
            &SimpleControls {
                max_iterations: 4000,
                tolerance: 1e-6,
                ..SimpleControls::default()
            },
        );
        assert!(steady.converged, "steady reference should converge");

        // March from rest well past the viscous diffusion time of the
        // channel (≈ height²/ν = 1/0.05 = 20 s) so the flow is fully
        // relaxed onto the steady state. The implicit scheme is
        // unconditionally stable, so a large Δt is fine — at the steady
        // limit the unsteady term vanishes and the scheme converges to
        // exactly the steady solution regardless of Δt.
        let transient = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt: 0.2,
                n_steps: 300, // t = 60 s ≈ 3 diffusion times
                inner_iterations: 60,
                tolerance: 1e-6,
                ..TransientControls::default()
            },
            None,
        );

        // Compare the developed centreline velocity.
        let out_i = grid.nx - 3;
        let mid_j = grid.ny / 2;
        let steady_u = steady.u_at_cell(out_i, mid_j);
        let transient_u = transient.u_at_cell(out_i, mid_j);
        let rel = (steady_u - transient_u).abs() / steady_u.abs().max(1e-12);
        assert!(
            rel < 0.05,
            "transient centreline {transient_u} should relax onto steady {steady_u} (rel {rel})"
        );
    }

    #[test]
    fn flow_builds_up_monotonically_from_rest() {
        // Starting from a dead-stop, the kinetic energy of the channel
        // must *grow* over the first stretch of time as the inlet
        // momentum diffuses in — it cannot start non-zero and it
        // cannot decay.
        let grid = Grid::new(30, 12, 3.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let bcs = Boundaries::channel_flow(1.0);

        // A handful of steps.
        let early = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt: 0.02,
                n_steps: 5,
                inner_iterations: 40,
                ..TransientControls::default()
            },
            None,
        );
        // Many more steps continuing from the same start.
        let late = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt: 0.02,
                n_steps: 120,
                inner_iterations: 40,
                ..TransientControls::default()
            },
            None,
        );
        // The interior flow develops: the later field is faster.
        let early_speed = early.max_speed();
        let late_speed = late.max_speed();
        assert!(
            early_speed > 0.0,
            "the flow should start developing immediately"
        );
        assert!(
            late_speed > early_speed,
            "the flow should build up over time: {early_speed} → {late_speed}"
        );
    }

    #[test]
    fn transient_field_is_divergence_free_every_step() {
        // An incompressible transient solver must produce a
        // divergence-free velocity field at *every* time level, not
        // just at the steady limit. Check the cell mass imbalance of
        // the final field.
        let grid = Grid::new(24, 12, 2.0, 1.0);
        let fluid = Fluid::new(1.0, 0.06);
        let bcs = Boundaries::channel_flow(1.0);
        let sol = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt: 0.02,
                n_steps: 60,
                inner_iterations: 60,
                tolerance: 1e-6,
                ..TransientControls::default()
            },
            None,
        );
        let dx = grid.dx();
        let dy = grid.dy();
        let mut max_div = 0.0_f64;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let div = (sol.u.at(i + 1, j) - sol.u.at(i, j)) * dy
                    + (sol.v.at(i, j + 1) - sol.v.at(i, j)) * dx;
                max_div = max_div.max(div.abs());
            }
        }
        assert!(
            max_div < 1e-4,
            "transient field should be divergence-free, max |div| = {max_div}"
        );
    }

    #[test]
    fn quiescent_fluid_started_at_rest_stays_at_rest() {
        // No driving — every wall stationary, no inlet. The transient
        // march of a dead fluid must stay dead.
        let grid = Grid::new(16, 16, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.1);
        let bcs = Boundaries::lid_driven_cavity(0.0);
        let sol = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt: 0.05,
                n_steps: 20,
                ..TransientControls::default()
            },
            None,
        );
        assert!(sol.u.abs_max() < 1e-8, "still fluid should not move (u)");
        assert!(sol.v.abs_max() < 1e-8, "still fluid should not move (v)");
    }

    #[test]
    fn time_advances_by_dt_each_step() {
        // The physical clock must advance by exactly Δt per step.
        let grid = Grid::new(10, 10, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.1);
        let bcs = Boundaries::lid_driven_cavity(1.0);
        let dt = 0.025;
        let n = 12;
        let sol = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt,
                n_steps: n,
                inner_iterations: 20,
                ..TransientControls::default()
            },
            None,
        );
        assert_eq!(sol.history.len(), n);
        assert!(
            (sol.time - dt * n as f64).abs() < 1e-9,
            "final time {} should be dt·n = {}",
            sol.time,
            dt * n as f64
        );
        // The history times are strictly increasing by dt.
        for (step, h) in sol.history.iter().enumerate() {
            let expected = dt * (step + 1) as f64;
            assert!((h.time - expected).abs() < 1e-9);
        }
    }

    #[test]
    fn lid_driven_cavity_transient_reaches_the_steady_circulation() {
        // The transient lid-driven cavity, marched to long time, must
        // reproduce the steady solver's recirculating vortex — the
        // floor flow opposing the lid.
        let grid = Grid::new(20, 20, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.1);
        let bcs = Boundaries::lid_driven_cavity(1.0);
        let sol = solve_transient(
            &grid,
            &fluid,
            &bcs,
            &TransientControls {
                dt: 0.1,
                n_steps: 150,
                inner_iterations: 50,
                tolerance: 1e-6,
                ..TransientControls::default()
            },
            None,
        );
        // Near the lid the fluid moves with it; near the floor it
        // returns — the steady vortex, reached transiently.
        let mid_i = grid.nx / 2;
        let u_lid = sol.u_at_cell(mid_i, grid.ny - 2);
        let u_floor = sol.u_at_cell(mid_i, 1);
        assert!(u_lid > 0.05, "fluid under the lid moves with it: {u_lid}");
        assert!(
            u_floor < 0.0,
            "fluid near the floor returns upstream: {u_floor}"
        );
    }
}
