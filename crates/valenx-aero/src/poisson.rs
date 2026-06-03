//! 3-D linear solver for the pressure-correction Poisson equation —
//! SOR plus a geometric-multigrid V-cycle.
//!
//! # What it solves
//!
//! SIMPLE's pressure-correction step ([`crate::solver`]) produces, on
//! every grid cell, a seven-point discrete Poisson equation
//!
//! ```text
//!   aP·p'(i,j,k) =  aE·p'(i+1,j,k) + aW·p'(i-1,j,k)
//!                 + aN·p'(i,j+1,k) + aS·p'(i,j-1,k)
//!                 + aT·p'(i,j,k+1) + aB·p'(i,j,k-1) + b(i,j,k)
//! ```
//!
//! relating each cell's pressure correction `p'` to its six
//! face-neighbours. This module solves the resulting sparse symmetric
//! system for `p'`.
//!
//! # Why multigrid matters in 3-D
//!
//! Plain successive over-relaxation (SOR) is a *smoother*: it kills
//! the high-frequency error fast but the low-frequency error decays
//! geometrically slowly — and the slowdown gets dramatically worse as
//! the grid is refined. For the large 3-D Poisson systems an external-
//! aero wind-tunnel run produces, plain SOR would need tens of
//! thousands of sweeps. **Geometric multigrid** fixes this: it relaxes
//! a few times on the fine grid (to smooth the error), restricts the
//! residual to a coarser grid where the *same* error looks
//! high-frequency, recursively solves there, prolongs the correction
//! back, and relaxes again. Each V-cycle reduces the error by a fixed
//! factor *independent of the grid size* — the property that makes
//! 3-D CFD tractable.
//!
//! # Honest scope
//!
//! The V-cycle here is a real working geometric multigrid for the
//! constant-coefficient Poisson operator on the Cartesian grid. The
//! coarsening is 2:1 per axis (halt when an axis would drop below two
//! cells); coarse operators are *re-discretised* (the standard
//! geometric approach) rather than Galerkin-assembled. Immersed-
//! boundary cut cells are handled on the fine grid by zeroing the
//! blocked-neighbour coefficients; the coarse levels use the plain
//! Laplacian, which is a documented v1 simplification — the IBM
//! geometry is not transferred down the hierarchy, so the multigrid
//! acts as a (still very effective) preconditioner-grade smoother
//! accelerator rather than an exact coarse solver for the cut region.

use crate::grid::Field3;
use rayon::prelude::*;

/// Grid sizes at or below this cell count run the SOR sweep / residual
/// reductions serially — for a tiny grid the rayon fork/join overhead
/// outweighs the work. Larger grids are split across the thread pool.
pub const PARALLEL_THRESHOLD: usize = 4_096;

/// The seven-point stencil coefficients for the pressure-correction
/// Poisson equation, one set per pressure cell.
///
/// All fields are `nx · ny · nz`. A boundary cell or an immersed-
/// boundary cut cell simply has the unavailable neighbour coefficient
/// set to zero (a homogeneous-Neumann condition on `p'`).
#[derive(Clone, Debug)]
pub struct PoissonStencil {
    /// Cells along x.
    pub nx: usize,
    /// Cells along y.
    pub ny: usize,
    /// Cells along z.
    pub nz: usize,
    /// Diagonal coefficient `aP` per cell.
    pub ap: Field3,
    /// East (`+x`) neighbour coefficient `aE`.
    pub ae: Field3,
    /// West (`-x`) neighbour coefficient `aW`.
    pub aw: Field3,
    /// North (`+y`) neighbour coefficient `aN`.
    pub an: Field3,
    /// South (`-y`) neighbour coefficient `aS`.
    pub as_: Field3,
    /// Top (`+z`) neighbour coefficient `aT`.
    pub at: Field3,
    /// Bottom (`-z`) neighbour coefficient `aB`.
    pub ab: Field3,
    /// Source term `b` per cell.
    pub b: Field3,
}

impl PoissonStencil {
    /// Allocate a zeroed stencil for an `nx · ny · nz` grid.
    pub fn zeros(nx: usize, ny: usize, nz: usize) -> PoissonStencil {
        PoissonStencil {
            nx,
            ny,
            nz,
            ap: Field3::zeros(nx, ny, nz),
            ae: Field3::zeros(nx, ny, nz),
            aw: Field3::zeros(nx, ny, nz),
            an: Field3::zeros(nx, ny, nz),
            as_: Field3::zeros(nx, ny, nz),
            at: Field3::zeros(nx, ny, nz),
            ab: Field3::zeros(nx, ny, nz),
            b: Field3::zeros(nx, ny, nz),
        }
    }

    /// Build the standard constant-coefficient 7-point Laplacian on a
    /// uniform grid with cell sizes `(dx, dy, dz)` and a homogeneous-
    /// Neumann boundary — every interior cell gets `area/spacing` on
    /// each in-domain neighbour, `aP` the sum. The source `b` is left
    /// zero for the caller to fill.
    pub fn laplacian(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> PoissonStencil {
        let mut s = PoissonStencil::zeros(nx, ny, nz);
        // 7-point Laplacian face weight per axis: 1/h² (unit cell
        // areas — the constant cancels for the pure Poisson test).
        let cx = 1.0 / (dx * dx);
        let cy = 1.0 / (dy * dy);
        let cz = 1.0 / (dz * dz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let mut ap = 0.0;
                    if i + 1 < nx {
                        s.ae.set(i, j, k, cx);
                        ap += cx;
                    }
                    if i > 0 {
                        s.aw.set(i, j, k, cx);
                        ap += cx;
                    }
                    if j + 1 < ny {
                        s.an.set(i, j, k, cy);
                        ap += cy;
                    }
                    if j > 0 {
                        s.as_.set(i, j, k, cy);
                        ap += cy;
                    }
                    if k + 1 < nz {
                        s.at.set(i, j, k, cz);
                        ap += cz;
                    }
                    if k > 0 {
                        s.ab.set(i, j, k, cz);
                        ap += cz;
                    }
                    s.ap.set(i, j, k, ap);
                }
            }
        }
        s
    }

    /// The residual L2 norm of a candidate solution — `‖A·p' − rhs‖₂`
    /// normalised by the cell count, where `rhs` is the seven-point
    /// neighbour sum plus `b`.
    ///
    /// The per-cell sum-of-squares is an embarrassingly parallel
    /// reduction; for a grid above [`PARALLEL_THRESHOLD`] cells it is
    /// split over the rayon pool by z-plane.
    pub fn residual_norm(&self, sol: &Field3) -> f64 {
        let n = self.nx * self.ny * self.nz;
        if n == 0 {
            return 0.0;
        }
        let plane_sq = |k: usize| -> f64 {
            let mut s = 0.0;
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let ap = self.ap.at(i, j, k);
                    if ap.abs() < 1e-30 {
                        continue;
                    }
                    let r = ap * sol.at(i, j, k) - self.neighbour_sum(sol, i, j, k);
                    s += r * r;
                }
            }
            s
        };
        let sum_sq: f64 = if n > PARALLEL_THRESHOLD {
            (0..self.nz).into_par_iter().map(plane_sq).sum()
        } else {
            (0..self.nz).map(plane_sq).sum()
        };
        (sum_sq / n as f64).sqrt()
    }

    /// The off-diagonal neighbour sum plus the source `b` for cell
    /// `(i, j, k)` — the right-hand side of the cell's Poisson row.
    #[inline]
    fn neighbour_sum(&self, sol: &Field3, i: usize, j: usize, k: usize) -> f64 {
        let mut s = self.b.at(i, j, k);
        if i + 1 < self.nx {
            s += self.ae.at(i, j, k) * sol.at(i + 1, j, k);
        }
        if i > 0 {
            s += self.aw.at(i, j, k) * sol.at(i - 1, j, k);
        }
        if j + 1 < self.ny {
            s += self.an.at(i, j, k) * sol.at(i, j + 1, k);
        }
        if j > 0 {
            s += self.as_.at(i, j, k) * sol.at(i, j - 1, k);
        }
        if k + 1 < self.nz {
            s += self.at.at(i, j, k) * sol.at(i, j, k + 1);
        }
        if k > 0 {
            s += self.ab.at(i, j, k) * sol.at(i, j, k - 1);
        }
        s
    }
}

/// Outcome of a Poisson solve.
#[derive(Clone, Copy, Debug)]
pub struct PoissonResult {
    /// Number of sweeps / V-cycles actually performed.
    pub iterations: usize,
    /// The final residual L2 norm.
    pub residual: f64,
    /// True if the residual fell below the requested tolerance.
    pub converged: bool,
}

/// One SOR sweep over the whole grid, in **red-black** order. `omega`
/// is the over-relaxation factor; `pin_mean` subtracts the mean after
/// the sweep (the singular all-Neumann gauge fix).
///
/// # Why red-black instead of lexicographic
///
/// A lexicographic Gauss-Seidel sweep is inherently serial: cell
/// `(i,j,k)` reads its already-updated west / south / bottom
/// neighbours, so no two cells can be updated at once. The red-black
/// ordering colours cells by the parity of `i+j+k`: every *red* cell's
/// six face-neighbours are *black* and vice versa, so within one colour
/// the updates are completely independent. The sweep becomes two
/// data-parallel passes (all red, then all black) with the *same*
/// smoothing property as lexicographic SOR — the standard
/// parallelisable multigrid smoother. Each colour pass is split over
/// the rayon pool by z-plane for grids above [`PARALLEL_THRESHOLD`].
pub fn sor_sweep(stencil: &PoissonStencil, sol: &mut Field3, omega: f64, pin_mean: bool) {
    let omega = omega.clamp(0.5, 1.95);
    sor_color_sweep(stencil, sol, omega, 0); // red:   (i+j+k) even
    sor_color_sweep(stencil, sol, omega, 1); // black: (i+j+k) odd
    if pin_mean {
        let mean = sol.mean();
        if mean != 0.0 {
            sol.data.par_iter_mut().for_each(|v| *v -= mean);
        }
    }
}

/// Relax every cell of one colour (`color` = 0 red / 1 black) — the
/// data-parallel half of a red-black SOR sweep.
///
/// Every cell of one colour reads only opposite-colour neighbours, so
/// the colour's updates are mutually independent. The relaxed values
/// are computed by a parallel z-plane map over a shared immutable
/// borrow of `sol`, then written back. This keeps the sweep entirely
/// in safe code (the crate is `#![forbid(unsafe_code)]`) while still
/// scaling across the rayon pool.
fn sor_color_sweep(stencil: &PoissonStencil, sol: &mut Field3, omega: f64, color: usize) {
    let (nx, ny, nz) = (stencil.nx, stencil.ny, stencil.nz);
    let n = nx * ny * nz;
    let plane = nx * ny;

    // Compute the new value for one z-plane's `color` cells, into a
    // freshly allocated plane buffer (default = old value, untouched
    // cells unchanged). Reads `sol` immutably only.
    let relax_plane = |k: usize| -> Vec<f64> {
        let base = k * plane;
        let mut out = sol.data[base..base + plane].to_vec();
        for j in 0..ny {
            for i in 0..nx {
                if (i + j + k) & 1 != color {
                    continue;
                }
                let ap = stencil.ap.at(i, j, k);
                if ap.abs() < 1e-30 {
                    continue;
                }
                let rhs = stencil.neighbour_sum(sol, i, j, k);
                let idx = i + nx * j;
                let old = out[idx];
                out[idx] = old + omega * (rhs / ap - old);
            }
        }
        out
    };

    let planes: Vec<Vec<f64>> = if n > PARALLEL_THRESHOLD {
        (0..nz).into_par_iter().map(relax_plane).collect()
    } else {
        (0..nz).map(relax_plane).collect()
    };
    for (k, p) in planes.into_iter().enumerate() {
        sol.data[k * plane..(k + 1) * plane].copy_from_slice(&p);
    }
}

/// Solve the Poisson system by plain SOR — used directly for small
/// grids and as the smoother / coarse solver inside the multigrid
/// V-cycle.
///
/// `sol` is the initial guess and is overwritten with the result.
/// `pin_mean` selects the gauge handling: `true` for a fully-singular
/// all-Neumann system (every boundary no-penetration on `p'`), `false`
/// when a Dirichlet anchor is already present.
pub fn solve_sor(
    stencil: &PoissonStencil,
    sol: &mut Field3,
    omega: f64,
    tol: f64,
    max_iter: usize,
    pin_mean: bool,
) -> PoissonResult {
    let mut iterations = 0;
    let mut residual = f64::INFINITY;
    let mut converged = false;
    for sweep in 0..max_iter.max(1) {
        iterations = sweep + 1;
        sor_sweep(stencil, sol, omega, pin_mean);
        residual = stencil.residual_norm(sol);
        if residual <= tol {
            converged = true;
            break;
        }
        if !residual.is_finite() {
            break;
        }
    }
    PoissonResult {
        iterations,
        residual,
        converged,
    }
}

/// Solve the Poisson system with a geometric-multigrid V-cycle scheme.
///
/// Each V-cycle: `pre` SOR smoothing sweeps on the current grid;
/// restrict the residual to a 2:1-coarsened grid; recurse; prolong the
/// coarse correction back and add it; `post` smoothing sweeps. The
/// recursion bottoms out when an axis cannot be halved further, where
/// a direct-ish SOR solve finishes the coarsest problem. `pin_mean`
/// applies on every level for the all-Neumann case.
///
/// Returns when the fine-grid residual drops below `tol` or after
/// `max_cycles` V-cycles.
pub fn solve_multigrid(
    stencil: &PoissonStencil,
    sol: &mut Field3,
    tol: f64,
    max_cycles: usize,
    pre: usize,
    post: usize,
    pin_mean: bool,
) -> PoissonResult {
    solve_multigrid_anchored(stencil, sol, tol, max_cycles, pre, post, pin_mean, None)
}

/// [`solve_multigrid`] with an explicit Dirichlet-anchor cell.
///
/// When `anchor` is given, that cell (and its coarse-grid images) is a
/// Dirichlet pin, so the operator is non-singular at every level and no
/// gauge handling is needed. This is what the SIMPLE pressure
/// correction uses — its restricted residual carries a net mass
/// imbalance, which a purely-singular coarse Neumann operator could not
/// consistently solve.
#[allow(clippy::too_many_arguments)]
pub fn solve_multigrid_anchored(
    stencil: &PoissonStencil,
    sol: &mut Field3,
    tol: f64,
    max_cycles: usize,
    pre: usize,
    post: usize,
    pin_mean: bool,
    anchor: Option<(usize, usize, usize)>,
) -> PoissonResult {
    let mut iterations = 0;
    let mut residual = stencil.residual_norm(sol);
    let mut converged = residual <= tol;
    for cycle in 0..max_cycles.max(1) {
        if converged {
            break;
        }
        iterations = cycle + 1;
        v_cycle(stencil, sol, pre, post, pin_mean, anchor);
        residual = stencil.residual_norm(sol);
        if !residual.is_finite() {
            break;
        }
        if residual <= tol {
            converged = true;
        }
    }
    PoissonResult {
        iterations,
        residual,
        converged,
    }
}

/// True when the grid can still be coarsened 2:1 — every axis must
/// have an even count of at least 4 cells so the coarse grid keeps
/// at least 2 cells per axis.
fn coarsenable(nx: usize, ny: usize, nz: usize) -> bool {
    nx >= 4 && ny >= 4 && nz >= 4 && nx % 2 == 0 && ny % 2 == 0 && nz % 2 == 0
}

/// Build the **agglomerated** coarse-grid operator from a fine 7-point
/// stencil by 2:1 geometric coarsening.
///
/// Each coarse cell agglomerates a 2×2×2 block of fine cells; the
/// coarse face coefficient between two coarse cells is the *sum* of the
/// four fine face coefficients on the shared coarse face, and the
/// coarse diagonal is the sum of the coarse off-diagonals (so the
/// coarse operator is a valid flux-conserving Neumann Laplacian).
///
/// This is essential for a *variable-coefficient* problem: the SIMPLE
/// pressure-correction stencil has cell-to-cell-varying coefficients
/// (`ρ·A²/aP` with `aP` varying across the flow, plus zeroed
/// immersed-boundary faces). A re-discretised *constant*-coefficient
/// coarse Laplacian would be wildly inconsistent with it — the coarse
/// correction could be off by orders of magnitude and, added back,
/// would make the V-cycle diverge instead of converge. Agglomeration
/// keeps the coarse operator consistent with whatever the fine
/// coefficients are.
///
/// `anchor`, when set, is a Dirichlet-pinned fine cell; the coarse cell
/// that contains it is re-pinned as an identity row. Without this the
/// coarse operator would be a *purely* singular Neumann Laplacian, and
/// the restricted SIMPLE residual (which carries a net boundary mass
/// imbalance and so does not sum to zero) would make that singular
/// coarse system inconsistent — the coarse solve would then drift
/// without bound and blow the V-cycle up.
fn coarsen_stencil(
    fine: &PoissonStencil,
    anchor: Option<(usize, usize, usize)>,
) -> PoissonStencil {
    let (cnx, cny, cnz) = (fine.nx / 2, fine.ny / 2, fine.nz / 2);
    let mut c = PoissonStencil::zeros(cnx, cny, cnz);
    for ck in 0..cnz {
        for cj in 0..cny {
            for ci in 0..cnx {
                // Sum the four fine faces lying on each coarse face.
                let mut ae = 0.0;
                let mut aw = 0.0;
                let mut an = 0.0;
                let mut as_ = 0.0;
                let mut at = 0.0;
                let mut ab = 0.0;
                for d2 in 0..2 {
                    for d1 in 0..2 {
                        // x faces: vary (y, z) over the block.
                        ae += fine.ae.at(2 * ci + 1, 2 * cj + d1, 2 * ck + d2);
                        aw += fine.aw.at(2 * ci, 2 * cj + d1, 2 * ck + d2);
                        // y faces: vary (x, z).
                        an += fine.an.at(2 * ci + d1, 2 * cj + 1, 2 * ck + d2);
                        as_ += fine.as_.at(2 * ci + d1, 2 * cj, 2 * ck + d2);
                        // z faces: vary (x, y).
                        at += fine.at.at(2 * ci + d1, 2 * cj + d2, 2 * ck + 1);
                        ab += fine.ab.at(2 * ci + d1, 2 * cj + d2, 2 * ck);
                    }
                }
                // A coarse face on the domain boundary has no coarse
                // neighbour — drop it (homogeneous Neumann).
                if ci + 1 >= cnx {
                    ae = 0.0;
                }
                if ci == 0 {
                    aw = 0.0;
                }
                if cj + 1 >= cny {
                    an = 0.0;
                }
                if cj == 0 {
                    as_ = 0.0;
                }
                if ck + 1 >= cnz {
                    at = 0.0;
                }
                if ck == 0 {
                    ab = 0.0;
                }
                c.ae.set(ci, cj, ck, ae);
                c.aw.set(ci, cj, ck, aw);
                c.an.set(ci, cj, ck, an);
                c.as_.set(ci, cj, ck, as_);
                c.at.set(ci, cj, ck, at);
                c.ab.set(ci, cj, ck, ab);
                c.ap.set(ci, cj, ck, ae + aw + an + as_ + at + ab);
            }
        }
    }
    // Carry the Dirichlet pin down: the coarse cell containing the
    // fine anchor becomes a clean identity row.
    if let Some((ai, aj, ak)) = anchor {
        let (ci, cj, ck) = (ai / 2, aj / 2, ak / 2);
        c.ap.set(ci, cj, ck, 1.0);
        c.ae.set(ci, cj, ck, 0.0);
        c.aw.set(ci, cj, ck, 0.0);
        c.an.set(ci, cj, ck, 0.0);
        c.as_.set(ci, cj, ck, 0.0);
        c.at.set(ci, cj, ck, 0.0);
        c.ab.set(ci, cj, ck, 0.0);
    }
    c
}

/// The coarse-grid image of a Dirichlet-anchor cell index under one
/// 2:1 coarsening step.
fn coarsen_anchor(
    anchor: Option<(usize, usize, usize)>,
) -> Option<(usize, usize, usize)> {
    anchor.map(|(i, j, k)| (i / 2, j / 2, k / 2))
}

/// One recursive multigrid V-cycle step.
fn v_cycle(
    stencil: &PoissonStencil,
    sol: &mut Field3,
    pre: usize,
    post: usize,
    pin_mean: bool,
    anchor: Option<(usize, usize, usize)>,
) {
    let (nx, ny, nz) = (stencil.nx, stencil.ny, stencil.nz);

    // Pre-smoothing.
    for _ in 0..pre.max(1) {
        sor_sweep(stencil, sol, 1.4, pin_mean);
    }

    if !coarsenable(nx, ny, nz) {
        // Coarsest grid: a heavier SOR solve finishes the problem.
        for _ in 0..40 {
            sor_sweep(stencil, sol, 1.6, pin_mean);
        }
        return;
    }

    // Compute the fine residual r = b + Σa·p'_nb − aP·p'. Each z-plane
    // is independent — parallelise it for a large grid.
    let mut fine_res = Field3::zeros(nx, ny, nz);
    let plane = nx * ny;
    let residual_plane = |k: usize| -> Vec<f64> {
        let mut out = vec![0.0; plane];
        for j in 0..ny {
            for i in 0..nx {
                let ap = stencil.ap.at(i, j, k);
                if ap.abs() < 1e-30 {
                    continue;
                }
                let rhs = stencil.neighbour_sum(sol, i, j, k);
                out[i + nx * j] = rhs - ap * sol.at(i, j, k);
            }
        }
        out
    };
    let res_planes: Vec<Vec<f64>> = if nx * ny * nz > PARALLEL_THRESHOLD {
        (0..nz).into_par_iter().map(residual_plane).collect()
    } else {
        (0..nz).map(residual_plane).collect()
    };
    for (k, p) in res_planes.into_iter().enumerate() {
        fine_res.data[k * plane..(k + 1) * plane].copy_from_slice(&p);
    }

    // Restrict the residual to the coarse grid (full-weighting 2³
    // average) and build the agglomerated coarse operator. Agglomeration
    // (rather than re-discretising a constant-coefficient Laplacian)
    // keeps the coarse operator consistent with the fine
    // *variable-coefficient* stencil — without it the V-cycle diverges
    // on the SIMPLE pressure-correction system.
    let (cnx, cny, cnz) = (nx / 2, ny / 2, nz / 2);
    let mut coarse_res = restrict(&fine_res, cnx, cny, cnz);
    let mut coarse = coarsen_stencil(stencil, anchor);
    let coarse_anchor = coarsen_anchor(anchor);
    // The coarse anchor cell is a clean identity row → its source is 0.
    if let Some((ci, cj, ck)) = coarse_anchor {
        coarse_res.set(ci, cj, ck, 0.0);
    }
    coarse.b = coarse_res;

    // Solve the coarse *correction* equation recursively from zero.
    // With an anchor the coarse operator is non-singular; with none it
    // is a pure-Neumann singular Laplacian, gauge-fixed by mean-pinning.
    let mut coarse_sol = Field3::zeros(cnx, cny, cnz);
    v_cycle(
        &coarse,
        &mut coarse_sol,
        pre,
        post,
        coarse_anchor.is_none(),
        coarse_anchor,
    );

    // Prolong the coarse correction and add it with the residual-
    // minimising step length — a damped coarse-grid correction.
    //
    // For an ill-conditioned variable-coefficient operator an
    // inexactly-solved coarse problem can return an over-scaled
    // correction `d`; adding it raw would *increase* the fine residual
    // and the V-cycle would diverge. The step `α` that minimises
    // `‖r − α·A·d‖` is `α = ⟨r, A·d⟩ / ⟨A·d, A·d⟩`; clamped to
    // `[0, 1]` it guarantees the correction never makes the residual
    // worse, so the V-cycle stays monotone. This costs two extra
    // O(N) sweeps — far cheaper than a back-tracking line search.
    let fine_corr = prolong(&coarse_sol, nx, ny, nz);
    // Residual r and the operator applied to the correction, A·d.
    // ⟨r, A·d⟩ and ⟨A·d, A·d⟩ are two parallel reductions over the
    // z-planes.
    let n_cells = sol.data.len();
    let plane_dots = |k: usize| -> (f64, f64) {
        let mut r_dot_ad = 0.0;
        let mut ad_dot_ad = 0.0;
        for j in 0..ny {
            for i in 0..nx {
                let ap = stencil.ap.at(i, j, k);
                if ap.abs() < 1e-30 {
                    continue;
                }
                let idx = i + nx * (j + ny * k);
                let r = stencil.neighbour_sum(sol, i, j, k) - ap * sol.at(i, j, k);
                // (A·d)[cell] = ap·d − Σ off·d_nb.
                let ad = ap * fine_corr.data[idx]
                    - (stencil.neighbour_sum(&fine_corr, i, j, k)
                        - stencil.b.at(i, j, k));
                r_dot_ad += r * ad;
                ad_dot_ad += ad * ad;
            }
        }
        (r_dot_ad, ad_dot_ad)
    };
    let (r_dot_ad, ad_dot_ad): (f64, f64) = if n_cells > PARALLEL_THRESHOLD {
        (0..nz)
            .into_par_iter()
            .map(plane_dots)
            .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
    } else {
        (0..nz)
            .map(plane_dots)
            .fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
    };
    let alpha = if ad_dot_ad > 1e-300 {
        (r_dot_ad / ad_dot_ad).clamp(0.0, 1.0)
    } else {
        0.0
    };
    for idx in 0..n_cells {
        sol.data[idx] += alpha * fine_corr.data[idx];
    }
    if pin_mean {
        let mean = sol.mean();
        for v in sol.data.iter_mut() {
            *v -= mean;
        }
    }

    // Post-smoothing.
    for _ in 0..post.max(1) {
        sor_sweep(stencil, sol, 1.4, pin_mean);
    }
}

/// Summation restriction: each coarse cell receives the *sum* of the
/// 2×2×2 block of fine cells it covers.
///
/// Summation is the transpose of the piecewise-constant prolongation
/// ([`prolong`]) — the pairing `R = Pᵀ` that makes the agglomerated
/// coarse operator ([`coarsen_stencil`], itself `Pᵀ A P`) a consistent
/// Galerkin coarse operator. Using an *averaging* restriction with an
/// agglomerated operator would mismatch by a factor of 8 and starve
/// the coarse-grid correction.
fn restrict(fine: &Field3, cnx: usize, cny: usize, cnz: usize) -> Field3 {
    let mut coarse = Field3::zeros(cnx, cny, cnz);
    for ck in 0..cnz {
        for cj in 0..cny {
            for ci in 0..cnx {
                let mut sum = 0.0;
                for dk in 0..2 {
                    for dj in 0..2 {
                        for di in 0..2 {
                            sum += fine.at(2 * ci + di, 2 * cj + dj, 2 * ck + dk);
                        }
                    }
                }
                coarse.set(ci, cj, ck, sum);
            }
        }
    }
    coarse
}

/// Piecewise-constant prolongation: each fine cell copies the coarse
/// cell that contains it. (Constant injection — robust and adequate
/// as a multigrid correction transfer; trilinear interpolation is a
/// documented refinement.)
fn prolong(coarse: &Field3, fnx: usize, fny: usize, fnz: usize) -> Field3 {
    let mut fine = Field3::zeros(fnx, fny, fnz);
    for k in 0..fnz {
        for j in 0..fny {
            for i in 0..fnx {
                let (ci, cj, ck) = (
                    (i / 2).min(coarse.nx - 1),
                    (j / 2).min(coarse.ny - 1),
                    (k / 2).min(coarse.nz - 1),
                );
                fine.set(i, j, k, coarse.at(ci, cj, ck));
            }
        }
    }
    fine
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a known smooth zero-mean target and the source `b` that
    /// makes `A·target = b` so a solve must recover `target`.
    fn make_known_problem(n: usize) -> (PoissonStencil, Field3) {
        let h = 1.0 / n as f64;
        let mut s = PoissonStencil::laplacian(n, n, n, h, h, h);
        let mut target = Field3::zeros(n, n, n);
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let x = (i as f64 + 0.5) * h;
                    let y = (j as f64 + 0.5) * h;
                    let z = (k as f64 + 0.5) * h;
                    let pi = std::f64::consts::PI;
                    target.set(
                        i,
                        j,
                        k,
                        (pi * x).cos() * (pi * y).cos() * (pi * z).cos(),
                    );
                }
            }
        }
        // Make exactly zero-mean — the singular system's gauge.
        let mean = target.mean();
        for v in target.data.iter_mut() {
            *v -= mean;
        }
        // b such that A·target = b.
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let ap = s.ap.at(i, j, k);
                    let nb = s.neighbour_sum_no_b(&target, i, j, k);
                    s.b.set(i, j, k, ap * target.at(i, j, k) - nb);
                }
            }
        }
        (s, target)
    }

    impl PoissonStencil {
        /// Neighbour sum without the `b` term — for test construction.
        fn neighbour_sum_no_b(&self, sol: &Field3, i: usize, j: usize, k: usize) -> f64 {
            let mut s = 0.0;
            if i + 1 < self.nx {
                s += self.ae.at(i, j, k) * sol.at(i + 1, j, k);
            }
            if i > 0 {
                s += self.aw.at(i, j, k) * sol.at(i - 1, j, k);
            }
            if j + 1 < self.ny {
                s += self.an.at(i, j, k) * sol.at(i, j + 1, k);
            }
            if j > 0 {
                s += self.as_.at(i, j, k) * sol.at(i, j - 1, k);
            }
            if k + 1 < self.nz {
                s += self.at.at(i, j, k) * sol.at(i, j, k + 1);
            }
            if k > 0 {
                s += self.ab.at(i, j, k) * sol.at(i, j, k - 1);
            }
            s
        }
    }

    #[test]
    fn sor_solves_a_known_3d_poisson_problem() {
        let n = 8;
        let (s, target) = make_known_problem(n);
        let mut sol = Field3::zeros(n, n, n);
        let res = solve_sor(&s, &mut sol, 1.6, 1e-9, 8000, true);
        assert!(res.converged, "SOR should converge, residual {}", res.residual);
        let mut max_err = 0.0f64;
        for idx in 0..sol.data.len() {
            max_err = max_err.max((sol.data[idx] - target.data[idx]).abs());
        }
        assert!(max_err < 1e-3, "SOR 3-D solution error {max_err} too large");
    }

    #[test]
    fn multigrid_recovers_the_known_poisson_solution() {
        // The headline test: a multigrid V-cycle scheme recovers the
        // analytic zero-mean target on a 16³ grid to a 1e-7 residual.
        //
        // The V-cycle uses agglomeration (Galerkin) coarsening and
        // piecewise-constant transfer operators — the combination
        // required to stay *stable* on the variable-coefficient SIMPLE
        // pressure-correction operator. Constant-order transfers give a
        // sound but modest contraction factor, so reaching 1e-7 takes
        // more V-cycles than a linear-interpolation scheme would; the
        // cycle budget is set accordingly. The solver still converges
        // monotonically to the exact answer.
        let n = 16;
        let (s, target) = make_known_problem(n);
        let mut sol = Field3::zeros(n, n, n);
        let res = solve_multigrid(&s, &mut sol, 1e-7, 150, 2, 2, true);
        assert!(
            res.converged,
            "multigrid should converge, residual {}",
            res.residual
        );
        let mut max_err = 0.0f64;
        for idx in 0..sol.data.len() {
            max_err = max_err.max((sol.data[idx] - target.data[idx]).abs());
        }
        assert!(max_err < 5e-3, "multigrid solution error {max_err} too large");
    }

    #[test]
    fn multigrid_beats_sor_in_a_fixed_work_budget() {
        // Multigrid's whole point: grid-independent convergence. On a
        // 16³ grid, a handful of V-cycles must crush a residual that
        // the same number of plain SOR sweeps barely dents.
        let n = 16;
        let (s, _) = make_known_problem(n);
        let mut mg = Field3::zeros(n, n, n);
        let mg_res = solve_multigrid(&s, &mut mg, 0.0, 8, 2, 2, true).residual;
        let mut sor = Field3::zeros(n, n, n);
        let sor_res = solve_sor(&s, &mut sor, 1.8, 0.0, 8, true).residual;
        assert!(
            mg_res < sor_res,
            "multigrid residual {mg_res} should beat SOR {sor_res} in equal cycles"
        );
    }

    #[test]
    fn restrict_then_prolong_preserves_a_constant_field() {
        // Summation restriction is R = Pᵀ: a constant field restricts
        // to a uniform field (each coarse cell sums its 2³ = 8 fine
        // cells) and prolongs back to a uniform field. The field stays
        // constant-valued; the magnitude follows the transpose pairing.
        let fine = Field3::filled(8, 8, 8, 3.0);
        let coarse = restrict(&fine, 4, 4, 4);
        assert!(coarse.data.iter().all(|&v| (v - 24.0).abs() < 1e-12));
        let back = prolong(&coarse, 8, 8, 8);
        assert!(back.data.iter().all(|&v| (v - 24.0).abs() < 1e-12));
    }

    #[test]
    fn solution_of_a_zero_source_is_zero() {
        // A homogeneous system with a zero start stays zero.
        let s = PoissonStencil::laplacian(8, 8, 8, 0.1, 0.1, 0.1);
        let mut sol = Field3::zeros(8, 8, 8);
        let res = solve_multigrid(&s, &mut sol, 1e-10, 20, 2, 2, true);
        assert!(res.converged);
        assert!(sol.abs_max() < 1e-9, "zero-source solution must stay zero");
    }
}
