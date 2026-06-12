//! Iterative linear solver for the pressure-correction Poisson
//! equation.
//!
//! # What it solves
//!
//! The SIMPLE algorithm's pressure-correction step ([`crate::solver`])
//! produces, on every grid cell, a five-point discrete Poisson
//! equation
//!
//! ```text
//!   aP·p'(i,j) = aE·p'(i+1,j) + aW·p'(i-1,j)
//!              + aN·p'(i,j+1) + aS·p'(i,j-1) + b(i,j)
//! ```
//!
//! relating each cell's pressure correction `p'` to its four
//! neighbours. The coefficients `aE … aS` and the source `b` come from
//! the momentum equation's discretisation; this module solves the
//! resulting sparse symmetric system for `p'`.
//!
//! # The method
//!
//! **Successive over-relaxation (SOR)** — a Gauss-Seidel sweep with an
//! over-relaxation factor `ω`. Each sweep visits every cell, computes
//! the Gauss-Seidel update (using the freshest available neighbour
//! values), and over-shoots toward it by `ω` (`1 < ω < 2`), which
//! roughly squares the convergence rate of plain Gauss-Seidel. A direct
//! sparse factorisation would also work, but SOR needs no factorisation
//! storage, is trivially correct, and the pressure-correction system
//! only has to be solved *approximately* each SIMPLE outer iteration
//! anyway — SIMPLE is itself an outer iteration.
//!
//! The Poisson system from an all-Neumann (closed-wall) pressure-
//! correction problem is singular up to an additive constant; the
//! solver pins the correction by subtracting its mean after each sweep,
//! which selects the unique zero-mean solution.

use crate::grid::Field;

/// The five-point stencil coefficients for the pressure-correction
/// Poisson equation, one set per pressure cell.
///
/// All fields are `nx × ny`. A boundary cell simply has the
/// out-of-domain neighbour coefficient set to zero (a homogeneous
/// Neumann condition on `p'`), which the assembler in [`crate::solver`]
/// handles.
#[derive(Clone, Debug)]
pub struct PoissonCoeffs {
    /// Grid width (pressure cells along x).
    pub nx: usize,
    /// Grid height (pressure cells along y).
    pub ny: usize,
    /// Diagonal coefficient `aP` per cell.
    pub ap: Field,
    /// East-neighbour coefficient `aE` per cell.
    pub ae: Field,
    /// West-neighbour coefficient `aW` per cell.
    pub aw: Field,
    /// North-neighbour coefficient `aN` per cell.
    pub an: Field,
    /// South-neighbour coefficient `aS` per cell.
    pub as_: Field,
    /// Source term `b` per cell — the cell's mass imbalance.
    pub b: Field,
}

impl PoissonCoeffs {
    /// Allocate a zeroed coefficient set for an `nx × ny` grid.
    pub fn zeros(nx: usize, ny: usize) -> PoissonCoeffs {
        PoissonCoeffs {
            nx,
            ny,
            ap: Field::zeros(nx, ny),
            ae: Field::zeros(nx, ny),
            aw: Field::zeros(nx, ny),
            an: Field::zeros(nx, ny),
            as_: Field::zeros(nx, ny),
            b: Field::zeros(nx, ny),
        }
    }
}

/// Outcome of an SOR solve.
#[derive(Clone, Copy, Debug)]
pub struct SorResult {
    /// Number of sweeps actually performed.
    pub iterations: usize,
    /// The final residual L2 norm (`‖A·p' − b‖`, normalised by the
    /// cell count).
    pub residual: f64,
    /// True if the residual fell below the requested tolerance before
    /// the iteration cap.
    pub converged: bool,
}

/// Solve the pressure-correction Poisson system into `solution` by SOR.
///
/// `solution` is used as the initial guess and overwritten with the
/// result — passing the previous iteration's correction warm-starts the
/// solve. `omega` is the over-relaxation factor (`1.0` = plain
/// Gauss-Seidel; `1.5–1.9` is the usual productive range). The sweep
/// stops when the residual L2 norm drops below `tol` or after
/// `max_iter` sweeps.
///
/// `pin_mean` selects the gauge handling:
///
/// - `true` — the system is the **fully singular** all-Neumann case (a
///   closed domain: every boundary has a no-penetration condition on
///   `p'`). Its solution is defined only up to an additive constant, so
///   after every sweep the mean is subtracted and the unique zero-mean
///   correction is returned.
/// - `false` — the system is **non-singular** because a Dirichlet
///   condition is already present (an outlet cell pinned to `p' = 0`).
///   No mean subtraction is applied — that would fight the Dirichlet
///   anchor.
pub fn solve_sor(
    coeffs: &PoissonCoeffs,
    solution: &mut Field,
    omega: f64,
    tol: f64,
    max_iter: usize,
    pin_mean: bool,
) -> SorResult {
    let nx = coeffs.nx;
    let ny = coeffs.ny;
    assert_eq!(
        (solution.width, solution.height),
        (nx, ny),
        "solution field must match the coefficient grid"
    );
    let omega = omega.clamp(0.5, 1.99);

    let mut iterations = 0;
    let mut residual = f64::INFINITY;
    let mut converged = false;

    for sweep in 0..max_iter.max(1) {
        iterations = sweep + 1;
        // One Gauss-Seidel / SOR sweep over every cell.
        for j in 0..ny {
            for i in 0..nx {
                let ap = coeffs.ap.at(i, j);
                if ap.abs() < 1e-30 {
                    // A degenerate cell (no stencil) — leave it alone.
                    continue;
                }
                // Gather the four neighbour contributions using the
                // freshest values (cells to the west / south were
                // already updated this sweep — that is the
                // Gauss-Seidel acceleration).
                let mut sum = coeffs.b.at(i, j);
                if i + 1 < nx {
                    sum += coeffs.ae.at(i, j) * solution.at(i + 1, j);
                }
                if i > 0 {
                    sum += coeffs.aw.at(i, j) * solution.at(i - 1, j);
                }
                if j + 1 < ny {
                    sum += coeffs.an.at(i, j) * solution.at(i, j + 1);
                }
                if j > 0 {
                    sum += coeffs.as_.at(i, j) * solution.at(i, j - 1);
                }
                let gs = sum / ap;
                // SOR: over-relax from the old value toward the
                // Gauss-Seidel target.
                let old = solution.at(i, j);
                solution.set(i, j, old + omega * (gs - old));
            }
        }

        // Pin the gauge of a fully-singular (all-Neumann) system by
        // removing the mean. Skipped when a Dirichlet anchor already
        // fixes the gauge — subtracting the mean would fight it.
        if pin_mean {
            let mean: f64 = solution.data.iter().sum::<f64>() / solution.data.len() as f64;
            for v in solution.data.iter_mut() {
                *v -= mean;
            }
        }

        // Residual: ‖A·p' − b‖₂ over the interior stencil.
        residual = poisson_residual(coeffs, solution);
        if residual <= tol {
            converged = true;
            break;
        }
    }

    SorResult {
        iterations,
        residual,
        converged,
    }
}

/// The residual L2 norm of a candidate solution to the Poisson system —
/// `‖aP·p' − (aE·p'_E + … + b)‖₂`, normalised by the cell count.
pub fn poisson_residual(coeffs: &PoissonCoeffs, solution: &Field) -> f64 {
    let nx = coeffs.nx;
    let ny = coeffs.ny;
    let mut sum_sq = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let ap = coeffs.ap.at(i, j);
            if ap.abs() < 1e-30 {
                continue;
            }
            let mut rhs = coeffs.b.at(i, j);
            if i + 1 < nx {
                rhs += coeffs.ae.at(i, j) * solution.at(i + 1, j);
            }
            if i > 0 {
                rhs += coeffs.aw.at(i, j) * solution.at(i - 1, j);
            }
            if j + 1 < ny {
                rhs += coeffs.an.at(i, j) * solution.at(i, j + 1);
            }
            if j > 0 {
                rhs += coeffs.as_.at(i, j) * solution.at(i, j - 1);
            }
            let r = ap * solution.at(i, j) - rhs;
            sum_sq += r * r;
        }
    }
    let n = (nx * ny) as f64;
    if n > 0.0 {
        (sum_sq / n).sqrt()
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build the coefficients of a standard 5-point Laplacian on an
    /// `n × n` grid with cell size `h` and a homogeneous-Neumann
    /// boundary on `p'` — every interior cell gets `1/h²` on each
    /// in-domain neighbour, `aP` the sum.
    fn laplacian(n: usize, h: f64) -> PoissonCoeffs {
        let mut c = PoissonCoeffs::zeros(n, n);
        let w = 1.0 / (h * h);
        for j in 0..n {
            for i in 0..n {
                let mut ap = 0.0;
                if i + 1 < n {
                    c.ae.set(i, j, w);
                    ap += w;
                }
                if i > 0 {
                    c.aw.set(i, j, w);
                    ap += w;
                }
                if j + 1 < n {
                    c.an.set(i, j, w);
                    ap += w;
                }
                if j > 0 {
                    c.as_.set(i, j, w);
                    ap += w;
                }
                c.ap.set(i, j, ap);
            }
        }
        c
    }

    #[test]
    fn sor_solves_a_known_poisson_problem() {
        // Build a Laplacian, pick a smooth zero-mean target field φ*,
        // compute the source b = A·φ*, then confirm SOR recovers φ*.
        let n = 12;
        let h = 1.0 / n as f64;
        let mut c = laplacian(n, h);
        // A zero-mean target: φ*(i,j) = cos(πx)·cos(πy)-ish, sampled.
        let mut target = Field::zeros(n, n);
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * h;
                let y = (j as f64 + 0.5) * h;
                target.set(
                    i,
                    j,
                    (std::f64::consts::PI * x).cos() * (std::f64::consts::PI * y).cos(),
                );
            }
        }
        // Make it exactly zero-mean (the singular system's gauge).
        let mean: f64 = target.data.iter().sum::<f64>() / target.data.len() as f64;
        for v in target.data.iter_mut() {
            *v -= mean;
        }
        // Source b such that A·target = b → b(i,j) = aP·t − Σ a·t_nb.
        for j in 0..n {
            for i in 0..n {
                let mut nb = 0.0;
                if i + 1 < n {
                    nb += c.ae.at(i, j) * target.at(i + 1, j);
                }
                if i > 0 {
                    nb += c.aw.at(i, j) * target.at(i - 1, j);
                }
                if j + 1 < n {
                    nb += c.an.at(i, j) * target.at(i, j + 1);
                }
                if j > 0 {
                    nb += c.as_.at(i, j) * target.at(i, j - 1);
                }
                c.b.set(i, j, c.ap.at(i, j) * target.at(i, j) - nb);
            }
        }
        // Solve from a zero start.
        let mut sol = Field::zeros(n, n);
        let res = solve_sor(&c, &mut sol, 1.7, 1e-10, 5000, true);
        assert!(
            res.converged,
            "SOR should converge, residual {}",
            res.residual
        );
        // The recovered field must match the zero-mean target.
        let mut max_err = 0.0f64;
        for k in 0..sol.data.len() {
            max_err = max_err.max((sol.data[k] - target.data[k]).abs());
        }
        assert!(max_err < 1e-4, "SOR solution error {max_err} too large");
    }

    #[test]
    fn sor_reduces_the_residual_monotonically_overall() {
        // A simple source; the residual after the solve must be far
        // below the residual of the zero initial guess.
        let n = 10;
        let h = 0.1;
        let mut c = laplacian(n, h);
        // A localised source bump (zero-sum so the Neumann system is
        // compatible).
        c.b.set(2, 2, 1.0);
        c.b.set(7, 7, -1.0);
        let mut sol = Field::zeros(n, n);
        let initial = poisson_residual(&c, &sol);
        let res = solve_sor(&c, &mut sol, 1.6, 1e-9, 3000, true);
        assert!(
            res.residual < initial,
            "solve must reduce the residual: {initial} → {}",
            res.residual
        );
        assert!(
            res.residual < 1e-6,
            "final residual {} not small",
            res.residual
        );
    }

    #[test]
    fn sor_solution_is_zero_mean() {
        // The singular all-Neumann system is pinned to zero mean.
        let n = 8;
        let mut c = laplacian(n, 0.125);
        c.b.set(1, 1, 0.5);
        c.b.set(6, 6, -0.5);
        let mut sol = Field::zeros(n, n);
        solve_sor(&c, &mut sol, 1.5, 1e-9, 2000, true);
        let mean: f64 = sol.data.iter().sum::<f64>() / sol.data.len() as f64;
        assert!(mean.abs() < 1e-9, "solution mean {mean} should be ~0");
    }

    #[test]
    fn sor_converges_faster_than_plain_gauss_seidel() {
        // Over-relaxation should beat ω = 1 on the same problem — the
        // whole point of SOR. Compare the residual after a fixed,
        // small sweep budget.
        let n = 20;
        let h = 1.0 / n as f64;
        let mut c = laplacian(n, h);
        c.b.set(5, 5, 1.0);
        c.b.set(14, 14, -1.0);

        let mut gs = Field::zeros(n, n);
        let gs_res = solve_sor(&c, &mut gs, 1.0, 0.0, 40, true).residual;
        let mut sor = Field::zeros(n, n);
        let sor_res = solve_sor(&c, &mut sor, 1.8, 0.0, 40, true).residual;
        assert!(
            sor_res < gs_res,
            "SOR residual {sor_res} should beat Gauss-Seidel {gs_res} in a fixed budget"
        );
    }
}
