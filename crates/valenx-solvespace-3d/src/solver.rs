//! Newton-Raphson solver with Levenberg-Marquardt damping — direct port
//! of `valenx_sketch::solver` to the 3D entity / constraint set.
//!
//! The math is identical: assemble `J` and `r`, solve
//! `(JᵀJ + λI)·Δx = -Jᵀr`, accept the step if it reduced the residual
//! norm, otherwise bump λ and try again.

use serde::Serialize;

use crate::error::Solve3DError;
use crate::sketch::Sketch3D;

/// User-facing solver outcome.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize)]
pub enum SolverStatus {
    /// All constraints satisfied within tolerance.
    Converged,
    /// Reached max iterations without converging.
    MaxIterations,
}

/// Tunable knobs.
#[derive(Copy, Clone, Debug)]
pub struct SolverConfig {
    /// Stop when residual L2 norm is below this.
    pub tol: f64,
    /// Initial Levenberg-Marquardt damping.
    pub lambda_init: f64,
    /// Maximum Newton iterations.
    pub max_iter: usize,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            tol: 1e-9,
            lambda_init: 1e-4,
            max_iter: 80,
        }
    }
}

/// Diagnostic summary.
#[derive(Copy, Clone, Debug, Serialize)]
pub struct SolverDiagnostics {
    /// Total scalar residuals.
    pub n_residuals: usize,
    /// Number of free variables.
    pub n_variables: usize,
    /// `n_residuals - n_variables`. Negative = under-constrained.
    pub dof_balance: i64,
}

impl SolverDiagnostics {
    fn from(s: &Sketch3D) -> Self {
        Self {
            n_residuals: s.total_residuals(),
            n_variables: s.vars.len(),
            dof_balance: s.total_residuals() as i64 - s.vars.len() as i64,
        }
    }
}

/// Final solver report.
#[derive(Clone, Debug, Serialize)]
pub struct SolverReport {
    /// Final status.
    pub status: SolverStatus,
    /// Iterations performed.
    pub iterations: usize,
    /// Final residual L2 norm.
    pub residual_norm: f64,
    /// Diagnostics snapshot.
    pub diagnostics: SolverDiagnostics,
}

/// Assemble the residual vector.
pub fn assemble_residuals(s: &Sketch3D) -> nalgebra::DVector<f64> {
    let n = s.total_residuals();
    let mut out = nalgebra::DVector::<f64>::zeros(n);
    let mut row_offset = 0;
    for c in &s.constraints {
        let k = c.n_residuals();
        c.residuals(s, &mut out.as_mut_slice()[row_offset..row_offset + k]);
        row_offset += k;
    }
    out
}

/// Assemble the dense Jacobian — small enough for sketches we expect
/// (Phase 53 v1 targets ≤ a few dozen entities; cubic dense solve is
/// fine).
pub fn assemble_jacobian(s: &Sketch3D) -> nalgebra::DMatrix<f64> {
    let n_rows = s.total_residuals();
    let n_cols = s.vars.len();
    let mut j = nalgebra::DMatrix::<f64>::zeros(n_rows, n_cols);
    let mut row_offset = 0;
    for c in &s.constraints {
        let mut tri = Vec::new();
        c.jacobian_triplets(s, &mut tri);
        for (row, col, val) in tri {
            j[(row_offset + row, col)] += val;
        }
        row_offset += c.n_residuals();
    }
    j
}

/// One Newton-LM step. Returns `None` if the linear solve is singular.
pub fn newton_step(s: &Sketch3D, lambda: f64) -> Option<Vec<f64>> {
    let n_vars = s.vars.len();
    if n_vars == 0 || s.constraints.is_empty() {
        return Some(vec![0.0; n_vars]);
    }
    let j = assemble_jacobian(s);
    let r = assemble_residuals(s);
    let jt = j.transpose();
    let mut jtj = &jt * &j;
    for i in 0..n_vars {
        jtj[(i, i)] += lambda;
    }
    let neg_jtr = -&jt * r;
    let chol = nalgebra::linalg::Cholesky::new(jtj)?;
    let step = chol.solve(&neg_jtr);
    Some(step.iter().copied().collect())
}

/// Drive Newton-LM iterations until converged or budget exhausted.
pub fn solve(s: &mut Sketch3D, cfg: SolverConfig) -> Result<SolverReport, Solve3DError> {
    let diagnostics = SolverDiagnostics::from(s);
    let mut lambda = cfg.lambda_init;
    let mut last_norm = assemble_residuals(s).norm();
    for iter in 1..=cfg.max_iter {
        if last_norm < cfg.tol {
            return Ok(SolverReport {
                status: SolverStatus::Converged,
                iterations: iter - 1,
                residual_norm: last_norm,
                diagnostics,
            });
        }
        let Some(step) = newton_step(s, lambda) else {
            lambda *= 10.0;
            if lambda > 1e12 {
                return Err(Solve3DError::Singular(iter));
            }
            continue;
        };
        let saved = s.vars.clone();
        for (i, dx) in step.iter().enumerate() {
            s.vars[i] += dx;
        }
        let new_norm = assemble_residuals(s).norm();
        if new_norm < last_norm {
            last_norm = new_norm;
            lambda = (lambda * 0.5).max(1e-12);
        } else {
            s.vars = saved;
            lambda *= 10.0;
            if lambda > 1e12 {
                return Ok(SolverReport {
                    status: SolverStatus::MaxIterations,
                    iterations: iter,
                    residual_norm: last_norm,
                    diagnostics,
                });
            }
        }
    }
    Ok(SolverReport {
        status: if last_norm < cfg.tol {
            SolverStatus::Converged
        } else {
            SolverStatus::MaxIterations
        },
        iterations: cfg.max_iter,
        residual_norm: last_norm,
        diagnostics,
    })
}
