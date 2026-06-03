//! Newton-Raphson constraint solver with Levenberg-Marquardt damping.

use nalgebra_sparse::{CooMatrix, CsrMatrix};
use serde::Serialize;

use crate::error::SketchError;
use crate::sketch::Sketch;

/// User-facing solver outcome.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize)]
pub enum SolverStatus {
    /// All constraints satisfied within tolerance.
    Converged,
    /// Reached max iterations without converging (sketch may be
    /// inconsistent or simply hard).
    MaxIterations,
}

/// Tunable solver knobs.
#[derive(Copy, Clone, Debug)]
pub struct SolverConfig {
    /// Stop when `||residual||₂ < tol`.
    pub tol: f64,
    /// Initial Levenberg-Marquardt damping.
    pub lambda_init: f64,
    /// Max Newton iterations.
    pub max_iter: usize,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            tol: 1e-9,
            lambda_init: 1e-4,
            max_iter: 50,
        }
    }
}

/// Solver outcome, including diagnostic snapshots.
#[derive(Clone, Debug, Serialize)]
pub struct SolverReport {
    /// Final status.
    pub status: SolverStatus,
    /// Iterations performed.
    pub iterations: usize,
    /// Final residual L2 norm.
    pub residual_norm: f64,
    /// Diagnostics (constraint vs DOF counts).
    pub diagnostics: SolverDiagnostics,
}

/// Constraint vs DOF accounting.
#[derive(Copy, Clone, Debug, Serialize)]
pub struct SolverDiagnostics {
    /// Total scalar residual equations.
    pub n_residuals: usize,
    /// Total variables (DOF).
    pub n_variables: usize,
    /// `n_residuals - n_variables`. Negative => under-constrained;
    /// zero => exactly constrained; positive => over-constrained
    /// (system may have no solution or be redundant).
    pub dof_balance: i64,
}

impl SolverDiagnostics {
    fn from(sketch: &Sketch) -> Self {
        let n_residuals = sketch.total_residuals();
        let n_variables = sketch.vars.len();
        Self {
            n_residuals,
            n_variables,
            dof_balance: n_residuals as i64 - n_variables as i64,
        }
    }
}

/// Compute one Newton-Raphson step with Levenberg-Marquardt damping.
///
/// Solves `(JᵀJ + λI)·Δx = -Jᵀr` for Δx where λ is the LM damping
/// factor. Returns the step vector (same length as `sketch.vars`).
///
/// Phase 12D: variables registered in [`Sketch::frozen_vars`] (typically
/// external-geometry references) are excluded — their corresponding
/// Jacobian columns are zeroed before the normal-equations solve so the
/// LM update leaves them at their initial value.
pub fn newton_step(sketch: &Sketch, lambda: f64) -> Option<Vec<f64>> {
    let j = assemble_jacobian(sketch);
    let r = assemble_residuals(sketch);
    let n_vars = sketch.vars.len();
    if n_vars == 0 || sketch.constraints.is_empty() {
        return Some(vec![0.0; n_vars]);
    }
    let mut j_dense = nalgebra::DMatrix::<f64>::zeros(j.nrows(), j.ncols());
    for (row, col, val) in j.triplet_iter() {
        // Phase 12D: drop columns of frozen variables.
        if sketch.is_var_frozen(col) {
            continue;
        }
        j_dense[(row, col)] += *val;
    }
    let jt = j_dense.transpose();
    let mut jtj = &jt * &j_dense;
    for i in 0..n_vars {
        if sketch.is_var_frozen(i) {
            // Pin the diagonal so the Cholesky stays PD and the row/col
            // map to a zero step.
            jtj[(i, i)] = 1.0;
        } else {
            jtj[(i, i)] += lambda;
        }
    }
    let neg_jtr = -&jt * r;
    let chol = nalgebra::linalg::Cholesky::new(jtj)?;
    let mut step = chol.solve(&neg_jtr);
    // Force frozen-var entries to exact zero (the all-zero row above
    // already drives this, but being explicit avoids any rounding).
    for i in 0..n_vars {
        if sketch.is_var_frozen(i) {
            step[i] = 0.0;
        }
    }
    Some(step.iter().copied().collect())
}

/// Iterate Newton steps with LM damping until converged or max_iter.
pub fn solve(sketch: &mut Sketch, cfg: SolverConfig) -> Result<SolverReport, SketchError> {
    let diagnostics = SolverDiagnostics::from(sketch);
    let mut lambda = cfg.lambda_init;
    let mut last_norm = assemble_residuals(sketch).norm();
    for iter in 1..=cfg.max_iter {
        if last_norm < cfg.tol {
            return Ok(SolverReport {
                status: SolverStatus::Converged,
                iterations: iter - 1,
                residual_norm: last_norm,
                diagnostics,
            });
        }
        let Some(step) = newton_step(sketch, lambda) else {
            // JᵀJ + λI was not PD; bump damping and retry.
            lambda *= 10.0;
            continue;
        };
        // Apply step and accept if residual decreased; else reject and
        // bump damping (classical LM strategy).
        let saved = sketch.vars.clone();
        for (i, dx) in step.iter().enumerate() {
            sketch.vars[i] += dx;
        }
        let new_norm = assemble_residuals(sketch).norm();
        if new_norm < last_norm {
            last_norm = new_norm;
            lambda = (lambda * 0.5).max(1e-12);
        } else {
            sketch.vars = saved;
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

/// Assemble the global Jacobian matrix in row-major sparse form.
/// Rows are constraint residual equations; columns are sketch
/// variables (indices into [`Sketch::vars`]).
pub fn assemble_jacobian(sketch: &Sketch) -> CsrMatrix<f64> {
    let n_rows = sketch.total_residuals();
    let n_cols = sketch.vars.len();
    let mut coo = CooMatrix::<f64>::new(n_rows, n_cols);
    let mut row_offset = 0;
    for c in &sketch.constraints {
        let mut triplets = Vec::new();
        c.jacobian_triplets(sketch, &mut triplets);
        for (row, col, val) in triplets {
            coo.push(row_offset + row, col, val);
        }
        row_offset += c.n_residuals();
    }
    CsrMatrix::from(&coo)
}

/// Assemble the residual vector.
pub fn assemble_residuals(sketch: &Sketch) -> nalgebra::DVector<f64> {
    let n = sketch.total_residuals();
    let mut out = nalgebra::DVector::<f64>::zeros(n);
    let mut row_offset = 0;
    for c in &sketch.constraints {
        let k = c.n_residuals();
        let slice = out.as_mut_slice();
        c.residuals(sketch, &mut slice[row_offset..row_offset + k]);
        row_offset += k;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::Constraint;

    #[test]
    fn coincident_jacobian_has_4_nonzero_entries() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 1.0);
        s.add_constraint(Constraint::Coincident { a, b });
        let j = assemble_jacobian(&s);
        assert_eq!(j.nrows(), 2);
        assert_eq!(j.ncols(), 4);
        assert_eq!(j.nnz(), 4);
    }

    #[test]
    fn residual_vector_has_expected_size_and_values() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        s.add_constraint(Constraint::Coincident { a, b });
        let r = assemble_residuals(&s);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], 3.0);
        assert_eq!(r[1], 4.0);
    }

    #[test]
    fn newton_step_reduces_residual_for_simple_coincident() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        s.add_constraint(Constraint::Coincident { a, b });
        let initial_norm = assemble_residuals(&s).norm();
        let step = newton_step(&s, 1e-6).unwrap();
        // Apply step
        for (i, dx) in step.iter().enumerate() {
            s.vars[i] += dx;
        }
        let final_norm = assemble_residuals(&s).norm();
        assert!(final_norm < initial_norm, "step failed to reduce residual");
    }

    #[test]
    fn solver_converges_for_coincident_pair() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        s.add_constraint(Constraint::Coincident { a, b });
        let report = solve(&mut s, Default::default()).unwrap();
        assert!(matches!(report.status, SolverStatus::Converged));
        // After solving, b should now coincide with a (at (0, 0)).
        let pa = s.point_at(a).unwrap();
        let pb = s.point_at(b).unwrap();
        let (ax, ay) = pa.read(&s.vars);
        let (bx, by) = pb.read(&s.vars);
        assert!((ax - bx).abs() < 1e-9);
        assert!((ay - by).abs() < 1e-9);
    }

    #[test]
    fn solver_handles_horizontal_plus_distance() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(2.0, 1.5); // not horizontal yet
        let line = s.add_line(a, b).unwrap();
        s.add_constraint(Constraint::Horizontal(line));
        s.add_constraint(Constraint::Distance { a, b, target: 5.0 });
        let report = solve(&mut s, Default::default()).unwrap();
        assert!(
            matches!(report.status, SolverStatus::Converged),
            "{:?}",
            report.status
        );
        // Line should be horizontal and length 5.
        let pa = s.point_at(a).unwrap();
        let pb = s.point_at(b).unwrap();
        let (ax, ay) = pa.read(&s.vars);
        let (bx, by) = pb.read(&s.vars);
        assert!((ay - by).abs() < 1e-6, "not horizontal: ay={ay}, by={by}");
        let dist = ((bx - ax).powi(2) + (by - ay).powi(2)).sqrt();
        assert!((dist - 5.0).abs() < 1e-6, "length wrong: {dist}");
    }

    #[test]
    fn underconstrained_sketch_reports_negative_dof_balance() {
        let mut s = Sketch::new();
        let _ = s.add_point(0.0, 0.0);
        let _ = s.add_point(5.0, 0.0);
        // No constraints — 4 vars, 0 residuals → DOF balance = -4.
        let report = solve(&mut s, Default::default()).unwrap();
        assert_eq!(report.diagnostics.dof_balance, -4);
        assert_eq!(report.status, SolverStatus::Converged); // 0 residuals trivially converge
    }

    #[test]
    fn inconsistent_constraints_dont_converge() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        s.add_constraint(Constraint::Coincident { a, b });
        s.add_constraint(Constraint::Distance { a, b, target: 5.0 });
        let report = solve(&mut s, Default::default()).unwrap();
        // Can't satisfy both — solver hits MaxIterations or stays
        // with positive residual.
        assert!(report.residual_norm > 1e-3 || report.status == SolverStatus::MaxIterations);
    }

    /// Phase 12D Task 37: variables in `frozen_vars` don't move during
    /// solve. Pin a point at (3, 4) using `frozen_vars`, give it a
    /// coincident constraint to another (movable) point — the
    /// movable point should travel to the frozen one, NOT the other
    /// way around.
    #[test]
    fn frozen_vars_do_not_move_during_solve() {
        let mut s = Sketch::new();
        let frozen = s.add_point(3.0, 4.0);
        let movable = s.add_point(0.0, 0.0);
        let pf = s.point_at(frozen).unwrap();
        s.frozen_vars.insert(pf.x_var);
        s.frozen_vars.insert(pf.y_var);
        s.add_constraint(Constraint::Coincident {
            a: frozen,
            b: movable,
        });
        let report = solve(&mut s, Default::default()).unwrap();
        assert!(
            matches!(report.status, SolverStatus::Converged),
            "{report:?}"
        );
        // Frozen point still at (3, 4).
        let pf = s.point_at(frozen).unwrap();
        let (fx, fy) = pf.read(&s.vars);
        assert!(
            (fx - 3.0).abs() < 1e-9 && (fy - 4.0).abs() < 1e-9,
            "frozen point moved: ({fx}, {fy})"
        );
        // Movable point landed near (3, 4).
        let pm = s.point_at(movable).unwrap();
        let (mx, my) = pm.read(&s.vars);
        assert!(
            (mx - 3.0).abs() < 1e-6 && (my - 4.0).abs() < 1e-6,
            "movable point didn't reach frozen: ({mx}, {my})"
        );
    }

    #[test]
    fn solver_handles_right_triangle() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(0.5, 0.5); // perturbed; should land at (0, 1)-ish after solve
        let ab = s.add_line(a, b).unwrap();
        let ac = s.add_line(a, c).unwrap();
        // Pin a at origin
        s.add_constraint(Constraint::Coincident { a, b: a }); // no-op but legal
        s.add_constraint(Constraint::Perpendicular { a: ab, b: ac });
        s.add_constraint(Constraint::EqualLength { a: ab, b: ac });
        s.add_constraint(Constraint::Distance { a, b, target: 1.0 });
        let report = solve(&mut s, Default::default()).unwrap();
        assert!(
            matches!(report.status, SolverStatus::Converged),
            "{report:?}"
        );
        // After solving: ab and ac perpendicular, equal length, ab length = 1.
        let line_ab = s.line_at(ab).unwrap();
        let line_ac = s.line_at(ac).unwrap();
        let (dx1, dy1) = line_ab.direction(&s.vars);
        let (dx2, dy2) = line_ac.direction(&s.vars);
        let dot = dx1 * dx2 + dy1 * dy2;
        assert!(dot.abs() < 1e-6, "not perpendicular: dot = {dot}");
        let len1 = line_ab.length(&s.vars);
        let len2 = line_ac.length(&s.vars);
        assert!((len1 - len2).abs() < 1e-6);
        assert!((len1 - 1.0).abs() < 1e-6);
    }
}
