//! A from-scratch simplex linear-programming solver.
//!
//! Flux-balance analysis is, at its core, a linear program: maximise
//! a linear objective subject to linear equality constraints
//! (`S·v = 0`) and per-variable bounds. Rather than add an LP-solver
//! dependency, this module implements the **two-phase primal simplex
//! method** over a dense tableau — enough to solve the
//! genome-scale-*toy* problems a v1 metabolic model poses.
//!
//! The public entry point [`solve_lp`] takes a problem in the
//! canonical form
//!
//! ```text
//!     maximise   cᵀ x
//!     subject to A x { ≤ | = | ≥ } b
//!                lo ≤ x ≤ hi
//! ```
//!
//! Variable bounds (including negative lower bounds, essential for
//! reversible reactions) are handled by a **variable-substitution
//! shift**: each `x_i` is replaced by `lo_i + x_i'` with `x_i' ≥ 0`,
//! and a finite upper bound becomes an extra `≤` row. Phase 1 drives
//! out artificial variables to find a basic feasible solution; phase 2
//! optimises. Bland's rule is used for pivot selection so the method
//! cannot cycle.
//!
//! ## v1 caveats
//!
//! Dense tableau, `O(rows·cols)` per pivot — fine for the hundreds of
//! reactions a v1 metabolic model has, not for a 10 000-reaction
//! genome-scale reconstruction (that wants a sparse revised simplex or
//! an interior-point method). The solver is exact arithmetic-wise only
//! up to floating-point round-off; a small tolerance governs
//! degeneracy and feasibility tests.

/// Sense of a linear constraint row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintSense {
    /// `A_row · x ≤ b`.
    Le,
    /// `A_row · x = b`.
    Eq,
    /// `A_row · x ≥ b`.
    Ge,
}

/// Outcome of an LP solve.
#[derive(Debug, Clone, PartialEq)]
pub enum LpStatus {
    /// An optimal solution was found.
    Optimal,
    /// The constraints admit no feasible point.
    Infeasible,
    /// The objective is unbounded above on the feasible region.
    Unbounded,
}

/// A solved (or diagnosed) linear program.
#[derive(Debug, Clone, PartialEq)]
pub struct LpSolution {
    /// Solve status.
    pub status: LpStatus,
    /// Optimal variable values (empty unless `status == Optimal`).
    pub x: Vec<f64>,
    /// Optimal objective value (`0.0` unless `status == Optimal`).
    pub objective: f64,
}

/// A linear program in canonical maximise form.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearProgram {
    /// Objective coefficients `c`, length = number of variables.
    pub c: Vec<f64>,
    /// Constraint matrix `A`, row-major (`rows × n_vars`).
    pub a: Vec<Vec<f64>>,
    /// Right-hand sides `b`, one per row.
    pub b: Vec<f64>,
    /// Constraint senses, one per row.
    pub sense: Vec<ConstraintSense>,
    /// Per-variable lower bounds.
    pub lo: Vec<f64>,
    /// Per-variable upper bounds (use `f64::INFINITY` for unbounded).
    pub hi: Vec<f64>,
}

const EPS: f64 = 1e-9;

impl LinearProgram {
    /// Solve the program by the two-phase primal simplex method.
    // Index loops are intentional: the tableau construction below
    // walks several parallel `Vec`s by a shared row / column index.
    #[allow(clippy::needless_range_loop)]
    pub fn solve(&self) -> LpSolution {
        let n = self.c.len();
        if self.a.iter().any(|r| r.len() != n)
            || self.b.len() != self.a.len()
            || self.sense.len() != self.a.len()
            || self.lo.len() != n
            || self.hi.len() != n
        {
            return LpSolution {
                status: LpStatus::Infeasible,
                x: Vec::new(),
                objective: 0.0,
            };
        }

        // --- Shift variables so every lower bound is zero. ----------
        // x_i = lo_i + x'_i,  x'_i >= 0.  An infinite lower bound is
        // split into a free variable as x'_i - x''_i.
        let mut col_of: Vec<(usize, bool)> = Vec::new(); // (orig var, is_neg_part)
        let mut shift: Vec<f64> = Vec::new(); // additive shift per orig var
        for i in 0..n {
            shift.push(if self.lo[i].is_finite() {
                self.lo[i]
            } else {
                0.0
            });
            col_of.push((i, false));
            if !self.lo[i].is_finite() {
                // free variable: add a negative part
                col_of.push((i, true));
            }
        }
        let m_cols = col_of.len();

        // Build rows: original constraints (shifted) + finite upper
        // bounds as <= rows.
        let mut rows: Vec<Vec<f64>> = Vec::new();
        let mut rhs: Vec<f64> = Vec::new();
        let mut senses: Vec<ConstraintSense> = Vec::new();

        let expand = |orig_row: &[f64]| -> Vec<f64> {
            col_of
                .iter()
                .map(|&(v, neg)| if neg { -orig_row[v] } else { orig_row[v] })
                .collect()
        };

        for (ri, row) in self.a.iter().enumerate() {
            let exp = expand(row);
            // shift adjusts the rhs: A(x'+shift) {<=,=,>=} b
            let shifted_b: f64 =
                self.b[ri] - row.iter().zip(&shift).map(|(a, s)| a * s).sum::<f64>();
            rows.push(exp);
            rhs.push(shifted_b);
            senses.push(self.sense[ri]);
        }
        // Finite upper bounds.
        for i in 0..n {
            if self.hi[i].is_finite() {
                let mut row = vec![0.0; m_cols];
                // x_i = shift_i + (pos - neg). Upper bound on x_i.
                for (c, &(v, neg)) in col_of.iter().enumerate() {
                    if v == i {
                        row[c] = if neg { -1.0 } else { 1.0 };
                    }
                }
                rows.push(row);
                rhs.push(self.hi[i] - shift[i]);
                senses.push(ConstraintSense::Le);
            }
        }

        // Objective in the shifted space (constant term dropped, added
        // back at the end).
        let obj_const: f64 = self.c.iter().zip(&shift).map(|(c, s)| c * s).sum();
        let obj: Vec<f64> = col_of
            .iter()
            .map(|&(v, neg)| if neg { -self.c[v] } else { self.c[v] })
            .collect();

        match simplex_two_phase(&rows, &rhs, &senses, &obj) {
            Some((xp, val)) => {
                // Map back to original variables.
                let mut x = shift.clone();
                for (c, &(v, neg)) in col_of.iter().enumerate() {
                    if neg {
                        x[v] -= xp[c];
                    } else {
                        x[v] += xp[c];
                    }
                }
                LpSolution {
                    status: LpStatus::Optimal,
                    x,
                    objective: val + obj_const,
                }
            }
            None => LpSolution {
                status: LpStatus::Infeasible,
                x: Vec::new(),
                objective: 0.0,
            },
        }
    }
}

/// Standalone convenience wrapper around [`LinearProgram::solve`].
pub fn solve_lp(lp: &LinearProgram) -> LpSolution {
    lp.solve()
}

/// Two-phase primal simplex over a dense tableau. All variables are
/// `>= 0`. Returns `Some((x, objective))` on optimality, `None` if
/// infeasible. Unboundedness is reported by an objective of
/// `+infinity`.
// Index loops are intentional: the tableau row operations index
// several parallel `Vec`s (`tab`, `obj`, `basis`) by a shared index.
#[allow(clippy::needless_range_loop)]
fn simplex_two_phase(
    a: &[Vec<f64>],
    b: &[f64],
    sense: &[ConstraintSense],
    c: &[f64],
) -> Option<(Vec<f64>, f64)> {
    let m = a.len();
    let n = c.len();
    if m == 0 {
        // No constraints: optimum is 0 at the origin unless some c>0
        // (unbounded). All x>=0.
        if c.iter().any(|&ci| ci > EPS) {
            return Some((vec![0.0; n], f64::INFINITY));
        }
        return Some((vec![0.0; n], 0.0));
    }

    // Normalise so every b >= 0 (flip the row & its sense if not).
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    let mut sense = sense.to_vec();
    for i in 0..m {
        if b[i] < 0.0 {
            for v in a[i].iter_mut() {
                *v = -*v;
            }
            b[i] = -b[i];
            sense[i] = match sense[i] {
                ConstraintSense::Le => ConstraintSense::Ge,
                ConstraintSense::Ge => ConstraintSense::Le,
                ConstraintSense::Eq => ConstraintSense::Eq,
            };
        }
    }

    // Tableau columns: n structural, then one slack/surplus per row,
    // then artificials where needed.
    let mut slack_count = 0;
    let mut needs_artificial = vec![false; m];
    for i in 0..m {
        match sense[i] {
            ConstraintSense::Le => slack_count += 1,
            ConstraintSense::Ge => {
                slack_count += 1;
                needs_artificial[i] = true;
            }
            ConstraintSense::Eq => needs_artificial[i] = true,
        }
    }
    let n_art = needs_artificial.iter().filter(|&&x| x).count();
    let total = n + slack_count + n_art;

    // Build the tableau: m rows x (total + 1) including rhs column.
    let mut tab = vec![vec![0.0; total + 1]; m];
    let mut basis = vec![0usize; m];
    let mut slack_idx = n;
    let mut art_idx = n + slack_count;
    for i in 0..m {
        for j in 0..n {
            tab[i][j] = a[i][j];
        }
        match sense[i] {
            ConstraintSense::Le => {
                tab[i][slack_idx] = 1.0;
                basis[i] = slack_idx;
                slack_idx += 1;
            }
            ConstraintSense::Ge => {
                tab[i][slack_idx] = -1.0; // surplus
                slack_idx += 1;
                tab[i][art_idx] = 1.0;
                basis[i] = art_idx;
                art_idx += 1;
            }
            ConstraintSense::Eq => {
                tab[i][art_idx] = 1.0;
                basis[i] = art_idx;
                art_idx += 1;
            }
        }
        tab[i][total] = b[i];
    }

    // --- Phase 1: minimise the sum of artificial variables. --------
    if n_art > 0 {
        // Phase-1 objective row (we minimise, so store as cost coeffs).
        let mut obj = vec![0.0; total + 1];
        for j in (n + slack_count)..total {
            obj[j] = 1.0;
        }
        // Price out the artificials currently in the basis.
        for (i, &bv) in basis.iter().enumerate() {
            if bv >= n + slack_count {
                for j in 0..=total {
                    obj[j] -= tab[i][j];
                }
            }
        }
        // Phase 1: every column (artificials included) may enter.
        run_phase(&mut tab, &mut basis, &mut obj, total, total);
        // Infeasible if residual artificial cost remains.
        if -obj[total] > 1e-7 {
            return None;
        }
        // Drive any artificial still basic (at zero) out of the basis.
        for i in 0..m {
            if basis[i] >= n + slack_count {
                if let Some(pcol) = (0..n + slack_count).find(|&j| tab[i][j].abs() > EPS) {
                    pivot(&mut tab, &mut basis, i, pcol, total);
                }
            }
        }
    }

    // --- Phase 2: maximise c. Store as a minimise-of-(-c) row. -----
    let mut obj = vec![0.0; total + 1];
    for (j, &cj) in c.iter().enumerate() {
        obj[j] = -cj;
    }
    // Price out basic structural variables.
    for (i, &bv) in basis.iter().enumerate() {
        if bv < total {
            let coeff = obj[bv];
            if coeff != 0.0 {
                for j in 0..=total {
                    obj[j] -= coeff * tab[i][j];
                }
            }
        }
    }
    // Phase 2: artificial columns (indices >= n + slack_count) are
    // *not* eligible to enter — they only seeded phase 1.
    let unbounded = run_phase(&mut tab, &mut basis, &mut obj, total, n + slack_count);
    if unbounded {
        return Some((vec![0.0; n], f64::INFINITY));
    }

    // Extract the structural solution.
    let mut x = vec![0.0; n];
    for (i, &bv) in basis.iter().enumerate() {
        if bv < n {
            x[bv] = tab[i][total];
        }
    }
    let objective = obj[total]; // = c·x at optimum
    Some((x, objective))
}

/// Run primal-simplex pivots until optimal (or unbounded). Returns
/// `true` if the problem is unbounded.
///
/// `enter_limit` is the number of columns eligible to *enter* the
/// basis: in phase 2 it excludes the artificial columns so an
/// artificial — having served only to seed phase 1 — can never
/// re-enter and silently re-introduce constraint infeasibility. Row
/// operations still span the full tableau width (`total`).
// Index loops index the objective row and tableau rows in lockstep.
#[allow(clippy::needless_range_loop)]
fn run_phase(
    tab: &mut [Vec<f64>],
    basis: &mut [usize],
    obj: &mut [f64],
    total: usize,
    enter_limit: usize,
) -> bool {
    let m = tab.len();
    for _ in 0..10_000 {
        // Bland's rule: smallest *eligible* index with a negative
        // reduced cost.
        let entering = (0..enter_limit).find(|&j| obj[j] < -EPS);
        let Some(col) = entering else {
            return false; // optimal
        };
        // Ratio test.
        let mut leave: Option<usize> = None;
        let mut best_ratio = f64::INFINITY;
        for i in 0..m {
            let aij = tab[i][col];
            if aij > EPS {
                let ratio = tab[i][total] / aij;
                if ratio < best_ratio - EPS
                    || (ratio < best_ratio + EPS
                        && leave.map(|l| basis[i] < basis[l]).unwrap_or(true))
                {
                    best_ratio = ratio;
                    leave = Some(i);
                }
            }
        }
        let Some(row) = leave else {
            return true; // unbounded
        };
        pivot(tab, basis, row, col, total);
        // Eliminate the entering column from the objective row.
        let coeff = obj[col];
        if coeff != 0.0 {
            for j in 0..=total {
                obj[j] -= coeff * tab[row][j];
            }
        }
    }
    false
}

/// Gauss-Jordan pivot on `tab[row][col]`, updating the basis.
fn pivot(tab: &mut [Vec<f64>], basis: &mut [usize], row: usize, col: usize, total: usize) {
    let piv = tab[row][col];
    for j in 0..=total {
        tab[row][j] /= piv;
    }
    let m = tab.len();
    for i in 0..m {
        if i != row {
            let factor = tab[i][col];
            if factor != 0.0 {
                for j in 0..=total {
                    tab[i][j] -= factor * tab[row][j];
                }
            }
        }
    }
    basis[row] = col;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maximise_simple_le_program() {
        // max x + y  s.t. x + 2y <= 14, 3x - y >= 0, x - y <= 2,
        // x,y >= 0. Classic textbook LP; optimum at (6, 4) -> 10.
        let lp = LinearProgram {
            c: vec![1.0, 1.0],
            a: vec![vec![1.0, 2.0], vec![3.0, -1.0], vec![1.0, -1.0]],
            b: vec![14.0, 0.0, 2.0],
            sense: vec![
                ConstraintSense::Le,
                ConstraintSense::Ge,
                ConstraintSense::Le,
            ],
            lo: vec![0.0, 0.0],
            hi: vec![f64::INFINITY, f64::INFINITY],
        };
        let sol = lp.solve();
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 10.0).abs() < 1e-6, "obj {}", sol.objective);
        assert!((sol.x[0] - 6.0).abs() < 1e-6);
        assert!((sol.x[1] - 4.0).abs() < 1e-6);
    }

    #[test]
    fn equality_constraint_is_respected() {
        // max x  s.t. x + y = 5, x,y in [0,10]. Optimum x=5.
        let lp = LinearProgram {
            c: vec![1.0, 0.0],
            a: vec![vec![1.0, 1.0]],
            b: vec![5.0],
            sense: vec![ConstraintSense::Eq],
            lo: vec![0.0, 0.0],
            hi: vec![10.0, 10.0],
        };
        let sol = lp.solve();
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.x[0] - 5.0).abs() < 1e-6);
        assert!((sol.x[0] + sol.x[1] - 5.0).abs() < 1e-6);
    }

    #[test]
    fn negative_lower_bound_is_handled() {
        // max x  s.t. x <= 3, x >= -5  (bound only). Optimum x = 3.
        let lp = LinearProgram {
            c: vec![1.0],
            a: vec![],
            b: vec![],
            sense: vec![],
            lo: vec![-5.0],
            hi: vec![3.0],
        };
        let sol = lp.solve();
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.x[0] - 3.0).abs() < 1e-6);
        // Minimising direction should hit the negative bound.
        let lp_min = LinearProgram {
            c: vec![-1.0],
            ..lp
        };
        let sol_min = lp_min.solve();
        assert!((sol_min.x[0] - (-5.0)).abs() < 1e-6);
    }

    #[test]
    fn detects_unbounded_objective() {
        // max x  with no upper bound and no constraint.
        let lp = LinearProgram {
            c: vec![1.0],
            a: vec![],
            b: vec![],
            sense: vec![],
            lo: vec![0.0],
            hi: vec![f64::INFINITY],
        };
        let sol = lp.solve();
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(sol.objective.is_infinite());
    }

    #[test]
    fn detects_infeasible_program() {
        // x >= 10 and x <= 1 simultaneously.
        let lp = LinearProgram {
            c: vec![1.0],
            a: vec![vec![1.0], vec![1.0]],
            b: vec![10.0, 1.0],
            sense: vec![ConstraintSense::Ge, ConstraintSense::Le],
            lo: vec![0.0],
            hi: vec![f64::INFINITY],
        };
        let sol = lp.solve();
        assert_eq!(sol.status, LpStatus::Infeasible);
    }

    #[test]
    fn equality_system_with_zero_rhs_is_feasible() {
        // A chain of zero-RHS equality constraints (the shape an FBA
        // steady-state S·v = 0 produces): v1 = v2 = v3, v1 <= 10,
        // maximise v3. The optimum pins every flux to 10. This guards
        // the regression where artificial variables could re-enter the
        // basis in phase 2 and return an infeasible "optimum".
        let lp = LinearProgram {
            c: vec![0.0, 0.0, 1.0],
            a: vec![vec![1.0, -1.0, 0.0], vec![0.0, 1.0, -1.0]],
            b: vec![0.0, 0.0],
            sense: vec![ConstraintSense::Eq, ConstraintSense::Eq],
            lo: vec![0.0, 0.0, 0.0],
            hi: vec![10.0, 1000.0, 1000.0],
        };
        let sol = lp.solve();
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 10.0).abs() < 1e-6, "obj {}", sol.objective);
        for (i, &xi) in sol.x.iter().enumerate() {
            assert!((xi - 10.0).abs() < 1e-6, "x[{i}] = {xi}");
        }
    }

    #[test]
    fn dimension_mismatch_reports_infeasible() {
        let lp = LinearProgram {
            c: vec![1.0, 1.0],
            a: vec![vec![1.0]], // wrong width
            b: vec![1.0],
            sense: vec![ConstraintSense::Le],
            lo: vec![0.0, 0.0],
            hi: vec![1.0, 1.0],
        };
        assert_eq!(lp.solve().status, LpStatus::Infeasible);
    }
}
