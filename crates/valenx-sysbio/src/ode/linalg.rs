//! Small dense linear-algebra helpers shared by the numerical layer.
//!
//! The implicit BDF integrator, the damped-Newton steady-state solver
//! and the bifurcation continuation all need to solve a small dense
//! linear system `A x = b`. Rather than convert to and from
//! `nalgebra` matrices on every Newton iteration, this module offers a
//! direct Gaussian-elimination solver over `Vec<Vec<f64>>` — the same
//! row-major representation the [`OdeSystem`](crate::ode::OdeSystem)
//! Jacobian produces.
//!
//! [`null_space`] (used by the conservation analysis) *does* go
//! through `nalgebra`'s SVD, because a robust rank-revealing
//! factorisation is exactly what `nalgebra` is in the dependency tree
//! for.

use nalgebra::{DMatrix, DVector};

/// Solve the dense linear system `a · x = b` by Gaussian elimination
/// with partial pivoting.
///
/// `a` is an `n × n` row-major matrix, `b` has length `n`. Returns
/// `None` if `a` is singular (a zero pivot survives pivoting) or if
/// the dimensions disagree.
pub fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    if a.len() != n || a.iter().any(|r| r.len() != n) {
        return None;
    }
    // Working copy of the augmented system.
    let mut m: Vec<Vec<f64>> = a
        .iter()
        .zip(b)
        .map(|(row, &bi)| {
            let mut r = row.clone();
            r.push(bi);
            r
        })
        .collect();

    for col in 0..n {
        // Partial pivot: largest magnitude in this column at or below
        // the diagonal.
        let mut pivot = col;
        let mut best = m[col][col].abs();
        for (r, mr) in m.iter().enumerate().skip(col + 1) {
            let v = mr[col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-300 {
            return None; // singular
        }
        m.swap(col, pivot);
        // Eliminate below.
        let pivot_val = m[col][col];
        for r in (col + 1)..n {
            let factor = m[r][col] / pivot_val;
            if factor != 0.0 {
                for c in col..=n {
                    m[r][c] -= factor * m[col][c];
                }
            }
        }
    }
    // Back substitution.
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let row = &m[i];
        let s: f64 = row[n]
            - x.iter()
                .enumerate()
                .skip(i + 1)
                .map(|(j, xj)| row[j] * xj)
                .sum::<f64>();
        x[i] = s / row[i];
    }
    if x.iter().any(|v| !v.is_finite()) {
        return None;
    }
    Some(x)
}

/// An orthonormal basis for the (right) null space of a `rows × cols`
/// matrix `a`.
///
/// `a` is given row-major. Directions whose singular value is below
/// `tol` (relative to the largest) span the null space; they are
/// returned each as a `Vec<f64>` of length `cols`. An empty result
/// means full column rank.
///
/// `nalgebra`'s `svd` is only the *thin* factorisation — its `Vᵀ` has
/// just `min(rows, cols)` rows, so for an under-determined system
/// (`rows < cols`, the usual stoichiometry / conservation case) it
/// cannot represent the null-space directions beyond that rank. To get
/// the **complete** `cols`-dimensional right-singular basis the matrix
/// is first padded with zero rows up to `cols × cols`: zero rows leave
/// `AᵀA` — hence every singular value and right singular vector —
/// unchanged, while making the factorisation square so the thin `Vᵀ`
/// becomes the full one. Because this routes through the SVD directly
/// (not `AᵀA`) the singular values are not squared, so the relative
/// `tol` keeps its usual singular-value meaning.
pub fn null_space(a: &[Vec<f64>], tol: f64) -> Vec<Vec<f64>> {
    let rows = a.len();
    if rows == 0 {
        return Vec::new();
    }
    let cols = a[0].len();
    if cols == 0 {
        return Vec::new();
    }
    // Pad with zero rows to `padded_rows × cols` so the SVD is square
    // (or tall) and its thin `Vᵀ` spans all `cols` right directions.
    let padded_rows = rows.max(cols);
    let mut flat: Vec<f64> = Vec::with_capacity(padded_rows * cols);
    for r in a {
        flat.extend(r.iter().copied());
    }
    flat.resize(padded_rows * cols, 0.0);
    let m = DMatrix::from_row_slice(padded_rows, cols, &flat);
    let svd = m.svd(false, true);
    let v_t = match &svd.v_t {
        Some(v) => v,
        None => return Vec::new(),
    };
    let sv = &svd.singular_values;
    let smax = sv.iter().cloned().fold(0.0_f64, f64::max).max(1e-300);
    let mut basis = Vec::new();
    // Right singular vectors are the rows of Vᵀ; with the square
    // padding there are exactly `cols` of them.
    for j in 0..cols {
        let small = if j < sv.len() {
            sv[j] / smax <= tol
        } else {
            true
        };
        if small {
            let row = v_t.row(j);
            basis.push(row.iter().copied().collect());
        }
    }
    basis
}

/// Solve a least-squares / possibly-singular system `a x = b` via the
/// `nalgebra` SVD pseudo-inverse. Used as a fallback by the
/// steady-state solver when the Newton Jacobian is rank-deficient.
pub fn solve_least_squares(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let rows = a.len();
    if rows == 0 || rows != b.len() {
        return None;
    }
    let cols = a[0].len();
    let flat: Vec<f64> = a.iter().flat_map(|r| r.iter().copied()).collect();
    let m = DMatrix::from_row_slice(rows, cols, &flat);
    let rhs = DVector::from_row_slice(b);
    let svd = m.svd(true, true);
    svd.solve(&rhs, 1e-12)
        .ok()
        .map(|x| x.iter().copied().collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_2x2_system() {
        // [2 1; 1 3] x = [3; 5]  =>  x = [4/5, 7/5].
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let b = vec![3.0, 5.0];
        let x = solve_linear(&a, &b).unwrap();
        assert!((x[0] - 0.8).abs() < 1e-12);
        assert!((x[1] - 1.4).abs() < 1e-12);
    }

    #[test]
    fn singular_system_returns_none() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(solve_linear(&a, &[1.0, 2.0]).is_none());
    }

    #[test]
    fn requires_a_pivot() {
        // Identity solve is the cheapest sanity check.
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let x = solve_linear(&a, &[7.0, -3.0, 2.0]).unwrap();
        assert_eq!(x, vec![7.0, -3.0, 2.0]);
    }

    #[test]
    fn null_space_of_rank_deficient_matrix() {
        // Rows [1 1 0] and [0 0 1] leave a 1-D null space along
        // (1,-1,0)/sqrt(2).
        let a = vec![vec![1.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]];
        let ns = null_space(&a, 1e-9);
        assert_eq!(ns.len(), 1);
        let v = &ns[0];
        // v . [1,1,0] ~ 0  and  v . [0,0,1] ~ 0.
        assert!((v[0] + v[1]).abs() < 1e-9);
        assert!(v[2].abs() < 1e-9);
    }

    #[test]
    fn full_rank_matrix_has_empty_null_space() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(null_space(&a, 1e-9).is_empty());
    }
}
