//! Discrete differential operators on triangle meshes — the
//! cotangent Laplacian and the mass matrix.
//!
//! These are the building blocks the libigl parameterisation /
//! deformation / geodesics algorithms all rest on. libigl ships them
//! as `igl::cotmatrix` and `igl::massmatrix`; this module is the
//! pure-Rust port.
//!
//! ## The cotangent Laplacian
//!
//! For a triangle mesh the standard discrete Laplace-Beltrami
//! operator `L` has, for an edge `(i, j)` shared by triangles with
//! opposite angles `α` and `β`:
//!
//! ```text
//! L[i][j] = (cot α + cot β) / 2          (off-diagonal)
//! L[i][i] = -Σ_{j ≠ i} L[i][j]            (diagonal — rows sum to 0)
//! ```
//!
//! A boundary edge contributes only the one cotangent it has. This
//! is the FEM stiffness matrix of the piecewise-linear hat basis.
//!
//! ## Solver
//!
//! The systems that arise (`L φ = b` for Poisson, `(M − tL) u = b`
//! for heat flow) are symmetric. For the SPD cases this module uses a
//! dense `nalgebra` Cholesky factorisation; for the indefinite /
//! pinned cases it uses a dense LU. Dense linear algebra is exact and
//! correct for the workshop-scale meshes libigl targets — a sparse
//! factorisation is a speed optimisation that yields identical
//! numbers, tracked as a follow-up for very large meshes.

use nalgebra::{DMatrix, DVector, Vector3};

use crate::triangle::TriMesh;

/// Cotangent of the angle at vertex `apex` in the triangle
/// `apex, a, b` — `cot θ = (e1 · e2) / |e1 × e2|`.
fn cotangent(apex: Vector3<f64>, a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    let e1 = a - apex;
    let e2 = b - apex;
    let cross = e1.cross(&e2).norm();
    if cross < 1e-12 {
        0.0
    } else {
        e1.dot(&e2) / cross
    }
}

/// Build the dense cotangent Laplacian `L` (n × n) for `mesh`.
///
/// Convention: `L` is positive-semidefinite with **non-negative**
/// off-diagonals' negation — i.e. `L[i][i] > 0`, `L[i][j] ≤ 0`, rows
/// sum to zero. This is the sign convention where `L` acts as a
/// discrete `−Δ` (so `L` is PSD), matching libigl's `cotmatrix`
/// negated — the form most solvers want.
pub fn cotangent_laplacian(mesh: &TriMesh) -> DMatrix<f64> {
    let n = mesh.vertices.len();
    let mut l = DMatrix::<f64>::zeros(n, n);
    for tri in &mesh.triangles {
        let (i, j, k) = (tri[0], tri[1], tri[2]);
        let (vi, vj, vk) = (mesh.vertices[i], mesh.vertices[j], mesh.vertices[k]);
        // Angle at each vertex is opposite the corresponding edge.
        // Edge (j, k) is opposite vertex i, etc.
        let cot_i = cotangent(vi, vj, vk); // angle at i, edge (j,k)
        let cot_j = cotangent(vj, vk, vi); // angle at j, edge (k,i)
        let cot_k = cotangent(vk, vi, vj); // angle at k, edge (i,j)
        // Each cotangent weights the *opposite* edge.
        let mut add = |a: usize, b: usize, w: f64| {
            let half = w * 0.5;
            l[(a, b)] -= half;
            l[(b, a)] -= half;
            l[(a, a)] += half;
            l[(b, b)] += half;
        };
        add(j, k, cot_i);
        add(k, i, cot_j);
        add(i, j, cot_k);
    }
    l
}

/// Build the lumped (diagonal) mass matrix as a vector of per-vertex
/// areas — one third of the area of each incident triangle (the
/// "barycentric" lumped mass, libigl's `MASSMATRIX_TYPE_BARYCENTRIC`).
pub fn lumped_mass(mesh: &TriMesh) -> DVector<f64> {
    let n = mesh.vertices.len();
    let mut m = DVector::<f64>::zeros(n);
    for tri in &mesh.triangles {
        let (vi, vj, vk) = (
            mesh.vertices[tri[0]],
            mesh.vertices[tri[1]],
            mesh.vertices[tri[2]],
        );
        let area = (vj - vi).cross(&(vk - vi)).norm() * 0.5;
        let third = area / 3.0;
        for &idx in tri {
            m[idx] += third;
        }
    }
    m
}

/// Mean edge length over every triangle edge — a natural time-step
/// scale for the heat method.
pub fn mean_edge_length(mesh: &TriMesh) -> f64 {
    let mut total = 0.0;
    let mut n = 0usize;
    for tri in &mesh.triangles {
        for k in 0..3 {
            let a = mesh.vertices[tri[k]];
            let b = mesh.vertices[tri[(k + 1) % 3]];
            total += (b - a).norm();
            n += 1;
        }
    }
    if n == 0 {
        1.0
    } else {
        total / n as f64
    }
}

/// Solve the symmetric system `A x = b` for a single right-hand side.
///
/// Tries a Cholesky factorisation first (valid + fastest when `A` is
/// SPD); falls back to LU for indefinite or singular-ish `A`. Returns
/// `None` only if both factorisations fail (a genuinely singular
/// system).
pub fn solve_symmetric(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    if let Some(chol) = a.clone().cholesky() {
        return Some(chol.solve(b));
    }
    let lu = a.clone().lu();
    lu.solve(b)
}

/// Solve `A X = B` for many right-hand sides (the columns of `B`).
pub fn solve_symmetric_multi(a: &DMatrix<f64>, b: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    if let Some(chol) = a.clone().cholesky() {
        return Some(chol.solve(b));
    }
    let lu = a.clone().lu();
    lu.solve(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_quad() -> TriMesh {
        TriMesh {
            vertices: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3]],
        }
    }

    #[test]
    fn laplacian_rows_sum_to_zero() {
        let l = cotangent_laplacian(&unit_quad());
        for i in 0..l.nrows() {
            let row_sum: f64 = (0..l.ncols()).map(|j| l[(i, j)]).sum();
            assert!(row_sum.abs() < 1e-9, "row {i} sum = {row_sum}");
        }
    }

    #[test]
    fn laplacian_is_symmetric() {
        let l = cotangent_laplacian(&unit_quad());
        for i in 0..l.nrows() {
            for j in 0..l.ncols() {
                assert!((l[(i, j)] - l[(j, i)]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn laplacian_diagonal_is_positive() {
        // With the −Δ convention every diagonal entry is > 0 for an
        // interior-or-boundary vertex with at least one incident
        // triangle.
        let l = cotangent_laplacian(&unit_quad());
        for i in 0..l.nrows() {
            assert!(l[(i, i)] > 0.0, "L[{i}][{i}] = {}", l[(i, i)]);
        }
    }

    #[test]
    fn lumped_mass_sums_to_total_area() {
        // The unit square has area 1; the lumped mass entries sum to it.
        let m = lumped_mass(&unit_quad());
        let total: f64 = m.iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "total mass = {total}");
    }

    #[test]
    fn solve_symmetric_recovers_known_solution() {
        // A simple SPD system: A = diag(2, 3), b = (4, 9) → x = (2, 3).
        let a = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 3.0]);
        let b = DVector::from_vec(vec![4.0, 9.0]);
        let x = solve_symmetric(&a, &b).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-9);
        assert!((x[1] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn solve_symmetric_handles_indefinite_via_lu() {
        // An indefinite symmetric matrix — Cholesky fails, LU succeeds.
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 1.0, 0.0]);
        let b = DVector::from_vec(vec![3.0, 5.0]);
        let x = solve_symmetric(&a, &b).unwrap();
        // A x = b ⇒ x = (5, 3).
        assert!((x[0] - 5.0).abs() < 1e-9);
        assert!((x[1] - 3.0).abs() < 1e-9);
    }
}
