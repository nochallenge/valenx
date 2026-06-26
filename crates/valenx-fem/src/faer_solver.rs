//! High-performance sparse-direct linear-solve backend for the native
//! FEA path, built on the pure-Rust [`faer`] linear-algebra crate.
//!
//! valenx-fem assembles the global stiffness matrix `K` (sparse,
//! symmetric positive-definite after the penalty boundary conditions)
//! and solves `K·u = f`. The original — and still default for small
//! systems — path factorises `K` with
//! [`nalgebra_sparse::factorization::CscCholesky`]. That solver does **no
//! fill-reducing reordering**, so on a large, badly-ordered mesh the
//! Cholesky factor fills in and the factorisation becomes both slow and
//! memory-hungry.
//!
//! [`faer`]'s sparse Cholesky ([`faer::sparse::linalg::solvers::Llt`])
//! computes a fill-reducing permutation (approximate minimum degree) as
//! part of its symbolic analysis and runs a cache-friendly,
//! optionally-parallel supernodal numeric factorisation — typically much
//! faster on the larger systems valenx solves. It is a pure-Rust crate
//! (no C/Fortran, no LAPACK link) and has no GUI/egui coupling, so it
//! slots into the headless solver cleanly.
//!
//! This module is **non-breaking**: it adds a *selectable* backend
//! ([`SolverBackend`]) and a single SPD-solve entry point
//! ([`solve_spd`]). The legacy `CscCholesky` path is retained as
//! [`SolverBackend::Legacy`] and remains the explicit fallback; the faer
//! path is [`SolverBackend::Faer`]. The two are validated to agree to
//! solver precision on the canonical patch-test / cantilever problems
//! (see the tests in [`crate::native_solver`]).

use nalgebra::DVector;
use nalgebra_sparse::CscMatrix;

/// Selectable linear-solve backend for the SPD stiffness system `K·u =
/// f`.
///
/// Both backends solve the *same* assembled system; they differ only in
/// the factorisation engine. The result agrees to solver precision (the
/// difference is rounding in two different — but both backward-stable —
/// Cholesky implementations), so the choice is purely a
/// performance/robustness trade-off, never a change in the physics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SolverBackend {
    /// The original [`nalgebra_sparse::factorization::CscCholesky`]
    /// sparse Cholesky. No fill-reducing reordering; simplest and the
    /// safe fallback. Best on small systems where reordering overhead
    /// would not pay off.
    Legacy,
    /// The [`faer`] supernodal sparse Cholesky
    /// ([`faer::sparse::linalg::solvers::Llt`]) with an approximate-
    /// minimum-degree fill-reducing permutation. The high-performance
    /// path; the default for large systems.
    #[default]
    Faer,
}

/// Degree-of-freedom count at or above which [`SolverBackend::default`]
/// (used by the convenience entry points) picks the faer supernodal
/// path. Below it the legacy `CscCholesky` is used — the faer symbolic
/// analysis / reordering overhead does not pay off on tiny systems, and
/// keeping the well-worn legacy path for small problems minimises any
/// behaviour change for the existing unit tests.
pub const FAER_DOF_THRESHOLD: usize = 1_000;

/// Pick the backend for a system of `n_dof` unknowns: faer for large
/// systems (`n_dof ≥ `[`FAER_DOF_THRESHOLD`]), the legacy Cholesky for
/// small ones.
#[must_use]
pub fn default_backend_for(n_dof: usize) -> SolverBackend {
    if n_dof >= FAER_DOF_THRESHOLD {
        SolverBackend::Faer
    } else {
        SolverBackend::Legacy
    }
}

/// Failure of the sparse SPD solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaerSolveError {
    /// The matrix was not positive-definite (the Cholesky factorisation
    /// failed) — typically an under-constrained / singular system the
    /// cheap up-front checks did not catch.
    NotPositiveDefinite,
}

/// Solve the symmetric-positive-definite system `a·x = b` and return `x`
/// as an `n`-vector, using the requested [`SolverBackend`].
///
/// `a` is the assembled global stiffness in CSC form (already stiffened
/// by the penalty boundary conditions, hence SPD); `b` is the load
/// vector. The two backends produce the same `x` to solver precision.
///
/// # Errors
///
/// [`FaerSolveError::NotPositiveDefinite`] if the factorisation fails.
pub fn solve_spd(
    a: &CscMatrix<f64>,
    b: &DVector<f64>,
    backend: SolverBackend,
) -> Result<DVector<f64>, FaerSolveError> {
    match backend {
        SolverBackend::Legacy => solve_spd_legacy(a, b),
        SolverBackend::Faer => solve_spd_faer(a, b),
    }
}

/// Legacy path: [`nalgebra_sparse::factorization::CscCholesky`].
fn solve_spd_legacy(a: &CscMatrix<f64>, b: &DVector<f64>) -> Result<DVector<f64>, FaerSolveError> {
    use nalgebra_sparse::factorization::CscCholesky;
    let chol = CscCholesky::factor(a).map_err(|_| FaerSolveError::NotPositiveDefinite)?;
    let x = chol.solve(b);
    Ok(x.column(0).into_owned())
}

/// faer path: convert the CSC matrix to a [`faer::sparse::SparseColMat`]
/// and solve with the supernodal sparse Cholesky
/// ([`faer::sparse::linalg::solvers::Llt`]).
///
/// `K` is symmetric, so only the lower triangle is handed to faer (the
/// upper is reconstructed from it via [`faer::Side::Lower`]); the
/// CscCholesky on the legacy side likewise only reads one triangle. The
/// fill-reducing reordering and the numeric factorisation are internal
/// to faer.
fn solve_spd_faer(a: &CscMatrix<f64>, b: &DVector<f64>) -> Result<DVector<f64>, FaerSolveError> {
    use faer::linalg::solvers::Solve;
    use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
    use faer::sparse::{SparseColMat, Triplet};
    use faer::{Mat, Side};

    let n = a.nrows();

    // Build the lower-triangular triplet list `(row, col, value)` from
    // the CSC entries. Taking only `row >= col` halves the work handed
    // to faer and matches the `Side::Lower` directive below; duplicate
    // coordinates (none here — the CSC has already summed them) would be
    // accumulated by `try_new_from_triplets`.
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::with_capacity(a.nnz());
    for (r, c, &v) in a.triplet_iter() {
        if r >= c {
            triplets.push(Triplet::new(r, c, v));
        }
    }

    let mat: SparseColMat<usize, f64> = SparseColMat::try_new_from_triplets(n, n, &triplets)
        .map_err(|_| FaerSolveError::NotPositiveDefinite)?;

    // Supernodal sparse Cholesky of the lower triangle: symbolic
    // analysis (fill-reducing ordering) first, then the numeric
    // factorisation. The symbolic step can only fail on a structurally
    // invalid pattern; the numeric step fails if `K` is not actually
    // positive-definite — both map to `NotPositiveDefinite`.
    let symbolic = SymbolicLlt::try_new(mat.symbolic(), Side::Lower)
        .map_err(|_| FaerSolveError::NotPositiveDefinite)?;
    let llt = Llt::try_new_with_symbolic(symbolic, mat.as_ref(), Side::Lower)
        .map_err(|_| FaerSolveError::NotPositiveDefinite)?;

    // RHS as a single-column dense matrix; `Llt::solve` returns the
    // owned solution column.
    let rhs = Mat::from_fn(n, 1, |i, _| b[i]);
    let sol = llt.solve(&rhs);

    let mut x = DVector::<f64>::zeros(n);
    for i in 0..n {
        x[i] = sol[(i, 0)];
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra_sparse::CooMatrix;

    /// Build a small SPD CSC matrix and a RHS with a known solution, and
    /// confirm the faer and legacy backends agree with each other and
    /// with the analytic answer.
    #[test]
    fn faer_matches_legacy_on_known_spd_system() {
        // A 3×3 SPD system:
        //   [ 4  1  0 ] [x]   [ 1 ]
        //   [ 1  3  1 ] [y] = [ 2 ]
        //   [ 0  1  2 ] [z]   [ 3 ]
        let mut coo = CooMatrix::<f64>::new(3, 3);
        for &(i, j, v) in &[
            (0, 0, 4.0),
            (0, 1, 1.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (1, 2, 1.0),
            (2, 1, 1.0),
            (2, 2, 2.0),
        ] {
            coo.push(i, j, v);
        }
        let csc = CscMatrix::from(&coo);
        let b = DVector::from_row_slice(&[1.0, 2.0, 3.0]);

        let x_legacy = solve_spd(&csc, &b, SolverBackend::Legacy).unwrap();
        let x_faer = solve_spd(&csc, &b, SolverBackend::Faer).unwrap();

        // Residual ‖A·x − b‖ must be ~0 for the faer solution.
        let ax = &csc * &x_faer;
        let resid = (&ax - &b).norm();
        assert!(resid < 1e-10, "faer residual too large: {resid}");

        // The two backends must agree to solver precision.
        let diff = (&x_faer - &x_legacy).norm();
        assert!(
            diff < 1e-10,
            "faer vs legacy disagree: {diff} (faer {x_faer:?}, legacy {x_legacy:?})"
        );
    }

    #[test]
    fn default_backend_switches_on_size() {
        assert_eq!(default_backend_for(10), SolverBackend::Legacy);
        assert_eq!(default_backend_for(FAER_DOF_THRESHOLD), SolverBackend::Faer);
        assert_eq!(SolverBackend::default(), SolverBackend::Faer);
    }
}
