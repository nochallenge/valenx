//! Linear-algebra helpers for the SCF — orthogonalisation and the
//! generalized symmetric eigenproblem.
//!
//! The Roothaan-Hall equations `FC = SCε` are a *generalized* symmetric
//! eigenproblem. The standard route is to transform them into an
//! ordinary one with an orthogonalising matrix `X` satisfying
//! `XᵀSX = 1`:
//!
//! ```text
//! F' = Xᵀ F X      F' C' = C' ε      C = X C'
//! ```
//!
//! Two choices of `X`:
//!
//! - **Symmetric (Löwdin)** orthogonalisation `X = S^{-1/2}` — keeps the
//!   transformed orbitals as close as possible to the originals.
//! - **Canonical** orthogonalisation `X = U s^{-1/2}` — drops the
//!   eigenvectors of `S` whose eigenvalue is below a threshold, which
//!   is the robust fallback when the basis is near-linearly-dependent
//!   and `S^{-1/2}` would blow up.
//!
//! Both are built from `nalgebra`'s symmetric eigendecomposition.

use crate::error::{QchemError, Result};
use nalgebra::{DMatrix, DVector, SymmetricEigen};

/// Eigenvalues of the overlap matrix below this threshold are treated as
/// linear dependencies by [`canonical_orthogonalizer`].
pub const LINEAR_DEPENDENCE_THRESHOLD: f64 = 1.0e-6;

/// The result of orthogonalising the overlap matrix.
#[derive(Clone, Debug)]
pub struct Orthogonalizer {
    /// The orthogonalising matrix `X` (`n × m`, `m ≤ n`).
    pub x: DMatrix<f64>,
    /// `true` when canonical orthogonalisation dropped one or more
    /// near-singular directions.
    pub linear_dependence_dropped: bool,
}

impl Orthogonalizer {
    /// Number of retained orthogonal functions `m`.
    #[inline]
    pub fn n_retained(&self) -> usize {
        self.x.ncols()
    }
}

/// Symmetric eigendecomposition of a real symmetric matrix, returned as
/// `(eigenvalues ascending, eigenvectors as columns)`.
pub fn symmetric_eigh(a: &DMatrix<f64>) -> (DVector<f64>, DMatrix<f64>) {
    let se = SymmetricEigen::new(a.clone());
    // nalgebra does not guarantee an order — sort ascending.
    let n = se.eigenvalues.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| se.eigenvalues[i].partial_cmp(&se.eigenvalues[j]).unwrap());
    let mut vals = DVector::zeros(n);
    let mut vecs = DMatrix::zeros(n, n);
    for (new, &old) in idx.iter().enumerate() {
        vals[new] = se.eigenvalues[old];
        vecs.set_column(new, &se.eigenvectors.column(old));
    }
    (vals, vecs)
}

/// Build the symmetric (Löwdin) orthogonaliser `X = S^{-1/2}`.
///
/// When the smallest overlap eigenvalue is below
/// [`LINEAR_DEPENDENCE_THRESHOLD`] the basis is near-linearly-dependent
/// and `S^{-1/2}` is numerically unsound; this routine then falls back
/// to [`canonical_orthogonalizer`] automatically.
///
/// # Errors
///
/// Returns [`QchemError::InvalidInput`] when `S` has a non-positive
/// eigenvalue (not a valid overlap matrix).
pub fn symmetric_orthogonalizer(s: &DMatrix<f64>) -> Result<Orthogonalizer> {
    let (vals, vecs) = symmetric_eigh(s);
    let min = vals[0];
    if min <= 0.0 {
        return Err(QchemError::invalid(format!(
            "overlap matrix has a non-positive eigenvalue ({min:.3e}) — \
             not a valid basis"
        )));
    }
    if min < LINEAR_DEPENDENCE_THRESHOLD {
        return canonical_orthogonalizer(s);
    }
    // X = U s^{-1/2} Uᵀ.
    let n = s.nrows();
    let mut inv_sqrt = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        inv_sqrt[(i, i)] = 1.0 / vals[i].sqrt();
    }
    let x = &vecs * inv_sqrt * vecs.transpose();
    Ok(Orthogonalizer {
        x,
        linear_dependence_dropped: false,
    })
}

/// Build the canonical orthogonaliser `X = U s^{-1/2}`, discarding the
/// eigenvectors of `S` whose eigenvalue is below
/// [`LINEAR_DEPENDENCE_THRESHOLD`].
///
/// The result has `m ≤ n` columns — the working dimension shrinks when
/// linear dependencies are removed.
///
/// # Errors
///
/// Returns [`QchemError::InvalidInput`] when every eigenvalue is below
/// threshold (a degenerate basis).
pub fn canonical_orthogonalizer(s: &DMatrix<f64>) -> Result<Orthogonalizer> {
    let (vals, vecs) = symmetric_eigh(s);
    let n = s.nrows();
    let keep: Vec<usize> = (0..n)
        .filter(|&i| vals[i] >= LINEAR_DEPENDENCE_THRESHOLD)
        .collect();
    if keep.is_empty() {
        return Err(QchemError::invalid(
            "overlap matrix is fully linearly dependent",
        ));
    }
    let m = keep.len();
    let mut x = DMatrix::<f64>::zeros(n, m);
    for (col, &i) in keep.iter().enumerate() {
        let scale = 1.0 / vals[i].sqrt();
        x.set_column(col, &(vecs.column(i) * scale));
    }
    Ok(Orthogonalizer {
        x,
        linear_dependence_dropped: m < n,
    })
}

/// Solve the generalized symmetric eigenproblem `F C = S C ε` given a
/// pre-built orthogonaliser `X`.
///
/// Returns `(ε ascending, C)` where `C` has the *full* `n` rows (one
/// per basis function) and `m = X.ncols()` columns (one per retained
/// orthogonal orbital).
pub fn solve_roothaan(
    f: &DMatrix<f64>,
    ortho: &Orthogonalizer,
) -> (DVector<f64>, DMatrix<f64>) {
    let x = &ortho.x;
    // F' = Xᵀ F X.
    let f_prime = x.transpose() * f * x;
    let (eps, c_prime) = symmetric_eigh(&f_prime);
    // C = X C'.
    let c = x * c_prime;
    (eps, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigh_sorted_ascending() {
        let a = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 5.0]);
        let (vals, _) = symmetric_eigh(&a);
        assert!(vals[0] <= vals[1]);
        assert!((vals[0] - 2.0).abs() < 1.0e-12);
        assert!((vals[1] - 5.0).abs() < 1.0e-12);
    }

    #[test]
    fn symmetric_orthogonalizer_satisfies_xt_s_x_identity() {
        // A well-conditioned 2x2 overlap matrix.
        let s = DMatrix::from_row_slice(2, 2, &[1.0, 0.3, 0.3, 1.0]);
        let o = symmetric_orthogonalizer(&s).unwrap();
        let id = o.x.transpose() * &s * &o.x;
        assert!((id[(0, 0)] - 1.0).abs() < 1.0e-10);
        assert!((id[(1, 1)] - 1.0).abs() < 1.0e-10);
        assert!(id[(0, 1)].abs() < 1.0e-10);
        assert!(!o.linear_dependence_dropped);
    }

    #[test]
    fn near_singular_overlap_drops_a_direction() {
        // S with one tiny eigenvalue → canonical fallback.
        let s = DMatrix::from_row_slice(2, 2, &[1.0, 0.999_999_9, 0.999_999_9, 1.0]);
        let o = symmetric_orthogonalizer(&s).unwrap();
        assert!(o.linear_dependence_dropped);
        assert_eq!(o.n_retained(), 1);
    }

    #[test]
    fn roothaan_recovers_eigenpairs() {
        // With S = identity, FC = SCε reduces to ordinary eigenproblem.
        let s = DMatrix::<f64>::identity(2, 2);
        let f = DMatrix::from_row_slice(2, 2, &[1.0, 0.5, 0.5, 2.0]);
        let o = symmetric_orthogonalizer(&s).unwrap();
        let (eps, c) = solve_roothaan(&f, &o);
        // Reconstruct F from C ε Cᵀ.
        let mut diag = DMatrix::<f64>::zeros(2, 2);
        diag[(0, 0)] = eps[0];
        diag[(1, 1)] = eps[1];
        let recon = &c * diag * c.transpose();
        for i in 0..2 {
            for j in 0..2 {
                assert!((recon[(i, j)] - f[(i, j)]).abs() < 1.0e-10);
            }
        }
    }

    #[test]
    fn non_positive_overlap_is_rejected() {
        let s = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 2.0, 1.0]);
        assert!(symmetric_orthogonalizer(&s).is_err());
    }
}
