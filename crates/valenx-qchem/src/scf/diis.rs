//! Pulay DIIS — direct inversion in the iterative subspace.
//!
//! DIIS (Pulay 1980, 1982) accelerates SCF convergence by replacing the
//! current Fock matrix with a linear combination of the Fock matrices
//! from recent iterations that minimises the norm of the corresponding
//! error vectors.
//!
//! ## The error vector
//!
//! At self-consistency the Fock and density matrices commute under the
//! overlap metric, so the off-commutator
//!
//! ```text
//! e = S D F - F D S   (then transformed by Xᵀ e X)
//! ```
//!
//! is exactly zero — it is the natural SCF error vector. DIIS stores
//! `(F, e)` pairs and, each cycle, solves the small linear system
//!
//! ```text
//! [ B  -1 ] [ c ]   [ 0 ]
//! [ -1ᵀ 0 ] [ λ ] = [ -1 ]      B_ij = ⟨e_i, e_j⟩
//! ```
//!
//! for the coefficients `c` that build the extrapolated Fock matrix
//! `F_diis = Σ_i c_i F_i`.

use crate::scf::linalg::Orthogonalizer;
use nalgebra::{DMatrix, DVector};

/// A rolling DIIS history of `(Fock, error)` pairs.
#[derive(Clone, Debug)]
pub struct Diis {
    /// Maximum number of `(F, e)` pairs retained.
    max_vectors: usize,
    /// Stored Fock matrices.
    fock_history: Vec<DMatrix<f64>>,
    /// Stored error vectors (flattened, transformed).
    error_history: Vec<DVector<f64>>,
}

impl Diis {
    /// Create a DIIS accelerator retaining up to `max_vectors` pairs
    /// (8 is the conventional default).
    pub fn new(max_vectors: usize) -> Self {
        Diis {
            max_vectors: max_vectors.max(2),
            fock_history: Vec::new(),
            error_history: Vec::new(),
        }
    }

    /// Number of `(F, e)` pairs currently stored.
    #[inline]
    pub fn len(&self) -> usize {
        self.fock_history.len()
    }

    /// `true` when no pairs are stored yet.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.fock_history.is_empty()
    }

    /// The SCF error vector `Xᵀ (S D F - F D S) X`, returned flattened.
    ///
    /// Transforming by the orthogonaliser `X` makes the error
    /// orthonormal-basis-consistent and removes the overlap metric.
    pub fn error_vector(
        fock: &DMatrix<f64>,
        density: &DMatrix<f64>,
        overlap: &DMatrix<f64>,
        ortho: &Orthogonalizer,
    ) -> DVector<f64> {
        let commutator = overlap * density * fock - fock * density * overlap;
        let transformed = ortho.x.transpose() * commutator * &ortho.x;
        flatten(&transformed)
    }

    /// The largest-magnitude component of an error vector — the scalar
    /// the SCF loop watches for convergence.
    pub fn error_norm(error: &DVector<f64>) -> f64 {
        error.iter().fold(0.0, |m, &e| m.max(e.abs()))
    }

    /// Record a `(Fock, error)` pair, evicting the oldest when the
    /// history is full.
    pub fn push(&mut self, fock: DMatrix<f64>, error: DVector<f64>) {
        self.fock_history.push(fock);
        self.error_history.push(error);
        while self.fock_history.len() > self.max_vectors {
            self.fock_history.remove(0);
            self.error_history.remove(0);
        }
    }

    /// The DIIS-extrapolated Fock matrix.
    ///
    /// With fewer than two stored pairs there is nothing to extrapolate
    /// and the most recent Fock matrix is returned unchanged. When the
    /// `B` system is singular (collinear error vectors) the routine
    /// likewise falls back to the latest Fock matrix.
    pub fn extrapolate(&self) -> Option<DMatrix<f64>> {
        let m = self.fock_history.len();
        if m == 0 {
            return None;
        }
        if m == 1 {
            return Some(self.fock_history[0].clone());
        }
        // Build the (m+1)x(m+1) DIIS system.
        let dim = m + 1;
        let mut b = DMatrix::<f64>::zeros(dim, dim);
        for i in 0..m {
            for j in 0..m {
                b[(i, j)] = self.error_history[i].dot(&self.error_history[j]);
            }
            b[(i, m)] = -1.0;
            b[(m, i)] = -1.0;
        }
        let mut rhs = DVector::<f64>::zeros(dim);
        rhs[m] = -1.0;

        let coeffs = match b.clone().lu().solve(&rhs) {
            Some(c) => c,
            None => return Some(self.fock_history[m - 1].clone()),
        };

        let mut f = DMatrix::<f64>::zeros(
            self.fock_history[0].nrows(),
            self.fock_history[0].ncols(),
        );
        for i in 0..m {
            f += coeffs[i] * &self.fock_history[i];
        }
        Some(f)
    }
}

/// Flatten a matrix into a column vector in row-major order.
fn flatten(m: &DMatrix<f64>) -> DVector<f64> {
    let mut v = DVector::zeros(m.nrows() * m.ncols());
    let mut k = 0;
    for i in 0..m.nrows() {
        for j in 0..m.ncols() {
            v[k] = m[(i, j)];
            k += 1;
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scf::linalg::symmetric_orthogonalizer;

    #[test]
    fn empty_diis_extrapolates_to_none() {
        let d = Diis::new(8);
        assert!(d.is_empty());
        assert!(d.extrapolate().is_none());
    }

    #[test]
    fn single_pair_returns_that_fock() {
        let mut d = Diis::new(8);
        let f = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]);
        d.push(f.clone(), DVector::from_vec(vec![0.1, 0.0, 0.0, 0.1]));
        let ex = d.extrapolate().unwrap();
        assert!((ex - f).norm() < 1.0e-13);
    }

    #[test]
    fn history_is_capped() {
        let mut d = Diis::new(3);
        for k in 0..6 {
            let f = DMatrix::from_element(2, 2, k as f64);
            d.push(f, DVector::from_element(4, k as f64));
        }
        assert_eq!(d.len(), 3);
    }

    #[test]
    fn extrapolation_coefficients_sum_to_one() {
        // The extrapolated Fock matrix of two identical Fock matrices
        // is that Fock matrix (coefficients sum to 1).
        let mut d = Diis::new(8);
        let f = DMatrix::from_row_slice(2, 2, &[1.0, 0.2, 0.2, 1.5]);
        d.push(f.clone(), DVector::from_vec(vec![0.2, 0.0, 0.0, 0.1]));
        d.push(f.clone(), DVector::from_vec(vec![0.05, 0.0, 0.0, 0.03]));
        let ex = d.extrapolate().unwrap();
        assert!((ex - f).norm() < 1.0e-9, "extrapolated F differs");
    }

    #[test]
    fn zero_error_vector_at_self_consistency() {
        // If [S, D, F] commute the error vector is zero.
        let s = DMatrix::<f64>::identity(2, 2);
        let f = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]);
        // A density that commutes with F (also diagonal).
        let d = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 0.0]);
        let o = symmetric_orthogonalizer(&s).unwrap();
        let e = Diis::error_vector(&f, &d, &s, &o);
        assert!(Diis::error_norm(&e) < 1.0e-12);
    }

    #[test]
    fn nonzero_error_when_not_converged() {
        let s = DMatrix::<f64>::identity(2, 2);
        let f = DMatrix::from_row_slice(2, 2, &[1.0, 0.7, 0.7, 2.0]);
        let d = DMatrix::from_row_slice(2, 2, &[1.0, 0.4, 0.4, 1.0]);
        let o = symmetric_orthogonalizer(&s).unwrap();
        let e = Diis::error_vector(&f, &d, &s, &o);
        assert!(Diis::error_norm(&e) > 1.0e-6);
    }
}
