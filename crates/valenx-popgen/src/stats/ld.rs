//! Linkage disequilibrium (LD).
//!
//! LD measures the non-random association of alleles at two sites. For
//! a pair of biallelic sites `A` and `B` with derived-allele
//! frequencies `pA`, `pB` and observed derived-derived haplotype
//! frequency `pAB`:
//!
//! - **D** ([`ld_d`]) — the raw covariance `pAB - pA * pB`. Zero under
//!   linkage equilibrium; its sign and magnitude depend on allele
//!   frequencies, so it is hard to compare across site pairs.
//! - **D'** ([`ld_d_prime`]) — `D` scaled by its theoretical maximum
//!   given the allele frequencies, so `D'` lies in `[-1, 1]`
//!   (Lewontin 1964).
//! - **r squared** ([`ld_r_squared`]) — `D^2 / (pA(1-pA) pB(1-pB))`,
//!   the squared correlation between the two sites' alleles; this is
//!   the quantity that relates directly to recombination distance and
//!   association-mapping power.
//!
//! [`ld_matrix`] computes pairwise `r^2` for every pair of sites in a
//! [`GenotypeMatrix`], the input to an LD-decay or LD-block analysis.

use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;

/// The four two-locus quantities for a site pair.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LdStats {
    /// Raw disequilibrium coefficient `D`.
    pub d: f64,
    /// Normalised `D'` in `[-1, 1]`.
    pub d_prime: f64,
    /// Squared correlation `r^2` in `[0, 1]`.
    pub r_squared: f64,
}

/// Computes `D`, `D'` and `r^2` for the pair of sites `(a, b)` of a
/// genotype matrix.
///
/// Both sites are treated as biallelic with the `1` allele derived.
///
/// # Errors
/// [`PopgenError::Invalid`] if either index is out of range or the
/// matrix has fewer than two samples.
pub fn ld_pair(matrix: &GenotypeMatrix, a: usize, b: usize) -> Result<LdStats> {
    let n = matrix.n_samples();
    if n < 2 {
        return Err(PopgenError::invalid(
            "matrix",
            "need at least two samples for LD",
        ));
    }
    if a >= matrix.n_sites() || b >= matrix.n_sites() {
        return Err(PopgenError::invalid("site", "index out of range"));
    }
    let nn = n as f64;
    let mut count_a = 0.0;
    let mut count_b = 0.0;
    let mut count_ab = 0.0;
    for row in matrix.rows() {
        let ha = row[a];
        let hb = row[b];
        count_a += ha as f64;
        count_b += hb as f64;
        if ha == 1 && hb == 1 {
            count_ab += 1.0;
        }
    }
    let p_a = count_a / nn;
    let p_b = count_b / nn;
    let p_ab = count_ab / nn;
    let d = p_ab - p_a * p_b;

    // D' normalisation: divide by D_max, which depends on the sign of
    // D and the marginal frequencies (Lewontin 1964).
    let d_max = if d < 0.0 {
        (p_a * p_b).min((1.0 - p_a) * (1.0 - p_b))
    } else {
        (p_a * (1.0 - p_b)).min((1.0 - p_a) * p_b)
    };
    let d_prime = if d_max.abs() < 1e-12 { 0.0 } else { d / d_max };

    // r^2 = D^2 / (pA qA pB qB).
    let denom = p_a * (1.0 - p_a) * p_b * (1.0 - p_b);
    let r_squared = if denom.abs() < 1e-12 {
        0.0
    } else {
        d * d / denom
    };

    Ok(LdStats {
        d,
        d_prime,
        r_squared,
    })
}

/// The raw disequilibrium coefficient `D` for a site pair.
///
/// # Errors
/// See [`ld_pair`].
pub fn ld_d(matrix: &GenotypeMatrix, a: usize, b: usize) -> Result<f64> {
    Ok(ld_pair(matrix, a, b)?.d)
}

/// The normalised `D'` for a site pair.
///
/// # Errors
/// See [`ld_pair`].
pub fn ld_d_prime(matrix: &GenotypeMatrix, a: usize, b: usize) -> Result<f64> {
    Ok(ld_pair(matrix, a, b)?.d_prime)
}

/// The squared correlation `r^2` for a site pair.
///
/// # Errors
/// See [`ld_pair`].
pub fn ld_r_squared(matrix: &GenotypeMatrix, a: usize, b: usize) -> Result<f64> {
    Ok(ld_pair(matrix, a, b)?.r_squared)
}

/// Computes the full pairwise `r^2` matrix for every pair of sites.
///
/// The result is a symmetric `n_sites x n_sites` matrix; the diagonal
/// is `1.0` for segregating sites and `0.0` for monomorphic ones.
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than two samples.
pub fn ld_matrix(matrix: &GenotypeMatrix) -> Result<Vec<Vec<f64>>> {
    if matrix.n_samples() < 2 {
        return Err(PopgenError::invalid(
            "matrix",
            "need at least two samples for LD",
        ));
    }
    let s = matrix.n_sites();
    let mut out = vec![vec![0.0; s]; s];
    for i in 0..s {
        let seg_i = matrix.is_segregating(i)?;
        out[i][i] = if seg_i { 1.0 } else { 0.0 };
        for j in (i + 1)..s {
            let r2 = ld_pair(matrix, i, j)?.r_squared;
            out[i][j] = r2;
            out[j][i] = r2;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: Vec<Vec<u8>>) -> GenotypeMatrix {
        let cols = rows[0].len();
        let pos: Vec<f64> = (0..cols).map(|i| i as f64).collect();
        GenotypeMatrix::from_rows(rows, pos).unwrap()
    }

    #[test]
    fn perfectly_correlated_sites_have_r2_one() {
        // Two sites with identical columns -> r^2 = 1, D' = 1.
        let m = matrix(vec![vec![1, 1], vec![1, 1], vec![0, 0], vec![0, 0]]);
        let ld = ld_pair(&m, 0, 1).unwrap();
        assert!((ld.r_squared - 1.0).abs() < 1e-9, "r2 = {}", ld.r_squared);
        assert!((ld.d_prime - 1.0).abs() < 1e-9, "D' = {}", ld.d_prime);
        assert!(ld.d > 0.0);
    }

    #[test]
    fn perfectly_anticorrelated_sites() {
        // Site B is the complement of site A -> r^2 = 1, D' = -1.
        let m = matrix(vec![vec![1, 0], vec![1, 0], vec![0, 1], vec![0, 1]]);
        let ld = ld_pair(&m, 0, 1).unwrap();
        assert!((ld.r_squared - 1.0).abs() < 1e-9);
        assert!((ld.d_prime + 1.0).abs() < 1e-9, "D' = {}", ld.d_prime);
        assert!(ld.d < 0.0);
    }

    #[test]
    fn independent_sites_have_low_ld() {
        // Two sites whose alleles are uncorrelated.
        let m = matrix(vec![vec![1, 1], vec![1, 0], vec![0, 1], vec![0, 0]]);
        let ld = ld_pair(&m, 0, 1).unwrap();
        assert!(ld.r_squared.abs() < 1e-9, "r2 = {}", ld.r_squared);
        assert!(ld.d.abs() < 1e-9);
    }

    #[test]
    fn ld_matrix_is_symmetric() {
        let m = matrix(vec![
            vec![1, 1, 0],
            vec![1, 0, 1],
            vec![0, 1, 0],
            vec![0, 0, 1],
        ]);
        let mat = ld_matrix(&m).unwrap();
        assert_eq!(mat.len(), 3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((mat[i][j] - mat[j][i]).abs() < 1e-12);
            }
        }
        // Diagonal of a segregating site is 1.
        assert!((mat[0][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn convenience_wrappers_agree_with_ld_pair() {
        let m = matrix(vec![vec![1, 1], vec![1, 1], vec![0, 0], vec![0, 0]]);
        let full = ld_pair(&m, 0, 1).unwrap();
        assert!((ld_d(&m, 0, 1).unwrap() - full.d).abs() < 1e-12);
        assert!((ld_d_prime(&m, 0, 1).unwrap() - full.d_prime).abs() < 1e-12);
        assert!((ld_r_squared(&m, 0, 1).unwrap() - full.r_squared).abs() < 1e-12);
    }

    #[test]
    fn rejects_out_of_range_sites() {
        let m = matrix(vec![vec![1], vec![0]]);
        assert!(ld_pair(&m, 0, 5).is_err());
    }
}
