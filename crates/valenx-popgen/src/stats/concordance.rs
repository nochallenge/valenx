//! Genotype concordance between two matrices.
//!
//! Concordance is the fraction of matching allelic calls across all
//! samples and sites — a standard QC metric for comparing two genotyping
//! experiments, imputed versus observed data, or two sequencing platforms
//! run on the same individuals.

use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;

/// Genotype concordance — the fraction of the `n_samples × n_sites` calls that match between
/// two equally-shaped matrices, in `[0, 1]` (1.0 = perfect agreement). Both matrices must have
/// identical dimensions.
///
/// # Errors
/// Returns [`PopgenError::Dimension`] if the matrices differ in sample count or site count.
pub fn genotype_concordance(a: &GenotypeMatrix, b: &GenotypeMatrix) -> Result<f64> {
    if a.n_samples() != b.n_samples() {
        return Err(PopgenError::dimension(
            a.n_samples(),
            b.n_samples(),
            "sample count",
        ));
    }
    if a.n_sites() != b.n_sites() {
        return Err(PopgenError::dimension(
            a.n_sites(),
            b.n_sites(),
            "site count",
        ));
    }
    let total = a.n_samples() * a.n_sites();
    if total == 0 {
        return Ok(0.0);
    }
    let matches: usize = a
        .rows()
        .iter()
        .zip(b.rows())
        .map(|(r1, r2)| r1.iter().zip(r2).filter(|(x, y)| x == y).count())
        .sum();
    Ok(matches as f64 / total as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_concordance_fraction_of_matching_calls() {
        // Identical 2×3 matrices → perfect concordance 1.0.
        let m =
            GenotypeMatrix::from_rows(vec![vec![0, 1, 0], vec![1, 0, 1]], vec![10.0, 20.0, 30.0])
                .unwrap();
        assert!((genotype_concordance(&m, &m).unwrap() - 1.0).abs() < 1e-12);
        // One mismatch out of 4 calls → 3/4 = 0.75.
        let a = GenotypeMatrix::from_rows(vec![vec![0, 1], vec![1, 0]], vec![10.0, 20.0]).unwrap();
        let b = GenotypeMatrix::from_rows(vec![vec![0, 0], vec![1, 0]], vec![10.0, 20.0]).unwrap();
        assert!((genotype_concordance(&a, &b).unwrap() - 0.75).abs() < 1e-12);
        // Fully discordant → 0.0.
        let z1 = GenotypeMatrix::from_rows(vec![vec![0, 0]], vec![10.0, 20.0]).unwrap();
        let z2 = GenotypeMatrix::from_rows(vec![vec![1, 1]], vec![10.0, 20.0]).unwrap();
        assert!((genotype_concordance(&z1, &z2).unwrap() - 0.0).abs() < 1e-12);
        // Dimension mismatch → Err.
        let c = GenotypeMatrix::from_rows(vec![vec![0, 1, 0]], vec![10.0, 20.0, 30.0]).unwrap();
        assert!(genotype_concordance(&a, &c).is_err());
    }
}
