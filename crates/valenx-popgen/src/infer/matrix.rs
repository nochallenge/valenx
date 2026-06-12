//! The [`GenotypeMatrix`] — the central exchange type for statistics
//! and inference.
//!
//! A genotype matrix is `n_samples` rows by `n_sites` columns of
//! biallelic `0`/`1` calls (ancestral / derived). Every summary
//! statistic in [`crate::stats`], the VCF / ms exporters in
//! [`crate::catalog`] and [`crate::infer::abc`] all consume a
//! `GenotypeMatrix`. A [`crate::model::Population`], a coalescent
//! genealogy and a tree sequence can each be projected into one, which
//! is what makes the statistics code model-agnostic.
//!
//! Rows are *haplotypes* (chromosome copies), not individuals — a
//! diploid sample of 10 individuals is 20 rows. Diploid genotype
//! statistics ([`crate::popstats`]) pair consecutive rows.

use crate::error::{PopgenError, Result};
use serde::{Deserialize, Serialize};

/// An `n_samples` x `n_sites` biallelic genotype matrix.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GenotypeMatrix {
    /// One row per haplotype; every row has length `n_sites`.
    rows: Vec<Vec<u8>>,
    /// Genomic position of each site (column). Length `n_sites`.
    positions: Vec<f64>,
}

impl GenotypeMatrix {
    /// Builds a matrix from `0`/`1` rows and a per-column position
    /// vector.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on an empty matrix or a non-`{0,1}`
    /// entry; [`PopgenError::Dimension`] if the rows are ragged or the
    /// position vector length disagrees with the column count.
    pub fn from_rows(rows: Vec<Vec<u8>>, positions: Vec<f64>) -> Result<Self> {
        if rows.is_empty() {
            return Err(PopgenError::invalid("rows", "matrix has no samples"));
        }
        let n_sites = rows[0].len();
        for (i, r) in rows.iter().enumerate() {
            if r.len() != n_sites {
                return Err(PopgenError::dimension(
                    n_sites,
                    r.len(),
                    "genotype matrix row",
                ));
            }
            if r.iter().any(|&v| v > 1) {
                return Err(PopgenError::invalid(
                    "rows",
                    format!("row {i} has a non-biallelic entry"),
                ));
            }
        }
        if positions.len() != n_sites {
            return Err(PopgenError::dimension(
                n_sites,
                positions.len(),
                "site positions",
            ));
        }
        Ok(GenotypeMatrix { rows, positions })
    }

    /// Number of haplotype rows.
    pub fn n_samples(&self) -> usize {
        self.rows.len()
    }

    /// Number of site columns.
    pub fn n_sites(&self) -> usize {
        self.positions.len()
    }

    /// `true` if there are no sites.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// The rows (haplotypes).
    pub fn rows(&self) -> &[Vec<u8>] {
        &self.rows
    }

    /// The per-site genomic positions.
    pub fn positions(&self) -> &[f64] {
        &self.positions
    }

    /// Allelic state of sample `row` at site `col`.
    ///
    /// # Panics
    /// If `row` or `col` is out of range.
    pub fn get(&self, row: usize, col: usize) -> u8 {
        self.rows[row][col]
    }

    /// Derived-allele count at site `col` (column sum).
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `col` is out of range.
    pub fn derived_count(&self, col: usize) -> Result<usize> {
        if col >= self.n_sites() {
            return Err(PopgenError::invalid("col", "site index out of range"));
        }
        Ok(self.rows.iter().map(|r| r[col] as usize).sum())
    }

    /// Derived-allele frequency at site `col` in `[0, 1]`.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `col` is out of range.
    pub fn frequency(&self, col: usize) -> Result<f64> {
        let d = self.derived_count(col)?;
        Ok(d as f64 / self.n_samples().max(1) as f64)
    }

    /// A site is *segregating* if both alleles are present (the derived
    /// count is neither `0` nor `n_samples`).
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `col` is out of range.
    pub fn is_segregating(&self, col: usize) -> Result<bool> {
        let d = self.derived_count(col)?;
        Ok(d != 0 && d != self.n_samples())
    }

    /// Number of segregating sites in the matrix.
    pub fn segregating_sites(&self) -> usize {
        (0..self.n_sites())
            .filter(|&c| self.is_segregating(c).unwrap_or(false))
            .count()
    }

    /// Returns a copy with all monomorphic (non-segregating) columns
    /// removed — the standard pre-processing for `S`-based statistics.
    pub fn drop_monomorphic(&self) -> GenotypeMatrix {
        let keep: Vec<usize> = (0..self.n_sites())
            .filter(|&c| self.is_segregating(c).unwrap_or(false))
            .collect();
        self.select_sites(&keep)
    }

    /// Returns a copy retaining only the listed site columns, in the
    /// given order. Out-of-range indices are skipped.
    pub fn select_sites(&self, cols: &[usize]) -> GenotypeMatrix {
        let cols: Vec<usize> = cols
            .iter()
            .copied()
            .filter(|&c| c < self.n_sites())
            .collect();
        let rows: Vec<Vec<u8>> = self
            .rows
            .iter()
            .map(|r| cols.iter().map(|&c| r[c]).collect())
            .collect();
        let positions: Vec<f64> = cols.iter().map(|&c| self.positions[c]).collect();
        GenotypeMatrix { rows, positions }
    }

    /// Returns a copy retaining only the listed sample rows, in the
    /// given order. Out-of-range indices are skipped — the natural way
    /// to slice out a sub-population for [`crate::stats::fst`].
    pub fn select_samples(&self, samples: &[usize]) -> Result<GenotypeMatrix> {
        let kept: Vec<&Vec<u8>> = samples
            .iter()
            .filter(|&&s| s < self.n_samples())
            .map(|&s| &self.rows[s])
            .collect();
        if kept.is_empty() {
            return Err(PopgenError::invalid(
                "samples",
                "selection retained no rows",
            ));
        }
        Ok(GenotypeMatrix {
            rows: kept.into_iter().cloned().collect(),
            positions: self.positions.clone(),
        })
    }

    /// Transposes to a site-major layout: `n_sites` rows of
    /// `n_samples` calls. Convenient for column-wise LD scans.
    pub fn transpose(&self) -> Vec<Vec<u8>> {
        (0..self.n_sites())
            .map(|c| self.rows.iter().map(|r| r[c]).collect())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> GenotypeMatrix {
        // 4 haplotypes, 3 sites. Site 0 derived count 2, site 1
        // monomorphic (0), site 2 derived count 4 (monomorphic).
        GenotypeMatrix::from_rows(
            vec![vec![1, 0, 1], vec![1, 0, 1], vec![0, 0, 1], vec![0, 0, 1]],
            vec![10.0, 20.0, 30.0],
        )
        .unwrap()
    }

    #[test]
    fn shape_and_counts() {
        let m = sample();
        assert_eq!(m.n_samples(), 4);
        assert_eq!(m.n_sites(), 3);
        assert_eq!(m.derived_count(0).unwrap(), 2);
        assert!((m.frequency(0).unwrap() - 0.5).abs() < 1e-12);
        assert!(m.is_segregating(0).unwrap());
        assert!(!m.is_segregating(1).unwrap());
        assert!(!m.is_segregating(2).unwrap());
        assert_eq!(m.segregating_sites(), 1);
    }

    #[test]
    fn rejects_ragged_and_nonbiallelic() {
        assert!(GenotypeMatrix::from_rows(vec![vec![0, 1], vec![0]], vec![1.0, 2.0]).is_err());
        assert!(GenotypeMatrix::from_rows(vec![vec![0, 2]], vec![1.0, 2.0]).is_err());
        assert!(GenotypeMatrix::from_rows(vec![vec![0, 1]], vec![1.0]).is_err());
    }

    #[test]
    fn drop_monomorphic_keeps_only_segregating() {
        let m = sample().drop_monomorphic();
        assert_eq!(m.n_sites(), 1);
        assert!((m.positions()[0] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn select_samples_and_sites() {
        let m = sample();
        let sub = m.select_samples(&[0, 1]).unwrap();
        assert_eq!(sub.n_samples(), 2);
        assert_eq!(sub.derived_count(0).unwrap(), 2);
        let cols = m.select_sites(&[2, 0]);
        assert_eq!(cols.n_sites(), 2);
        assert!((cols.positions()[0] - 30.0).abs() < 1e-12);
    }

    #[test]
    fn transpose_is_site_major() {
        let t = sample().transpose();
        assert_eq!(t.len(), 3);
        assert_eq!(t[0], vec![1, 1, 0, 0]);
    }
}
