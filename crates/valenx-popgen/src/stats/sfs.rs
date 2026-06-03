//! The site-frequency spectrum (SFS).
//!
//! The SFS is the workhorse summary of a sample's variation: for a
//! sample of `n` chromosomes, `xi[i]` is the number of segregating
//! sites at which the derived allele is present in exactly `i` copies.
//!
//! Two flavours:
//!
//! - **Unfolded** ([`site_frequency_spectrum`]) — needs the ancestral
//!   allele known (which `valenx-popgen`'s `0`/`1` matrices encode);
//!   `xi` runs over derived counts `1..n`.
//! - **Folded** ([`folded_spectrum`]) — when the ancestral state is
//!   unknown, minor-allele counts are used: `eta[i]` for `i =
//!   1..floor(n/2)`, with `eta[i] = xi[i] + xi[n-i]`.
//!
//! Every other diversity statistic in [`crate::stats`] can be written
//! as a linear functional of the SFS, which is why it gets its own
//! module.

use crate::error::Result;
use crate::infer::GenotypeMatrix;

/// A site-frequency spectrum.
#[derive(Clone, Debug, PartialEq)]
pub struct Sfs {
    /// `counts[i]` = number of sites with frequency class `i`.
    ///
    /// For an unfolded spectrum index `i` (`1..=n-1`) is the derived
    /// allele count; index `0` and `n` (fixed classes) are kept for a
    /// uniform layout and are usually zero. For a folded spectrum index
    /// `i` (`1..=floor(n/2)`) is the minor-allele count.
    counts: Vec<usize>,
    /// `true` if this is a folded (minor-allele) spectrum.
    folded: bool,
    /// Sample size `n` (chromosomes).
    n: usize,
}

impl Sfs {
    /// The raw frequency-class counts.
    pub fn counts(&self) -> &[usize] {
        &self.counts
    }

    /// `true` for a folded spectrum.
    pub fn is_folded(&self) -> bool {
        self.folded
    }

    /// Sample size `n`.
    pub fn sample_size(&self) -> usize {
        self.n
    }

    /// Total number of segregating sites the spectrum summarises.
    pub fn segregating_sites(&self) -> usize {
        // Sum over the polymorphic classes only (skip the two fixed
        // classes 0 and n of an unfolded spectrum).
        if self.folded {
            // counts[0] is the monomorphic (minor-allele-count 0) class
            // of a folded spectrum; skip it so monomorphic sites are not
            // miscounted as segregating.
            self.counts.iter().skip(1).sum()
        } else {
            self.counts
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != 0 && i != self.n)
                .map(|(_, &c)| c)
                .sum()
        }
    }

    /// Count of singletons (frequency-class-1 sites).
    pub fn singletons(&self) -> usize {
        self.counts.get(1).copied().unwrap_or(0)
    }
}

/// Computes the **unfolded** site-frequency spectrum of a genotype
/// matrix.
///
/// The matrix's `1` allele is taken to be derived. The returned
/// [`Sfs`] has length `n + 1`; `counts[i]` is the number of sites with
/// `i` derived copies.
///
/// # Errors
/// Propagates [`crate::error::PopgenError`] from column access — the
/// internal calls cannot actually fail for an in-range column.
pub fn site_frequency_spectrum(matrix: &GenotypeMatrix) -> Result<Sfs> {
    let n = matrix.n_samples();
    let mut counts = vec![0usize; n + 1];
    for col in 0..matrix.n_sites() {
        let d = matrix.derived_count(col)?;
        counts[d] += 1;
    }
    Ok(Sfs {
        counts,
        folded: false,
        n,
    })
}

/// Computes the **folded** (minor-allele) site-frequency spectrum.
///
/// The returned [`Sfs`] has length `floor(n/2) + 1`; `counts[i]` is the
/// number of sites whose *minor* allele is in `i` copies.
///
/// # Errors
/// Propagates [`crate::error::PopgenError`] from column access.
pub fn folded_spectrum(matrix: &GenotypeMatrix) -> Result<Sfs> {
    let n = matrix.n_samples();
    let half = n / 2;
    let mut counts = vec![0usize; half + 1];
    for col in 0..matrix.n_sites() {
        let d = matrix.derived_count(col)?;
        let minor = d.min(n - d);
        if minor <= half {
            counts[minor] += 1;
        }
    }
    Ok(Sfs {
        counts,
        folded: true,
        n,
    })
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
    fn unfolded_spectrum_counts_derived_classes() {
        // 4 samples, 3 sites. Site 0 -> 1 derived, site 1 -> 2 derived,
        // site 2 -> 1 derived.
        let m = matrix(vec![
            vec![1, 1, 1],
            vec![0, 1, 0],
            vec![0, 0, 0],
            vec![0, 0, 0],
        ]);
        let sfs = site_frequency_spectrum(&m).unwrap();
        assert_eq!(sfs.sample_size(), 4);
        // Two singletons (sites 0 and 2), one doubleton (site 1).
        assert_eq!(sfs.counts()[1], 2);
        assert_eq!(sfs.counts()[2], 1);
        assert_eq!(sfs.singletons(), 2);
        assert_eq!(sfs.segregating_sites(), 3);
        assert!(!sfs.is_folded());
    }

    #[test]
    fn folded_spectrum_uses_minor_allele() {
        // 4 samples; site with 3 derived copies folds to minor count 1.
        let m = matrix(vec![
            vec![1, 1],
            vec![1, 0],
            vec![1, 0],
            vec![0, 0],
        ]);
        // Site 0: 3 derived -> minor 1. Site 1: 1 derived -> minor 1.
        let sfs = folded_spectrum(&m).unwrap();
        assert!(sfs.is_folded());
        assert_eq!(sfs.counts()[1], 2);
        assert_eq!(sfs.counts().len(), 4 / 2 + 1);
    }

    #[test]
    fn fixed_sites_land_in_the_edge_classes() {
        // A monomorphic-derived site goes to class n; it is not a
        // segregating site.
        let m = matrix(vec![vec![1], vec![1], vec![1], vec![1]]);
        let sfs = site_frequency_spectrum(&m).unwrap();
        assert_eq!(sfs.counts()[4], 1);
        assert_eq!(sfs.segregating_sites(), 0);
    }

    #[test]
    fn folded_segregating_sites_excludes_monomorphic() {
        // 4 samples, 2 sites: site 0 is monomorphic (all ancestral),
        // site 1 is segregating (one derived copy -> minor 1).
        let m = matrix(vec![
            vec![0, 1],
            vec![0, 0],
            vec![0, 0],
            vec![0, 0],
        ]);
        let sfs = folded_spectrum(&m).unwrap();
        assert!(sfs.is_folded());
        // counts[0] = monomorphic class (site 0); counts[1] = single
        // minor-allele class (site 1).
        assert_eq!(sfs.counts()[0], 1);
        assert_eq!(sfs.counts()[1], 1);
        // Only the polymorphic site is segregating — the monomorphic
        // site in counts[0] must NOT be counted.
        assert_eq!(sfs.segregating_sites(), 1);
    }

    #[test]
    fn empty_matrix_has_an_empty_spectrum() {
        let m = GenotypeMatrix::from_rows(vec![vec![], vec![], vec![]], vec![])
            .unwrap();
        let sfs = site_frequency_spectrum(&m).unwrap();
        assert_eq!(sfs.segregating_sites(), 0);
    }
}
