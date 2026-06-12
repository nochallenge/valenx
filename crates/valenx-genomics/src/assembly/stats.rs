//! k-mer-based genome / contig assembly statistics.
//!
//! After an assembly run, the quality of the result is summarised by a
//! handful of length statistics — what QUAST and `assembly-stats`
//! report: N50 / L50 (and the N75 / N90 family), the total assembly
//! size, the contig-count, the longest contig, GC content and a
//! length distribution. This module computes them, plus a k-mer
//! spectrum useful for estimating genome size and coverage.

use crate::error::{GenomicsError, Result};
use std::collections::HashMap;

/// The full assembly statistics summary.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyStats {
    /// Number of contigs.
    pub n_contigs: usize,
    /// Total assembly length (sum of contig lengths).
    pub total_length: u64,
    /// Length of the longest contig.
    pub longest: u64,
    /// Length of the shortest contig.
    pub shortest: u64,
    /// Mean contig length.
    pub mean_length: f64,
    /// **N50** — the length such that contigs of at least this length
    /// contain half the assembly.
    pub n50: u64,
    /// **L50** — the smallest number of contigs whose combined length
    /// reaches half the assembly.
    pub l50: usize,
    /// **N75** — the N50 statistic at the 75 % threshold.
    pub n75: u64,
    /// **N90** — the N50 statistic at the 90 % threshold.
    pub n90: u64,
    /// GC fraction across the whole assembly.
    pub gc_content: f64,
    /// Count of `N` (ambiguous / gap) bases.
    pub n_count: u64,
}

/// Computes the N-`fraction` statistic of a *descending-sorted* length
/// list — the length at which the cumulative sum first reaches
/// `fraction · total`. Returns `0` for an empty list.
fn nx(sorted_desc: &[u64], total: u64, fraction: f64) -> u64 {
    if sorted_desc.is_empty() || total == 0 {
        return 0;
    }
    let threshold = (total as f64 * fraction).ceil() as u64;
    let mut cum = 0u64;
    for &len in sorted_desc {
        cum += len;
        if cum >= threshold {
            return len;
        }
    }
    *sorted_desc.last().unwrap()
}

/// Computes the L-`fraction` statistic — the count of contigs needed
/// to reach `fraction · total`.
fn lx(sorted_desc: &[u64], total: u64, fraction: f64) -> usize {
    if sorted_desc.is_empty() || total == 0 {
        return 0;
    }
    let threshold = (total as f64 * fraction).ceil() as u64;
    let mut cum = 0u64;
    for (i, &len) in sorted_desc.iter().enumerate() {
        cum += len;
        if cum >= threshold {
            return i + 1;
        }
    }
    sorted_desc.len()
}

/// Computes [`AssemblyStats`] from a slice of contig sequences.
///
/// Returns [`GenomicsError::Invalid`] for an empty contig set.
pub fn assembly_stats(contigs: &[&[u8]]) -> Result<AssemblyStats> {
    if contigs.is_empty() {
        return Err(GenomicsError::invalid("contigs", "no contigs supplied"));
    }
    let mut lengths: Vec<u64> = contigs.iter().map(|c| c.len() as u64).collect();
    let total: u64 = lengths.iter().sum();
    let mut gc = 0u64;
    let mut ns = 0u64;
    for c in contigs {
        for &b in *c {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => gc += 1,
                b'A' | b'T' | b'U' => {}
                _ => ns += 1,
            }
        }
    }

    lengths.sort_unstable_by(|a, b| b.cmp(a)); // descending
    let longest = *lengths.first().unwrap();
    let shortest = *lengths.last().unwrap();
    let mean = total as f64 / contigs.len() as f64;
    let gc_content = if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64
    };

    Ok(AssemblyStats {
        n_contigs: contigs.len(),
        total_length: total,
        longest,
        shortest,
        mean_length: mean,
        n50: nx(&lengths, total, 0.5),
        l50: lx(&lengths, total, 0.5),
        n75: nx(&lengths, total, 0.75),
        n90: nx(&lengths, total, 0.90),
        gc_content,
        n_count: ns,
    })
}

/// A k-mer-spectrum summary — the histogram of k-mer multiplicities,
/// the workhorse of genome-size estimation (GenomeScope, kmergenie).
#[derive(Clone, Debug, PartialEq)]
pub struct KmerSpectrum {
    /// The k-mer length.
    pub k: usize,
    /// Number of *distinct* k-mers.
    pub distinct: usize,
    /// Total k-mer observations (with multiplicity).
    pub total: u64,
    /// Multiplicity histogram: `multiplicity -> number of distinct
    /// k-mers seen exactly that many times`.
    pub histogram: Vec<(u32, usize)>,
}

impl KmerSpectrum {
    /// The most common non-trivial multiplicity — an estimate of the
    /// sequencing-coverage peak (k-mers at multiplicity 1 are excluded
    /// as likely errors). Returns `None` when no such peak exists.
    pub fn coverage_peak(&self) -> Option<u32> {
        self.histogram
            .iter()
            .filter(|(mult, _)| *mult > 1)
            .max_by_key(|(_, count)| *count)
            .map(|(mult, _)| *mult)
    }

    /// A rough genome-size estimate: `total k-mer observations /
    /// coverage peak`. Returns `None` without a coverage peak.
    pub fn estimated_genome_size(&self) -> Option<u64> {
        self.coverage_peak()
            .filter(|&peak| peak > 0)
            .map(|peak| self.total / peak as u64)
    }
}

/// Builds a k-mer spectrum from a set of reads.
///
/// k-mers are counted **canonically** — a k-mer and its reverse
/// complement are folded to the lexicographically-smaller of the two,
/// so the strand a read came from does not matter. Returns
/// [`GenomicsError::Invalid`] for `k == 0`.
pub fn kmer_spectrum(reads: &[&[u8]], k: usize) -> Result<KmerSpectrum> {
    if k == 0 {
        return Err(GenomicsError::invalid("k", "k must be positive"));
    }
    let mut counts: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut total = 0u64;
    for read in reads {
        let upper: Vec<u8> = read.iter().map(|b| b.to_ascii_uppercase()).collect();
        if upper.len() < k {
            continue;
        }
        for window in upper.windows(k) {
            // Skip k-mers with ambiguous bases.
            if window
                .iter()
                .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
            {
                continue;
            }
            let canon = canonical_kmer(window);
            *counts.entry(canon).or_insert(0) += 1;
            total += 1;
        }
    }
    let mut hist: HashMap<u32, usize> = HashMap::new();
    for &c in counts.values() {
        *hist.entry(c).or_insert(0) += 1;
    }
    let mut histogram: Vec<(u32, usize)> = hist.into_iter().collect();
    histogram.sort_unstable();
    Ok(KmerSpectrum {
        k,
        distinct: counts.len(),
        total,
        histogram,
    })
}

/// The canonical form of a k-mer — the lexicographically-smaller of the
/// k-mer and its reverse complement.
pub fn canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc: Vec<u8> = kmer
        .iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect();
    let fwd: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();
    if fwd <= rc {
        fwd
    } else {
        rc
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n50_simple() {
        // Contigs 2,3,4,5,6 (total 20, half 10). Descending: 6,5,4,3,2.
        // Cumulative: 6, 11 (>=10) -> N50 = 5, L50 = 2.
        let contigs: Vec<Vec<u8>> = [6, 5, 4, 3, 2].iter().map(|&n| vec![b'A'; n]).collect();
        let refs: Vec<&[u8]> = contigs.iter().map(|c| c.as_slice()).collect();
        let s = assembly_stats(&refs).unwrap();
        assert_eq!(s.total_length, 20);
        assert_eq!(s.n50, 5);
        assert_eq!(s.l50, 2);
        assert_eq!(s.longest, 6);
        assert_eq!(s.shortest, 2);
    }

    #[test]
    fn single_contig() {
        let c = vec![b'A'; 100];
        let s = assembly_stats(&[c.as_slice()]).unwrap();
        assert_eq!(s.n_contigs, 1);
        assert_eq!(s.n50, 100);
        assert_eq!(s.l50, 1);
    }

    #[test]
    fn gc_and_n_counts() {
        let contig = b"GGCCNNAATT"; // 4 GC, 2 N, 4 AT
        let s = assembly_stats(&[contig.as_slice()]).unwrap();
        assert!((s.gc_content - 0.4).abs() < 1e-9);
        assert_eq!(s.n_count, 2);
    }

    #[test]
    fn rejects_empty_assembly() {
        let empty: Vec<&[u8]> = vec![];
        assert!(assembly_stats(&empty).is_err());
    }

    #[test]
    fn n75_n90_ordering() {
        let contigs: Vec<Vec<u8>> = (1..=20).map(|n| vec![b'A'; n * 10]).collect();
        let refs: Vec<&[u8]> = contigs.iter().map(|c| c.as_slice()).collect();
        let s = assembly_stats(&refs).unwrap();
        // N50 >= N75 >= N90 by definition.
        assert!(s.n50 >= s.n75);
        assert!(s.n75 >= s.n90);
    }

    #[test]
    fn canonical_kmer_folds_strands() {
        // "ACG" rc is "CGT"; ACG < CGT so canonical = ACG.
        assert_eq!(canonical_kmer(b"ACG"), b"ACG".to_vec());
        // "CGT" should also canonicalize to "ACG".
        assert_eq!(canonical_kmer(b"CGT"), b"ACG".to_vec());
    }

    #[test]
    fn kmer_spectrum_counts() {
        // Two identical reads -> every k-mer seen twice.
        let reads: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGTACGT"];
        let spec = kmer_spectrum(&reads, 4).unwrap();
        assert!(spec.distinct > 0);
        // All distinct k-mers appear at an even multiplicity.
        let odd = spec.histogram.iter().any(|(mult, _)| mult % 2 == 1);
        assert!(!odd, "identical reads must give even multiplicities");
    }

    #[test]
    fn kmer_spectrum_coverage_peak() {
        // Repeat a read 10x -> coverage peak should be 10. The read is
        // non-repetitive at k = 5 *under canonical folding* (every
        // canonical 5-mer occurs exactly once per read), so each one is
        // seen 10 times total. A tandem-repeat read like ACGTACGTACGT
        // would instead inflate the peak (each 5-mer recurs within the
        // read, and a strand-symmetric k-mer folds onto its rc).
        let reads: Vec<&[u8]> = (0..10).map(|_| b"GCTAAAGACAATTACA".as_slice()).collect();
        let spec = kmer_spectrum(&reads, 5).unwrap();
        assert_eq!(spec.coverage_peak(), Some(10));
    }

    #[test]
    fn kmer_rejects_zero_k() {
        let reads: Vec<&[u8]> = vec![b"ACGT"];
        assert!(kmer_spectrum(&reads, 0).is_err());
    }
}
