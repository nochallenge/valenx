//! FASTA / FASTQ subsampling and downsampling utilities.
//!
//! Down-sampling a read set to a fixed count or a target depth is a
//! routine step before a quick assembly or a coverage-controlled
//! benchmark (`seqtk sample`, `samtools view -s`, `rasusa`). Every
//! routine here is **seeded** and deterministic — the same `(seed,
//! input)` pair always yields the same subsample, on every platform —
//! built on [`crate::util::rng::Rng`].

use crate::error::{GenomicsError, Result};
use crate::util::rng::Rng;
use valenx_bioseq::io::fastq::FastqRecord;
use valenx_bioseq::record::SeqRecord;

/// Draws exactly `n` records uniformly at random *without replacement*
/// from `records`, using `seed` for reproducibility.
///
/// Returns [`GenomicsError::Invalid`] when `n` exceeds the input
/// length. The order of the returned records is the order in which a
/// partial Fisher-Yates shuffle surfaced them (not the input order).
pub fn sample_n<T: Clone>(records: &[T], n: usize, seed: u64) -> Result<Vec<T>> {
    if n > records.len() {
        return Err(GenomicsError::invalid(
            "n",
            format!("cannot sample {n} records from {}", records.len()),
        ));
    }
    let mut idx: Vec<usize> = (0..records.len()).collect();
    let mut rng = Rng::new(seed);
    // Partial Fisher-Yates: surface the first `n` slots.
    for i in 0..n {
        let j = i + rng.below(idx.len() - i);
        idx.swap(i, j);
    }
    Ok(idx[..n].iter().map(|&i| records[i].clone()).collect())
}

/// Keeps each record independently with probability `fraction`.
///
/// This is the `samtools view -s` model: the output size is
/// *approximately* `fraction * len`, not exact. `fraction` must lie in
/// `[0, 1]`.
pub fn sample_fraction<T: Clone>(records: &[T], fraction: f64, seed: u64) -> Result<Vec<T>> {
    if !(0.0..=1.0).contains(&fraction) {
        return Err(GenomicsError::invalid(
            "fraction",
            "fraction must be in [0, 1]",
        ));
    }
    let mut rng = Rng::new(seed);
    Ok(records
        .iter()
        .filter(|_| rng.chance(fraction))
        .cloned()
        .collect())
}

/// Down-samples FASTQ reads to a target mean depth over a genome of
/// `genome_size` bases.
///
/// The retained read count is `genome_size * target_depth / mean_read_length`,
/// capped at the input size; reads are then drawn with [`sample_n`].
/// Returns [`GenomicsError::Invalid`] for a zero `genome_size` or a
/// non-positive `target_depth`.
pub fn downsample_to_depth(
    reads: &[FastqRecord],
    genome_size: u64,
    target_depth: f64,
    seed: u64,
) -> Result<Vec<FastqRecord>> {
    if genome_size == 0 {
        return Err(GenomicsError::invalid("genome_size", "must be positive"));
    }
    if target_depth <= 0.0 {
        return Err(GenomicsError::invalid("target_depth", "must be positive"));
    }
    if reads.is_empty() {
        return Ok(Vec::new());
    }
    let total_bases: u64 = reads.iter().map(|r| r.len() as u64).sum();
    let mean_len = total_bases as f64 / reads.len() as f64;
    if mean_len <= 0.0 {
        return Ok(Vec::new());
    }
    let wanted = ((genome_size as f64 * target_depth) / mean_len).round() as usize;
    let n = wanted.min(reads.len());
    sample_n(reads, n, seed)
}

/// Subsamples plain FASTA records (id + sequence, no quality) — a thin
/// wrapper over [`sample_n`] for [`SeqRecord`].
pub fn sample_fasta(records: &[SeqRecord], n: usize, seed: u64) -> Result<Vec<SeqRecord>> {
    sample_n(records, n, seed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_bioseq::alphabet::SeqKind;
    use valenx_bioseq::seq::Seq;

    fn fq(id: &str, len: usize) -> FastqRecord {
        let s = Seq::new(SeqKind::Dna, "A".repeat(len)).unwrap();
        FastqRecord {
            record: SeqRecord::new(id, s),
            quality: vec![40u8; len],
        }
    }

    #[test]
    fn sample_n_exact_count() {
        let items: Vec<i32> = (0..100).collect();
        let s = sample_n(&items, 10, 7).unwrap();
        assert_eq!(s.len(), 10);
        // All drawn indices are distinct.
        let mut sorted = s.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), 10);
    }

    #[test]
    fn sample_n_is_deterministic() {
        let items: Vec<i32> = (0..100).collect();
        let a = sample_n(&items, 20, 42).unwrap();
        let b = sample_n(&items, 20, 42).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn sample_n_rejects_oversize() {
        let items: Vec<i32> = (0..5).collect();
        assert!(sample_n(&items, 10, 1).is_err());
    }

    #[test]
    fn sample_fraction_approximate() {
        let items: Vec<i32> = (0..1000).collect();
        let s = sample_fraction(&items, 0.3, 99).unwrap();
        // Roughly 300, allow a wide tolerance for a small sample.
        assert!((200..400).contains(&s.len()), "got {}", s.len());
    }

    #[test]
    fn sample_fraction_bounds() {
        let items: Vec<i32> = (0..10).collect();
        assert_eq!(sample_fraction(&items, 0.0, 1).unwrap().len(), 0);
        assert_eq!(sample_fraction(&items, 1.0, 1).unwrap().len(), 10);
        assert!(sample_fraction(&items, 1.5, 1).is_err());
    }

    #[test]
    fn downsample_to_depth_caps_at_input() {
        // 100 reads of length 100 = 10_000 bases. Genome 1000, depth 5
        // wants 50 reads.
        let reads: Vec<FastqRecord> = (0..100).map(|i| fq(&format!("r{i}"), 100)).collect();
        let s = downsample_to_depth(&reads, 1000, 5.0, 3).unwrap();
        assert_eq!(s.len(), 50);
    }

    #[test]
    fn downsample_depth_above_available_returns_all() {
        let reads: Vec<FastqRecord> = (0..10).map(|i| fq(&format!("r{i}"), 100)).collect();
        // Asking for depth 1000 over a tiny genome wants more reads
        // than exist -> capped.
        let s = downsample_to_depth(&reads, 100_000, 1000.0, 3).unwrap();
        assert_eq!(s.len(), 10);
    }

    #[test]
    fn downsample_validation() {
        let reads = vec![fq("r", 100)];
        assert!(downsample_to_depth(&reads, 0, 5.0, 1).is_err());
        assert!(downsample_to_depth(&reads, 100, 0.0, 1).is_err());
    }
}
