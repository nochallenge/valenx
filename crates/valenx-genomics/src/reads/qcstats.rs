//! FASTQ quality handling and FastQC-class read statistics.
//!
//! This module turns a set of
//! [`valenx_bioseq::io::fastq::FastqRecord`] values into a
//! [`FastqcReport`] — the per-base and per-read summary FastQC popularised:
//! per-position quality quartiles, the per-read mean-quality
//! distribution, per-position base composition, the GC distribution
//! and a length histogram.
//!
//! It works on `valenx-bioseq`'s already-decoded Phred scores, so the
//! quality-codec offset is handled upstream.

use valenx_bioseq::io::fastq::FastqRecord;

/// A simple five-number summary of a numeric distribution.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Quartiles {
    /// Minimum value.
    pub min: f64,
    /// First quartile (25th percentile).
    pub q1: f64,
    /// Median (50th percentile).
    pub median: f64,
    /// Third quartile (75th percentile).
    pub q3: f64,
    /// Maximum value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
}

impl Quartiles {
    /// Computes the five-number summary of `values`. Returns an
    /// all-zero summary for an empty slice.
    pub fn of(values: &[f64]) -> Self {
        if values.is_empty() {
            return Quartiles {
                min: 0.0,
                q1: 0.0,
                median: 0.0,
                q3: 0.0,
                max: 0.0,
                mean: 0.0,
            };
        }
        let mut v = values.to_vec();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let pct = |p: f64| -> f64 {
            // Nearest-rank percentile.
            let rank = (p * (v.len() as f64 - 1.0)).round() as usize;
            v[rank.min(v.len() - 1)]
        };
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        Quartiles {
            min: v[0],
            q1: pct(0.25),
            median: pct(0.50),
            q3: pct(0.75),
            max: v[v.len() - 1],
            mean,
        }
    }
}

/// Per-position base-composition counts (the FastQC "per-base sequence
/// content" plot).
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct BaseComposition {
    /// Count of `A`.
    pub a: u64,
    /// Count of `C`.
    pub c: u64,
    /// Count of `G`.
    pub g: u64,
    /// Count of `T`.
    pub t: u64,
    /// Count of `N` (and any other ambiguity code).
    pub n: u64,
}

impl BaseComposition {
    /// Total observations at this position.
    pub fn total(&self) -> u64 {
        self.a + self.c + self.g + self.t + self.n
    }

    /// GC fraction at this position (`0.0` when no observations).
    pub fn gc_fraction(&self) -> f64 {
        let t = self.total();
        if t == 0 {
            0.0
        } else {
            (self.g + self.c) as f64 / t as f64
        }
    }

    fn observe(&mut self, base: u8) {
        match base.to_ascii_uppercase() {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' | b'U' => self.t += 1,
            _ => self.n += 1,
        }
    }
}

/// The full FastQC-class report for a FASTQ dataset.
#[derive(Clone, Debug, PartialEq)]
pub struct FastqcReport {
    /// Number of reads analysed.
    pub n_reads: usize,
    /// Total bases across all reads.
    pub total_bases: u64,
    /// Shortest and longest read length.
    pub min_len: usize,
    /// Longest read length.
    pub max_len: usize,
    /// Overall mean read length.
    pub mean_len: f64,
    /// Overall GC fraction across all bases.
    pub gc_content: f64,
    /// Per-position quality quartiles; index `i` is read position `i`.
    pub per_position_quality: Vec<Quartiles>,
    /// Per-position base composition; index `i` is read position `i`.
    pub per_position_composition: Vec<BaseComposition>,
    /// Per-read mean-quality summary (the "per-sequence quality score"
    /// distribution).
    pub per_read_quality: Quartiles,
    /// Histogram of read lengths: `length -> count`, sorted by length.
    pub length_histogram: Vec<(usize, usize)>,
    /// Fraction of all bases that are `N`.
    pub n_fraction: f64,
}

impl FastqcReport {
    /// Builds the report from a slice of FASTQ records.
    pub fn analyze(records: &[FastqRecord]) -> Self {
        let n_reads = records.len();
        if n_reads == 0 {
            return FastqcReport {
                n_reads: 0,
                total_bases: 0,
                min_len: 0,
                max_len: 0,
                mean_len: 0.0,
                gc_content: 0.0,
                per_position_quality: Vec::new(),
                per_position_composition: Vec::new(),
                per_read_quality: Quartiles::of(&[]),
                length_histogram: Vec::new(),
                n_fraction: 0.0,
            };
        }

        let max_len = records.iter().map(|r| r.len()).max().unwrap_or(0);
        let min_len = records.iter().map(|r| r.len()).min().unwrap_or(0);

        // Per-position accumulators.
        let mut pos_quals: Vec<Vec<f64>> = vec![Vec::new(); max_len];
        let mut pos_comp: Vec<BaseComposition> = vec![BaseComposition::default(); max_len];
        let mut read_means: Vec<f64> = Vec::with_capacity(n_reads);
        let mut total_bases: u64 = 0;
        let mut gc_bases: u64 = 0;
        let mut n_bases: u64 = 0;
        let mut len_counts: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();

        for rec in records {
            let seq = rec.record.seq.as_bytes();
            *len_counts.entry(seq.len()).or_insert(0) += 1;
            total_bases += seq.len() as u64;
            for (i, &b) in seq.iter().enumerate() {
                pos_comp[i].observe(b);
                match b.to_ascii_uppercase() {
                    b'G' | b'C' => gc_bases += 1,
                    b'A' | b'T' | b'U' => {}
                    _ => n_bases += 1,
                }
            }
            for (i, &q) in rec.quality.iter().enumerate() {
                if i < max_len {
                    pos_quals[i].push(q as f64);
                }
            }
            if !rec.quality.is_empty() {
                read_means.push(
                    rec.quality.iter().map(|&q| q as f64).sum::<f64>() / rec.quality.len() as f64,
                );
            }
        }

        let per_position_quality: Vec<Quartiles> =
            pos_quals.iter().map(|v| Quartiles::of(v)).collect();
        let mean_len = total_bases as f64 / n_reads as f64;
        let gc_content = if total_bases == 0 {
            0.0
        } else {
            gc_bases as f64 / total_bases as f64
        };
        let n_fraction = if total_bases == 0 {
            0.0
        } else {
            n_bases as f64 / total_bases as f64
        };

        FastqcReport {
            n_reads,
            total_bases,
            min_len,
            max_len,
            mean_len,
            gc_content,
            per_position_quality,
            per_position_composition: pos_comp,
            per_read_quality: Quartiles::of(&read_means),
            length_histogram: len_counts.into_iter().collect(),
            n_fraction,
        }
    }

    /// A simple pass / warn / fail verdict on the per-base quality, in
    /// the spirit of FastQC's module flags: **fail** when any position's
    /// lower quartile drops below 5 or the median below 20; **warn**
    /// when the lower quartile drops below 10 or the median below 25;
    /// **pass** otherwise.
    pub fn per_base_quality_verdict(&self) -> QcVerdict {
        let mut verdict = QcVerdict::Pass;
        for q in &self.per_position_quality {
            if q.q1 < 5.0 || q.median < 20.0 {
                return QcVerdict::Fail;
            }
            if q.q1 < 10.0 || q.median < 25.0 {
                verdict = QcVerdict::Warn;
            }
        }
        verdict
    }

    /// Renders a compact multi-line text summary.
    pub fn summary_text(&self) -> String {
        format!(
            "reads={} total_bases={} len={}..{} mean_len={:.1} GC={:.1}% N={:.2}% \
             per_read_quality(median={:.1}) verdict={:?}",
            self.n_reads,
            self.total_bases,
            self.min_len,
            self.max_len,
            self.mean_len,
            self.gc_content * 100.0,
            self.n_fraction * 100.0,
            self.per_read_quality.median,
            self.per_base_quality_verdict(),
        )
    }
}

/// A FastQC-style module verdict.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum QcVerdict {
    /// The module passed.
    Pass,
    /// The module raised a warning.
    Warn,
    /// The module failed.
    Fail,
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_bioseq::alphabet::SeqKind;
    use valenx_bioseq::record::SeqRecord;
    use valenx_bioseq::seq::Seq;

    fn fq(id: &str, seq: &str, quals: &[u8]) -> FastqRecord {
        let s = Seq::new(SeqKind::Dna, seq).unwrap();
        FastqRecord {
            record: SeqRecord::new(id, s),
            quality: quals.to_vec(),
        }
    }

    #[test]
    fn quartiles_basic() {
        let q = Quartiles::of(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert_eq!(q.min, 1.0);
        assert_eq!(q.max, 5.0);
        assert_eq!(q.median, 3.0);
        assert_eq!(q.mean, 3.0);
    }

    #[test]
    fn empty_dataset_is_safe() {
        let r = FastqcReport::analyze(&[]);
        assert_eq!(r.n_reads, 0);
        assert_eq!(r.gc_content, 0.0);
    }

    #[test]
    fn counts_reads_and_bases() {
        let recs = vec![
            fq("r1", "ACGT", &[40, 40, 40, 40]),
            fq("r2", "GGCC", &[30, 30, 30, 30]),
        ];
        let r = FastqcReport::analyze(&recs);
        assert_eq!(r.n_reads, 2);
        assert_eq!(r.total_bases, 8);
        assert_eq!(r.min_len, 4);
        assert_eq!(r.max_len, 4);
    }

    #[test]
    fn gc_content_computed() {
        // ACGT = 50% GC, GGCC = 100% GC -> overall 75%.
        let recs = vec![
            fq("r1", "ACGT", &[40, 40, 40, 40]),
            fq("r2", "GGCC", &[40, 40, 40, 40]),
        ];
        let r = FastqcReport::analyze(&recs);
        assert!((r.gc_content - 0.75).abs() < 1e-9);
    }

    #[test]
    fn per_position_quality_tracked() {
        let recs = vec![
            fq("r1", "ACGT", &[40, 30, 20, 10]),
            fq("r2", "ACGT", &[40, 30, 20, 10]),
        ];
        let r = FastqcReport::analyze(&recs);
        assert_eq!(r.per_position_quality.len(), 4);
        assert_eq!(r.per_position_quality[0].median, 40.0);
        assert_eq!(r.per_position_quality[3].median, 10.0);
    }

    #[test]
    fn verdict_flags_low_quality() {
        // All positions quality 2 -> fail.
        let recs = vec![fq("r1", "ACGT", &[2, 2, 2, 2])];
        let r = FastqcReport::analyze(&recs);
        assert_eq!(r.per_base_quality_verdict(), QcVerdict::Fail);
        // All positions quality 40 -> pass.
        let recs = vec![fq("r2", "ACGT", &[40, 40, 40, 40])];
        let r = FastqcReport::analyze(&recs);
        assert_eq!(r.per_base_quality_verdict(), QcVerdict::Pass);
    }

    #[test]
    fn length_histogram_built() {
        let recs = vec![
            fq("r1", "ACGT", &[40; 4]),
            fq("r2", "ACG", &[40; 3]),
            fq("r3", "ACGT", &[40; 4]),
        ];
        let r = FastqcReport::analyze(&recs);
        assert_eq!(r.length_histogram, vec![(3, 1), (4, 2)]);
    }

    #[test]
    fn n_fraction_counts_ambiguity() {
        let recs = vec![fq("r1", "ACGN", &[40; 4])];
        let r = FastqcReport::analyze(&recs);
        assert!((r.n_fraction - 0.25).abs() < 1e-9);
    }
}
