//! Read filtering — length, mean quality and low-complexity removal.
//!
//! The third pre-processing pillar: dropping whole reads that fail a
//! quality gate, the way Trimmomatic `MINLEN`, fastp's quality filter
//! and `prinseq` low-complexity filtering work.
//!
//! [`filter_reads`] applies a [`ReadFilter`] to a slice of FASTQ
//! records and returns the survivors plus a [`FilterStats`] breakdown
//! of why each dropped read was dropped.

use valenx_bioseq::io::fastq::FastqRecord;

/// A set of read-level pass/fail gates.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ReadFilter {
    /// Drop reads shorter than this many bases (`0` disables).
    pub min_length: usize,
    /// Drop reads longer than this many bases (`None` disables).
    pub max_length: Option<usize>,
    /// Drop reads whose mean Phred quality is below this (`0.0`
    /// disables).
    pub min_mean_quality: f64,
    /// Drop reads whose `N` fraction exceeds this (`1.0` disables).
    pub max_n_fraction: f64,
    /// Drop reads whose linguistic-complexity score (see
    /// [`complexity_score`]) is below this (`0.0` disables). The score
    /// is in `[0, 1]`; a homopolymer scores near `0`.
    pub min_complexity: f64,
}

impl Default for ReadFilter {
    /// A permissive default: only drop empty reads.
    fn default() -> Self {
        ReadFilter {
            min_length: 1,
            max_length: None,
            min_mean_quality: 0.0,
            max_n_fraction: 1.0,
            min_complexity: 0.0,
        }
    }
}

/// Why a read failed the filter (the first failing reason wins).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum FilterReason {
    /// The read is shorter than `min_length`.
    TooShort,
    /// The read is longer than `max_length`.
    TooLong,
    /// The read's mean quality is below `min_mean_quality`.
    LowQuality,
    /// The read's `N` fraction exceeds `max_n_fraction`.
    TooManyNs,
    /// The read's complexity is below `min_complexity`.
    LowComplexity,
}

/// A per-reason tally of dropped reads.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FilterStats {
    /// Reads that passed every gate.
    pub kept: usize,
    /// Reads dropped for being too short.
    pub too_short: usize,
    /// Reads dropped for being too long.
    pub too_long: usize,
    /// Reads dropped for low mean quality.
    pub low_quality: usize,
    /// Reads dropped for too many `N`s.
    pub too_many_ns: usize,
    /// Reads dropped for low complexity.
    pub low_complexity: usize,
}

impl FilterStats {
    /// Total reads dropped (all reasons summed).
    pub fn dropped(&self) -> usize {
        self.too_short + self.too_long + self.low_quality + self.too_many_ns + self.low_complexity
    }

    fn record(&mut self, reason: FilterReason) {
        match reason {
            FilterReason::TooShort => self.too_short += 1,
            FilterReason::TooLong => self.too_long += 1,
            FilterReason::LowQuality => self.low_quality += 1,
            FilterReason::TooManyNs => self.too_many_ns += 1,
            FilterReason::LowComplexity => self.low_complexity += 1,
        }
    }
}

/// The result of [`filter_reads`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FilterOutput {
    /// The reads that survived.
    pub kept: Vec<FastqRecord>,
    /// The drop breakdown.
    pub stats: FilterStats,
}

/// A linguistic-complexity score in `[0, 1]` based on the diversity of
/// 3-mers in the read.
///
/// The score is the count of *distinct* 3-mers observed divided by the
/// number of distinct 3-mers *possible* given the read length (capped
/// at the 64-element DNA 3-mer alphabet). A homopolymer (`AAAA…`) has
/// one distinct 3-mer and scores near `0`; a random read approaches
/// `1`. Reads shorter than 3 bp score `1.0` (nothing to penalise).
pub fn complexity_score(seq: &[u8]) -> f64 {
    if seq.len() < 3 {
        return 1.0;
    }
    use std::collections::HashSet;
    let mut seen: HashSet<[u8; 3]> = HashSet::new();
    for w in seq.windows(3) {
        seen.insert([
            w[0].to_ascii_uppercase(),
            w[1].to_ascii_uppercase(),
            w[2].to_ascii_uppercase(),
        ]);
    }
    let windows = seq.len() - 2;
    let possible = windows.min(64);
    if possible == 0 {
        1.0
    } else {
        (seen.len() as f64 / possible as f64).min(1.0)
    }
}

/// Tests one read against the filter, returning the first failing
/// [`FilterReason`] or `None` when the read passes.
pub fn evaluate(read: &FastqRecord, filter: &ReadFilter) -> Option<FilterReason> {
    let seq = read.record.seq.as_bytes();
    let len = seq.len();
    if filter.min_length > 0 && len < filter.min_length {
        return Some(FilterReason::TooShort);
    }
    if let Some(max) = filter.max_length {
        if len > max {
            return Some(FilterReason::TooLong);
        }
    }
    if filter.min_mean_quality > 0.0 && !read.quality.is_empty() {
        let mean = read.quality.iter().map(|&q| q as f64).sum::<f64>() / read.quality.len() as f64;
        if mean < filter.min_mean_quality {
            return Some(FilterReason::LowQuality);
        }
    }
    if filter.max_n_fraction < 1.0 && len > 0 {
        let ns = seq
            .iter()
            .filter(|&&b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U'))
            .count();
        if ns as f64 / len as f64 > filter.max_n_fraction {
            return Some(FilterReason::TooManyNs);
        }
    }
    if filter.min_complexity > 0.0 && complexity_score(seq) < filter.min_complexity {
        return Some(FilterReason::LowComplexity);
    }
    None
}

/// Applies a [`ReadFilter`] to a slice of reads.
pub fn filter_reads(reads: &[FastqRecord], filter: &ReadFilter) -> FilterOutput {
    let mut kept = Vec::new();
    let mut stats = FilterStats::default();
    for r in reads {
        match evaluate(r, filter) {
            None => {
                stats.kept += 1;
                kept.push(r.clone());
            }
            Some(reason) => stats.record(reason),
        }
    }
    FilterOutput { kept, stats }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_bioseq::alphabet::SeqKind;
    use valenx_bioseq::record::SeqRecord;
    use valenx_bioseq::seq::Seq;

    fn fq(seq: &str, quals: &[u8]) -> FastqRecord {
        let s = Seq::new(SeqKind::Dna, seq).unwrap();
        FastqRecord {
            record: SeqRecord::new("r", s),
            quality: quals.to_vec(),
        }
    }

    #[test]
    fn complexity_homopolymer_is_low() {
        assert!(complexity_score(b"AAAAAAAAAA") < 0.2);
    }

    #[test]
    fn complexity_diverse_is_high() {
        assert!(complexity_score(b"ACGTACGTGCATGCATCGAT") > 0.5);
    }

    #[test]
    fn length_filter() {
        let f = ReadFilter {
            min_length: 5,
            ..ReadFilter::default()
        };
        assert_eq!(
            evaluate(&fq("ACG", &[40; 3]), &f),
            Some(FilterReason::TooShort)
        );
        assert_eq!(evaluate(&fq("ACGTACGT", &[40; 8]), &f), None);
    }

    #[test]
    fn max_length_filter() {
        let f = ReadFilter {
            max_length: Some(4),
            ..ReadFilter::default()
        };
        assert_eq!(
            evaluate(&fq("ACGTACGT", &[40; 8]), &f),
            Some(FilterReason::TooLong)
        );
    }

    #[test]
    fn quality_filter() {
        let f = ReadFilter {
            min_mean_quality: 20.0,
            ..ReadFilter::default()
        };
        assert_eq!(
            evaluate(&fq("ACGT", &[5, 5, 5, 5]), &f),
            Some(FilterReason::LowQuality)
        );
        assert_eq!(evaluate(&fq("ACGT", &[40, 40, 40, 40]), &f), None);
    }

    #[test]
    fn n_fraction_filter() {
        let f = ReadFilter {
            max_n_fraction: 0.25,
            ..ReadFilter::default()
        };
        // 2 of 4 = 50% N -> fail.
        assert_eq!(
            evaluate(&fq("ANNT", &[40; 4]), &f),
            Some(FilterReason::TooManyNs)
        );
    }

    #[test]
    fn complexity_filter() {
        let f = ReadFilter {
            min_complexity: 0.5,
            ..ReadFilter::default()
        };
        assert_eq!(
            evaluate(&fq("AAAAAAAAAA", &[40; 10]), &f),
            Some(FilterReason::LowComplexity)
        );
    }

    #[test]
    fn filter_reads_tallies_reasons() {
        let f = ReadFilter {
            min_length: 4,
            min_mean_quality: 20.0,
            ..ReadFilter::default()
        };
        let reads = vec![
            fq("ACGTACGT", &[40; 8]), // keep
            fq("AC", &[40; 2]),       // too short
            fq("ACGTAC", &[5; 6]),    // low quality
        ];
        let out = filter_reads(&reads, &f);
        assert_eq!(out.stats.kept, 1);
        assert_eq!(out.stats.too_short, 1);
        assert_eq!(out.stats.low_quality, 1);
        assert_eq!(out.stats.dropped(), 2);
        assert_eq!(out.kept.len(), 1);
    }
}
