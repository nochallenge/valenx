//! Adapter trimming and quality trimming.
//!
//! Two of the three pillars of read pre-processing (Trimmomatic /
//! cutadapt / fastp):
//!
//! - **Adapter trimming** ([`trim_adapter`]) — removes a known adapter
//!   sequence from the 3′ end of a read, including a partial adapter
//!   that runs off the read end, via a mismatch-tolerant overlap scan.
//! - **Quality trimming** ([`trim_quality`]) — removes low-quality
//!   bases with a leading/trailing trim and a sliding-window trim, the
//!   Trimmomatic `LEADING` / `TRAILING` / `SLIDINGWINDOW` steps.
//!
//! Both operate on `valenx-bioseq`'s
//! [`valenx_bioseq::io::fastq::FastqRecord`] and return a
//! new trimmed record so the input is never mutated.

use crate::error::{GenomicsError, Result};
use valenx_bioseq::io::fastq::FastqRecord;
use valenx_bioseq::record::SeqRecord;
use valenx_bioseq::seq::Seq;

/// A handful of common Illumina adapter sequences (the 3′ adapter is
/// what most trimmers strip).
pub mod adapters {
    /// The TruSeq / generic Illumina universal 3′ adapter prefix.
    pub const ILLUMINA_TRUSEQ: &str = "AGATCGGAAGAGC";
    /// The Nextera transposase 3′ adapter prefix.
    pub const NEXTERA: &str = "CTGTCTCTTATACACATCT";
    /// The small-RNA 3′ adapter.
    pub const SMALL_RNA: &str = "TGGAATTCTCGG";
}

/// Outcome of an adapter-trimming call.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterTrimResult {
    /// The trimmed read.
    pub record: FastqRecord,
    /// Number of bases removed from the 3′ end (`0` when no adapter
    /// matched).
    pub trimmed: usize,
    /// `true` when an adapter (full or partial) was found and removed.
    pub adapter_found: bool,
}

/// Trims a 3′ adapter from a read.
///
/// The scan tries every alignment offset `o` of the adapter against the
/// read suffix starting at `o`. A match needs at least `min_overlap`
/// overlapping bases and an identity of at least `1 - max_error_rate`.
/// The earliest qualifying offset is taken (so the *most* of the read
/// is trimmed — standard cutadapt behaviour for a 3′ adapter). A
/// partial adapter running off the read end is matched too.
///
/// Returns [`GenomicsError::Invalid`] for an empty adapter, a
/// non-positive `min_overlap`, or `max_error_rate` outside `[0, 1]`.
pub fn trim_adapter(
    read: &FastqRecord,
    adapter: &[u8],
    min_overlap: usize,
    max_error_rate: f64,
) -> Result<AdapterTrimResult> {
    if adapter.is_empty() {
        return Err(GenomicsError::invalid(
            "adapter",
            "adapter must be non-empty",
        ));
    }
    if min_overlap == 0 {
        return Err(GenomicsError::invalid(
            "min_overlap",
            "min_overlap must be positive",
        ));
    }
    if !(0.0..=1.0).contains(&max_error_rate) {
        return Err(GenomicsError::invalid(
            "max_error_rate",
            "max_error_rate must be in [0, 1]",
        ));
    }

    let seq = read.record.seq.as_bytes();
    let n = seq.len();
    let adapter: Vec<u8> = adapter.iter().map(|b| b.to_ascii_uppercase()).collect();

    let mut best_cut: Option<usize> = None;
    for offset in 0..n {
        // Overlap length = min(adapter, read suffix from offset).
        let overlap = adapter.len().min(n - offset);
        if overlap < min_overlap {
            // Further offsets only shrink the overlap.
            break;
        }
        let mut mismatches = 0usize;
        for k in 0..overlap {
            if seq[offset + k].to_ascii_uppercase() != adapter[k] {
                mismatches += 1;
            }
        }
        let allowed = (overlap as f64 * max_error_rate).floor() as usize;
        if mismatches <= allowed {
            best_cut = Some(offset);
            break; // earliest offset = longest trim
        }
    }

    match best_cut {
        Some(cut) => {
            let kept = &seq[..cut];
            let new_seq = Seq::new(read.record.seq.kind(), kept)
                .map_err(|e| GenomicsError::invalid("read", e.to_string()))?;
            let mut new_rec = SeqRecord::new(read.record.id.clone(), new_seq);
            new_rec.description = read.record.description.clone();
            Ok(AdapterTrimResult {
                record: FastqRecord {
                    record: new_rec,
                    quality: read.quality[..cut].to_vec(),
                },
                trimmed: n - cut,
                adapter_found: true,
            })
        }
        None => Ok(AdapterTrimResult {
            record: read.clone(),
            trimmed: 0,
            adapter_found: false,
        }),
    }
}

/// Sliding-window and leading/trailing quality-trim parameters.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct QualityTrimParams {
    /// Trim leading bases below this Phred quality.
    pub leading: u8,
    /// Trim trailing bases below this Phred quality.
    pub trailing: u8,
    /// Sliding-window size (`0` disables the window step).
    pub window: usize,
    /// Cut once a window's mean quality drops below this threshold.
    pub window_quality: f64,
}

impl Default for QualityTrimParams {
    /// Trimmomatic-like defaults: `LEADING:3 TRAILING:3
    /// SLIDINGWINDOW:4:15`.
    fn default() -> Self {
        QualityTrimParams {
            leading: 3,
            trailing: 3,
            window: 4,
            window_quality: 15.0,
        }
    }
}

/// Outcome of a quality-trimming call.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QualityTrimResult {
    /// The trimmed read (may be empty if every base failed).
    pub record: FastqRecord,
    /// Bases removed from the 5′ end.
    pub trimmed_start: usize,
    /// Bases removed from the 3′ end.
    pub trimmed_end: usize,
}

/// Quality-trims a read.
///
/// The steps run in order: (1) trim leading bases below
/// [`leading`](QualityTrimParams::leading); (2) trim trailing bases
/// below [`trailing`](QualityTrimParams::trailing); (3) if
/// [`window`](QualityTrimParams::window) is non-zero, scan 5′→3′ and as
/// soon as a window's mean quality drops below
/// [`window_quality`](QualityTrimParams::window_quality), cut the read
/// at the start of that window.
pub fn trim_quality(read: &FastqRecord, params: QualityTrimParams) -> Result<QualityTrimResult> {
    let seq = read.record.seq.as_bytes();
    let qual = &read.quality;
    if seq.len() != qual.len() {
        return Err(GenomicsError::invalid_record(
            "fastq",
            "SEQ and QUAL lengths disagree",
        ));
    }
    let n = seq.len();

    // (1) leading trim
    let mut start = 0usize;
    while start < n && qual[start] < params.leading {
        start += 1;
    }
    // (2) trailing trim
    let mut end = n;
    while end > start && qual[end - 1] < params.trailing {
        end -= 1;
    }
    // (3) sliding window
    if params.window > 0 && end > start {
        let mut i = start;
        while i + params.window <= end {
            let mean = qual[i..i + params.window]
                .iter()
                .map(|&q| q as f64)
                .sum::<f64>()
                / params.window as f64;
            if mean < params.window_quality {
                end = i;
                break;
            }
            i += 1;
        }
    }
    let end = end.max(start);

    let kept = &seq[start..end];
    let new_seq = Seq::new(read.record.seq.kind(), kept)
        .map_err(|e| GenomicsError::invalid("read", e.to_string()))?;
    let mut new_rec = SeqRecord::new(read.record.id.clone(), new_seq);
    new_rec.description = read.record.description.clone();
    Ok(QualityTrimResult {
        record: FastqRecord {
            record: new_rec,
            quality: qual[start..end].to_vec(),
        },
        trimmed_start: start,
        trimmed_end: n - end,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_bioseq::alphabet::SeqKind;
    use valenx_bioseq::record::SeqRecord;

    fn fq(seq: &str, quals: &[u8]) -> FastqRecord {
        let s = Seq::new(SeqKind::Dna, seq).unwrap();
        FastqRecord {
            record: SeqRecord::new("r", s),
            quality: quals.to_vec(),
        }
    }

    #[test]
    fn trims_full_3prime_adapter() {
        // Read = AAAAAA + adapter.
        let adapter = adapters::ILLUMINA_TRUSEQ;
        let seq = format!("AAAAAA{adapter}");
        let q = vec![40u8; seq.len()];
        let read = fq(&seq, &q);
        let r = trim_adapter(&read, adapter.as_bytes(), 3, 0.1).unwrap();
        assert!(r.adapter_found);
        assert_eq!(r.record.record.seq.as_str(), "AAAAAA");
        assert_eq!(r.trimmed, adapter.len());
    }

    #[test]
    fn trims_partial_adapter_at_read_end() {
        // Only the first 5 bases of the adapter are present.
        let adapter = b"AGATCGGAAGAGC";
        let seq = "TTTTTTTTAGATC"; // 8 T + 5-base adapter prefix
        let q = vec![40u8; seq.len()];
        let read = fq(seq, &q);
        let r = trim_adapter(&read, adapter, 4, 0.1).unwrap();
        assert!(r.adapter_found);
        assert_eq!(r.record.record.seq.as_str(), "TTTTTTTT");
    }

    #[test]
    fn tolerates_a_mismatch() {
        // Adapter with one substituted base inside.
        let adapter = b"AGATCGGAAGAGC";
        let seq = "GGGGGGAGATCGGTAGAGC"; // one mismatch (A->T) in adapter
        let q = vec![40u8; seq.len()];
        let read = fq(seq, &q);
        let r = trim_adapter(&read, adapter, 6, 0.15).unwrap();
        assert!(r.adapter_found);
        assert_eq!(r.record.record.seq.as_str(), "GGGGGG");
    }

    #[test]
    fn no_adapter_leaves_read_intact() {
        let read = fq("ACGTACGTACGT", &[40u8; 12]);
        let r = trim_adapter(&read, b"TTTTTTTTTTTT", 5, 0.1).unwrap();
        assert!(!r.adapter_found);
        assert_eq!(r.record.record.seq.as_str(), "ACGTACGTACGT");
    }

    #[test]
    fn adapter_validation() {
        let read = fq("ACGT", &[40; 4]);
        assert!(trim_adapter(&read, b"", 3, 0.1).is_err());
        assert!(trim_adapter(&read, b"AC", 0, 0.1).is_err());
        assert!(trim_adapter(&read, b"AC", 3, 1.5).is_err());
    }

    #[test]
    fn quality_leading_trailing_trim() {
        // Low quality at both ends.
        let read = fq("ACGTAC", &[2, 2, 40, 40, 2, 2]);
        let p = QualityTrimParams {
            leading: 10,
            trailing: 10,
            window: 0,
            window_quality: 0.0,
        };
        let r = trim_quality(&read, p).unwrap();
        assert_eq!(r.record.record.seq.as_str(), "GT");
        assert_eq!(r.trimmed_start, 2);
        assert_eq!(r.trimmed_end, 2);
    }

    #[test]
    fn sliding_window_cuts_at_drop() {
        // Quality 40,40,40,40,2,2,2,2 with a 4-base window, threshold
        // 20. The first window whose mean falls below 20 is the one
        // *starting at index 3* (40,2,2,2 -> mean 11.5; the window
        // starting at index 2 is 40,40,2,2 -> mean 21, still passing).
        // The read is cut at the start of that failing window, the
        // standard Trimmomatic SLIDINGWINDOW behaviour, leaving "ACG".
        let read = fq("ACGTACGT", &[40, 40, 40, 40, 2, 2, 2, 2]);
        let p = QualityTrimParams {
            leading: 0,
            trailing: 0,
            window: 4,
            window_quality: 20.0,
        };
        let r = trim_quality(&read, p).unwrap();
        assert_eq!(r.record.record.seq.as_str(), "ACG");
    }

    #[test]
    fn all_low_quality_yields_empty_read() {
        let read = fq("ACGT", &[2, 2, 2, 2]);
        let p = QualityTrimParams {
            leading: 10,
            trailing: 10,
            window: 0,
            window_quality: 0.0,
        };
        let r = trim_quality(&read, p).unwrap();
        assert!(r.record.is_empty());
    }
}
