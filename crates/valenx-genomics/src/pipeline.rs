//! Top-level batch pipeline helpers and the bundled [`GenomicsReport`].
//!
//! The other modules each cover one stage of an NGS workflow. This
//! module wires them together: a few batch convenience functions and a
//! single [`GenomicsReport`] that bundles read QC, mapping statistics
//! and variant counts — the one-glance summary an analyst wants after
//! a resequencing run.
//!
//! Nothing here adds new science; it is orchestration over
//! [`crate::reads`], [`crate::format`] and [`crate::variant`].

use crate::error::Result;
use crate::format::pileup::{build_pileup, Reference};
use crate::format::sam::{SamFile, SamRecord};
use crate::format::vcf::VcfFile;
use crate::reads::coverage::{compute_depth, ContigLengths};
use crate::reads::qcstats::FastqcReport;
use crate::variant::call::{call_variants, CallParams};
use crate::variant::stats::{vcf_stats, VcfStats};
use valenx_bioseq::io::fastq::FastqRecord;

/// Read-set quality-control summary — the headline numbers of a
/// [`FastqcReport`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadQcSummary {
    /// Number of reads.
    pub n_reads: usize,
    /// Total bases.
    pub total_bases: u64,
    /// Mean read length.
    pub mean_length: f64,
    /// Overall GC fraction.
    pub gc_content: f64,
    /// Median per-read mean quality.
    pub median_read_quality: f64,
    /// Fraction of bases that are `N`.
    pub n_fraction: f64,
}

impl ReadQcSummary {
    /// Distils a [`FastqcReport`] into the headline summary.
    pub fn from_report(report: &FastqcReport) -> Self {
        ReadQcSummary {
            n_reads: report.n_reads,
            total_bases: report.total_bases,
            mean_length: report.mean_len,
            gc_content: report.gc_content,
            median_read_quality: report.per_read_quality.median,
            n_fraction: report.n_fraction,
        }
    }
}

/// Alignment / mapping summary over a SAM record set — the numbers
/// `samtools flagstat` reports.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MappingSummary {
    /// Total alignment records.
    pub total_records: usize,
    /// Mapped records.
    pub mapped: usize,
    /// Unmapped records.
    pub unmapped: usize,
    /// Records flagged as duplicates.
    pub duplicates: usize,
    /// Records flagged secondary or supplementary.
    pub secondary_supplementary: usize,
    /// Paired records.
    pub paired: usize,
    /// Properly-paired records.
    pub properly_paired: usize,
    /// Mean MAPQ over the mapped records.
    pub mean_mapq: f64,
    /// Mapping rate `mapped / total`.
    pub mapping_rate: f64,
}

impl MappingSummary {
    /// Computes the mapping summary from a slice of SAM records.
    pub fn from_records(records: &[SamRecord]) -> Self {
        let total = records.len();
        let mut mapped = 0usize;
        let mut unmapped = 0usize;
        let mut duplicates = 0usize;
        let mut sec_supp = 0usize;
        let mut paired = 0usize;
        let mut proper = 0usize;
        let mut mapq_sum = 0u64;
        for r in records {
            if r.is_unmapped() {
                unmapped += 1;
            } else {
                mapped += 1;
                mapq_sum += r.mapq as u64;
            }
            if r.flags.is_duplicate() {
                duplicates += 1;
            }
            if r.flags.is_secondary_or_supplementary() {
                sec_supp += 1;
            }
            if r.flags.is_paired() {
                paired += 1;
            }
            if r.flags.has(crate::format::sam::SamFlags::PROPER_PAIR) {
                proper += 1;
            }
        }
        MappingSummary {
            total_records: total,
            mapped,
            unmapped,
            duplicates,
            secondary_supplementary: sec_supp,
            paired,
            properly_paired: proper,
            mean_mapq: if mapped == 0 {
                0.0
            } else {
                mapq_sum as f64 / mapped as f64
            },
            mapping_rate: if total == 0 {
                0.0
            } else {
                mapped as f64 / total as f64
            },
        }
    }
}

/// The bundled end-of-run report: read QC, mapping stats, variant
/// counts and coverage in one struct.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomicsReport {
    /// Read-set QC (`None` when no reads were supplied).
    pub read_qc: Option<ReadQcSummary>,
    /// Mapping statistics (`None` when no alignments were supplied).
    pub mapping: Option<MappingSummary>,
    /// Variant-set statistics (`None` when no variants were supplied).
    pub variants: Option<VcfStats>,
    /// Mean genome-wide depth across the alignments (`0.0` when no
    /// alignments).
    pub mean_depth: f64,
    /// Genome-wide breadth of coverage at depth >= 1.
    pub breadth: f64,
}

impl GenomicsReport {
    /// An empty report.
    pub fn new() -> Self {
        GenomicsReport::default()
    }

    /// A compact, human-readable multi-line summary.
    pub fn summary_text(&self) -> String {
        let mut lines = Vec::new();
        match &self.read_qc {
            Some(q) => lines.push(format!(
                "reads: n={} bases={} mean_len={:.1} GC={:.1}% median_Q={:.1}",
                q.n_reads,
                q.total_bases,
                q.mean_length,
                q.gc_content * 100.0,
                q.median_read_quality,
            )),
            None => lines.push("reads: (none)".to_string()),
        }
        match &self.mapping {
            Some(m) => lines.push(format!(
                "mapping: {}/{} mapped ({:.1}%) dup={} mean_MAPQ={:.1}",
                m.mapped,
                m.total_records,
                m.mapping_rate * 100.0,
                m.duplicates,
                m.mean_mapq,
            )),
            None => lines.push("mapping: (none)".to_string()),
        }
        lines.push(format!(
            "coverage: mean_depth={:.1}x breadth={:.1}%",
            self.mean_depth,
            self.breadth * 100.0,
        ));
        match &self.variants {
            Some(v) => lines.push(format!(
                "variants: {} total ({} SNV, {} indel) Ts/Tv={:.2} pass={}",
                v.total,
                v.snvs,
                v.indels,
                v.ts_tv_ratio(),
                v.passing,
            )),
            None => lines.push("variants: (none)".to_string()),
        }
        lines.join("\n")
    }
}

/// Builds a [`GenomicsReport`] from the raw inputs of a run.
///
/// Any input may be empty; the corresponding report field is then
/// `None` / `0.0`. `contig_lengths` lets the depth computation extend
/// coverage to declared contig lengths.
pub fn build_report(
    reads: &[FastqRecord],
    alignments: &[SamRecord],
    variants: Option<&VcfFile>,
    contig_lengths: &ContigLengths,
) -> GenomicsReport {
    let read_qc = if reads.is_empty() {
        None
    } else {
        Some(ReadQcSummary::from_report(&FastqcReport::analyze(reads)))
    };
    let mapping = if alignments.is_empty() {
        None
    } else {
        Some(MappingSummary::from_records(alignments))
    };
    let (mean_depth, breadth) = if alignments.is_empty() {
        (0.0, 0.0)
    } else {
        let profile = compute_depth(alignments, 0, contig_lengths);
        (profile.genome_mean_depth(), profile.genome_breadth())
    };
    let variants = variants.map(vcf_stats);
    GenomicsReport {
        read_qc,
        mapping,
        variants,
        mean_depth,
        breadth,
    }
}

/// A convenience batch pipeline: take a SAM file and a reference, build
/// the pileup, call variants, and return both the calls and a
/// [`GenomicsReport`].
///
/// This is the "align → call → summarise" backbone an analyst runs.
/// The reads-QC field of the report is left `None` (this entry point
/// starts from alignments, not raw FASTQ).
pub fn call_and_report(
    sam: &SamFile,
    reference: &Reference,
    params: &CallParams,
) -> Result<(Vec<crate::variant::call::Variant>, GenomicsReport)> {
    let columns = build_pileup(&sam.records, reference, 0)?;
    let variants = call_variants(&columns, params)?;

    // Turn the calls into a small VCF so the report can summarise them.
    let mut vcf = VcfFile::new();
    for v in &variants {
        let mut rec = crate::format::vcf::VcfRecord::snv(&v.chrom, v.pos, &v.reference, &v.alt);
        rec.qual = Some(v.qual);
        rec.filter = vec!["PASS".to_string()];
        vcf.records.push(rec);
    }

    let report = build_report(&[], &sam.records, Some(&vcf), &ContigLengths::new());
    Ok((variants, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::sam::{Cigar, SamFlags};
    use valenx_bioseq::alphabet::SeqKind;
    use valenx_bioseq::record::SeqRecord;
    use valenx_bioseq::seq::Seq;

    fn fq(seq: &str) -> FastqRecord {
        let s = Seq::new(SeqKind::Dna, seq).unwrap();
        FastqRecord {
            record: SeqRecord::new("r", s),
            quality: vec![40u8; seq.len()],
        }
    }

    fn aln(pos: i64, mapped: bool) -> SamRecord {
        let mut r = SamRecord::unmapped("r");
        if mapped {
            r.flags = SamFlags(0);
            r.rname = "chr1".to_string();
            r.pos = pos;
            r.mapq = 60;
            r.cigar = Cigar::parse("10M").unwrap();
            r.seq = "ACGTACGTAC".to_string();
            r.qual = "IIIIIIIIII".to_string();
        }
        r
    }

    #[test]
    fn read_qc_summary() {
        let reads = vec![fq("ACGTACGT"), fq("GGCCGGCC")];
        let report = FastqcReport::analyze(&reads);
        let summary = ReadQcSummary::from_report(&report);
        assert_eq!(summary.n_reads, 2);
        assert_eq!(summary.total_bases, 16);
    }

    #[test]
    fn mapping_summary_counts() {
        let recs = vec![aln(1, true), aln(20, true), aln(0, false)];
        let m = MappingSummary::from_records(&recs);
        assert_eq!(m.total_records, 3);
        assert_eq!(m.mapped, 2);
        assert_eq!(m.unmapped, 1);
        assert!((m.mapping_rate - 2.0 / 3.0).abs() < 1e-9);
        assert_eq!(m.mean_mapq, 60.0);
    }

    #[test]
    fn empty_report_has_no_sections() {
        let report = build_report(&[], &[], None, &ContigLengths::new());
        assert!(report.read_qc.is_none());
        assert!(report.mapping.is_none());
        assert!(report.variants.is_none());
        assert_eq!(report.mean_depth, 0.0);
    }

    #[test]
    fn full_report_populated() {
        let reads = vec![fq("ACGTACGT"), fq("GGCCGGCC")];
        let aligns = vec![aln(1, true), aln(5, true)];
        let report = build_report(&reads, &aligns, None, &ContigLengths::new());
        assert!(report.read_qc.is_some());
        assert!(report.mapping.is_some());
        assert!(report.mean_depth > 0.0);
        assert!(report.breadth > 0.0);
        // The summary text mentions every section.
        let text = report.summary_text();
        assert!(text.contains("reads:"));
        assert!(text.contains("mapping:"));
        assert!(text.contains("coverage:"));
    }

    #[test]
    fn call_and_report_end_to_end() {
        // Build a SAM file: 20 reads at the same locus, half showing a
        // variant base at one position.
        let mut sam = SamFile::new();
        for i in 0..20 {
            let mut r = SamRecord::unmapped(format!("r{i}"));
            r.flags = SamFlags(0);
            r.rname = "chr1".to_string();
            r.pos = 1;
            r.mapq = 60;
            r.cigar = Cigar::parse("8M").unwrap();
            // Position 4 (0-based 3) is A for half, G for half.
            r.seq = if i % 2 == 0 {
                "ACGAACGT".to_string()
            } else {
                "ACGGACGT".to_string()
            };
            r.qual = "IIIIIIII".to_string();
            sam.records.push(r);
        }
        let mut refr = Reference::new();
        refr.add("chr1", "ACGAACGT");
        let (variants, report) = call_and_report(&sam, &refr, &CallParams::default()).unwrap();
        // One het SNV expected at position 4.
        assert_eq!(variants.len(), 1);
        assert!(report.variants.is_some());
        assert_eq!(report.variants.unwrap().total, 1);
    }
}
