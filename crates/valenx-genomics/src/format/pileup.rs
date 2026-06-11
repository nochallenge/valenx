//! Pileup generation from aligned reads.
//!
//! A *pileup* turns a set of [`crate::format::sam::SamRecord`]
//! alignments into a per-reference-position summary of the read bases
//! covering that position. It is the data structure every classical
//! pileup-based variant caller (samtools / bcftools `mpileup`,
//! VarScan) consumes, and the input to this crate's own caller in
//! [`crate::variant`].
//!
//! This module produces a structured [`PileupColumn`] per position
//! (carrying typed [`PileupBase`] evidence — base, quality, strand,
//! indels) and can render the classic samtools text-pileup format.
//!
//! ## v1 scope
//!
//! The CIGAR walk handles `M` / `=` / `X` / `I` / `D` / `S` / `N` /
//! `H`. `P` (padding) is skipped. Per-base qualities are read from the
//! SAM `QUAL` string (Phred+33); a read with `QUAL == "*"` contributes
//! bases with a default quality. Reference skips (`N`, used by spliced
//! RNA-seq aligners) leave a gap, not a deletion.

use crate::error::{GenomicsError, Result};
use crate::format::sam::{CigarKind, SamRecord};
use std::collections::BTreeMap;

/// Default Phred quality assigned to a base when the read's `QUAL`
/// field is the `*` placeholder.
pub const DEFAULT_BASE_QUALITY: u8 = 30;

/// Defensive per-record cap on the total reference span the builder
/// will materialise into pileup columns from a single SAM record.
/// 100 Mb — comfortably past any real chromosome arm — bounds the
/// `BTreeMap` insert count even when an upstream CIGAR cap is
/// bypassed (a future caller that hand-builds `SamRecord` values
/// without going through `Cigar::parse` still gets this guard).
pub const MAX_PILEUP_SPAN: usize = 100_000_000;

/// One read base observed at a pileup position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PileupBase {
    /// The read base (uppercased ASCII; `*` marks a deletion placeholder).
    pub base: u8,
    /// Phred quality of the base.
    pub quality: u8,
    /// `true` when the read maps to the reverse strand.
    pub reverse: bool,
    /// Mapping quality of the read carrying this base.
    pub mapq: u8,
    /// The 0-based read offset this base came from (useful for
    /// read-position bias tests).
    pub read_pos: usize,
    /// An insertion that begins *immediately after* this position, as
    /// the inserted bases (empty when none).
    pub insertion: Vec<u8>,
}

impl PileupBase {
    /// `true` when this entry is a deletion placeholder (`*`).
    pub fn is_deletion(&self) -> bool {
        self.base == b'*'
    }
}

/// One pileup column — every read base covering a single 1-based
/// reference position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PileupColumn {
    /// The reference contig name.
    pub chrom: String,
    /// The 1-based reference position.
    pub pos: i64,
    /// The reference base, if a reference was supplied (`N` otherwise).
    pub ref_base: u8,
    /// Every read base at this position.
    pub bases: Vec<PileupBase>,
}

impl PileupColumn {
    /// Read depth at this column (every covering read, deletions
    /// included).
    pub fn depth(&self) -> usize {
        self.bases.len()
    }

    /// Depth counting only non-deletion bases.
    pub fn base_depth(&self) -> usize {
        self.bases.iter().filter(|b| !b.is_deletion()).count()
    }

    /// Per-allele observation counts (A / C / G / T / N / `*`), keyed by
    /// the uppercased base byte.
    pub fn allele_counts(&self) -> BTreeMap<u8, usize> {
        let mut m = BTreeMap::new();
        for b in &self.bases {
            *m.entry(b.base).or_insert(0) += 1;
        }
        m
    }

    /// Count of forward-strand observations of a given base.
    pub fn strand_counts(&self, base: u8) -> (usize, usize) {
        let mut fwd = 0;
        let mut rev = 0;
        for b in &self.bases {
            if b.base == base {
                if b.reverse {
                    rev += 1;
                } else {
                    fwd += 1;
                }
            }
        }
        (fwd, rev)
    }

    /// Mean base quality across all non-deletion observations.
    pub fn mean_quality(&self) -> f64 {
        let quals: Vec<u8> = self
            .bases
            .iter()
            .filter(|b| !b.is_deletion())
            .map(|b| b.quality)
            .collect();
        if quals.is_empty() {
            return 0.0;
        }
        quals.iter().map(|&q| q as f64).sum::<f64>() / quals.len() as f64
    }

    /// Renders this column as one samtools-style pileup text line:
    /// `chrom  pos  refbase  depth  read-bases  base-qualities`.
    pub fn to_pileup_line(&self) -> String {
        let mut read_str = String::new();
        let mut qual_str = String::new();
        for b in &self.bases {
            if b.is_deletion() {
                read_str.push('*');
            } else if b.base == self.ref_base.to_ascii_uppercase() {
                // Match: `.` forward, `,` reverse.
                read_str.push(if b.reverse { ',' } else { '.' });
            } else {
                // Mismatch: uppercase forward, lowercase reverse.
                let c = b.base as char;
                read_str.push(if b.reverse {
                    c.to_ascii_lowercase()
                } else {
                    c.to_ascii_uppercase()
                });
            }
            // Insertion suffix `+<len><bases>`.
            if !b.insertion.is_empty() {
                read_str.push('+');
                read_str.push_str(&b.insertion.len().to_string());
                for &ib in &b.insertion {
                    read_str.push(ib as char);
                }
            }
            // Phred+33 quality.
            qual_str.push((b.quality.min(93) + 33) as char);
        }
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.pos,
            self.ref_base as char,
            self.depth(),
            read_str,
            qual_str,
        )
    }
}

/// A reference sequence keyed by contig name — what the pileup engine
/// needs to fill in `ref_base`.
#[derive(Clone, Debug, Default)]
pub struct Reference {
    contigs: BTreeMap<String, Vec<u8>>,
}

impl Reference {
    /// An empty reference set.
    pub fn new() -> Self {
        Reference::default()
    }

    /// Inserts a contig (sequence uppercased on insert).
    pub fn add(&mut self, name: impl Into<String>, seq: impl AsRef<[u8]>) {
        self.contigs.insert(
            name.into(),
            seq.as_ref().iter().map(|b| b.to_ascii_uppercase()).collect(),
        );
    }

    /// The 0-based byte of a contig, or `b'N'` when out of range or the
    /// contig is unknown.
    pub fn base_at(&self, contig: &str, idx0: usize) -> u8 {
        self.contigs
            .get(contig)
            .and_then(|c| c.get(idx0).copied())
            .unwrap_or(b'N')
    }

    /// Length of a contig (`0` when unknown).
    pub fn contig_len(&self, contig: &str) -> usize {
        self.contigs.get(contig).map(|c| c.len()).unwrap_or(0)
    }

    /// `true` when the contig is present.
    pub fn has(&self, contig: &str) -> bool {
        self.contigs.contains_key(contig)
    }
}

/// Builds the full set of [`PileupColumn`]s spanned by `records`.
///
/// Records must be mapped (an unmapped record is skipped). Columns are
/// returned sorted by `(chrom, pos)`. `min_mapq` filters out reads
/// whose mapping quality is below the threshold. When `reference` has a
/// matching contig, `ref_base` is filled; otherwise it is `N`.
pub fn build_pileup(
    records: &[SamRecord],
    reference: &Reference,
    min_mapq: u8,
) -> Result<Vec<PileupColumn>> {
    // (chrom, pos) -> bases
    let mut columns: BTreeMap<(String, i64), Vec<PileupBase>> = BTreeMap::new();

    for rec in records {
        if rec.is_unmapped() || rec.pos <= 0 || rec.cigar.is_empty() {
            continue;
        }
        if rec.mapq < min_mapq {
            continue;
        }
        if rec.seq.is_empty() {
            continue;
        }
        // Round-6 defensive guard: even with `Cigar::parse`'s
        // `MAX_CIGAR_OP_LEN` cap in place, a record assembled by a
        // hand-rolled caller could chain enough mid-cap ops to spill
        // past 100 M reference bases. Reject the whole record before
        // we walk the CIGAR and start filling `columns`.
        let total_ref_span = rec.cigar.ref_len();
        if total_ref_span > MAX_PILEUP_SPAN {
            return Err(GenomicsError::invalid_record(
                "sam",
                format!(
                    "CIGAR reference span {total_ref_span} exceeds the \
                     {MAX_PILEUP_SPAN}-base pileup cap (amplification-DoS guard)"
                ),
            ));
        }
        let seq = rec.seq.as_bytes();
        let qual: Vec<u8> = if rec.qual.is_empty() {
            vec![DEFAULT_BASE_QUALITY; seq.len()]
        } else {
            rec.qual
                .as_bytes()
                .iter()
                .map(|&q| q.saturating_sub(33))
                .collect()
        };
        if qual.len() != seq.len() {
            return Err(GenomicsError::invalid_record(
                "sam",
                "QUAL length disagrees with SEQ during pileup",
            ));
        }
        // The CIGAR's query-consuming length must equal SEQ, or the Match arm
        // below would index `seq` / `qual` out of bounds. `SamFile::parse`
        // enforces this through `SamRecord::validate`, but the public
        // `build_pileup` / `call_haplotype_variants` entry points also accept
        // caller-built records that never went through it.
        let qlen = rec.cigar.query_len();
        if qlen != seq.len() {
            return Err(GenomicsError::invalid_record(
                "sam",
                format!("CIGAR query length {qlen} != SEQ length {}", seq.len()),
            ));
        }

        let mut ref_pos = rec.pos; // 1-based current reference position
        let mut read_pos = 0usize; // 0-based current read offset

        let ops = &rec.cigar.ops;
        for (oi, op) in ops.iter().enumerate() {
            let n = op.len as usize;
            match op.kind {
                CigarKind::Match | CigarKind::Equal | CigarKind::Diff => {
                    for k in 0..n {
                        let rp = read_pos + k;
                        let base = seq[rp].to_ascii_uppercase();
                        // An insertion that begins right after this
                        // aligned base is attached to it.
                        let insertion = if k + 1 == n {
                            next_insertion(ops, oi, seq, read_pos + n)
                        } else {
                            Vec::new()
                        };
                        columns
                            .entry((rec.rname.clone(), ref_pos + k as i64))
                            .or_default()
                            .push(PileupBase {
                                base,
                                quality: qual[rp],
                                reverse: rec.flags.is_reverse(),
                                mapq: rec.mapq,
                                read_pos: rp,
                                insertion,
                            });
                    }
                    ref_pos += n as i64;
                    read_pos += n;
                }
                CigarKind::Del => {
                    // Each deleted reference base gets a `*` placeholder.
                    for k in 0..n {
                        columns
                            .entry((rec.rname.clone(), ref_pos + k as i64))
                            .or_default()
                            .push(PileupBase {
                                base: b'*',
                                quality: 0,
                                reverse: rec.flags.is_reverse(),
                                mapq: rec.mapq,
                                read_pos,
                                insertion: Vec::new(),
                            });
                    }
                    ref_pos += n as i64;
                }
                CigarKind::Ins => {
                    // Consumed by `next_insertion` on the preceding M.
                    read_pos += n;
                }
                CigarKind::SoftClip => {
                    read_pos += n;
                }
                CigarKind::Skip => {
                    // Intron: advance the reference, leave no column.
                    ref_pos += n as i64;
                }
                CigarKind::HardClip | CigarKind::Pad => {
                    // Consume nothing from SEQ; nothing to record.
                }
            }
        }
    }

    // Flatten into sorted columns, filling reference bases.
    let mut out = Vec::with_capacity(columns.len());
    for ((chrom, pos), bases) in columns {
        let ref_base = if pos > 0 {
            reference.base_at(&chrom, (pos - 1) as usize)
        } else {
            b'N'
        };
        out.push(PileupColumn {
            chrom,
            pos,
            ref_base,
            bases,
        });
    }
    Ok(out)
}

/// If the CIGAR op immediately following index `oi` is an insertion,
/// returns the inserted bases read from `seq` starting at `ins_start`.
fn next_insertion(
    ops: &[crate::format::sam::CigarOp],
    oi: usize,
    seq: &[u8],
    ins_start: usize,
) -> Vec<u8> {
    if let Some(next) = ops.get(oi + 1) {
        if next.kind == CigarKind::Ins {
            let end = (ins_start + next.len as usize).min(seq.len());
            return seq[ins_start.min(seq.len())..end]
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect();
        }
    }
    Vec::new()
}

/// Renders a slice of columns as a complete samtools-style pileup
/// document (one line per column, trailing newline).
pub fn to_pileup_text(columns: &[PileupColumn]) -> String {
    let mut s = String::new();
    for c in columns {
        s.push_str(&c.to_pileup_line());
        s.push('\n');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::sam::{Cigar, SamFlags};

    fn mapped(name: &str, pos: i64, cigar: &str, seq: &str, qual: &str) -> SamRecord {
        let mut r = SamRecord::unmapped(name);
        r.flags = SamFlags(0);
        r.rname = "chr1".to_string();
        r.pos = pos;
        r.mapq = 60;
        r.cigar = Cigar::parse(cigar).unwrap();
        r.seq = seq.to_string();
        r.qual = qual.to_string();
        r
    }

    #[test]
    fn simple_match_pileup() {
        // Two reads, both 4M at pos 1, fully overlapping.
        let recs = vec![
            mapped("r1", 1, "4M", "ACGT", "IIII"),
            mapped("r2", 1, "4M", "ACGT", "IIII"),
        ];
        let mut refr = Reference::new();
        refr.add("chr1", "ACGTACGT");
        let cols = build_pileup(&recs, &refr, 0).unwrap();
        assert_eq!(cols.len(), 4);
        assert_eq!(cols[0].pos, 1);
        assert_eq!(cols[0].depth(), 2);
        assert_eq!(cols[0].ref_base, b'A');
    }

    #[test]
    fn mismatch_recorded() {
        let recs = vec![mapped("r1", 1, "4M", "AGGT", "IIII")];
        let mut refr = Reference::new();
        refr.add("chr1", "ACGT");
        let cols = build_pileup(&recs, &refr, 0).unwrap();
        // Position 2 reference is C, read says G.
        let col2 = &cols[1];
        assert_eq!(col2.ref_base, b'C');
        assert_eq!(col2.bases[0].base, b'G');
        let counts = col2.allele_counts();
        assert_eq!(counts.get(&b'G'), Some(&1));
    }

    #[test]
    fn deletion_makes_placeholders() {
        // 2M2D2M: positions 3,4 are deletions.
        let recs = vec![mapped("r1", 1, "2M2D2M", "ACGT", "IIII")];
        let mut refr = Reference::new();
        refr.add("chr1", "ACGTAC");
        let cols = build_pileup(&recs, &refr, 0).unwrap();
        assert_eq!(cols.len(), 6);
        assert!(cols[2].bases[0].is_deletion());
        assert!(cols[3].bases[0].is_deletion());
        assert!(!cols[4].bases[0].is_deletion());
    }

    #[test]
    fn insertion_attached_to_preceding_base() {
        // 2M2I2M: insertion of "GG" after position 2.
        let recs = vec![mapped("r1", 1, "2M2I2M", "ACGGTA", "IIIIII")];
        let mut refr = Reference::new();
        refr.add("chr1", "ACTA");
        let cols = build_pileup(&recs, &refr, 0).unwrap();
        // Reference span is 4 (2M + 2M), so 4 columns.
        assert_eq!(cols.len(), 4);
        assert_eq!(cols[1].bases[0].insertion, b"GG".to_vec());
    }

    #[test]
    fn mapq_filter_excludes_reads() {
        let mut low = mapped("r1", 1, "4M", "ACGT", "IIII");
        low.mapq = 5;
        let cols = build_pileup(&[low], &Reference::new(), 20).unwrap();
        assert!(cols.is_empty());
    }

    #[test]
    fn pileup_text_format() {
        let recs = vec![
            mapped("r1", 1, "4M", "ACGT", "IIII"),
            mapped("r2", 1, "4M", "AGGT", "IIII"),
        ];
        let mut refr = Reference::new();
        refr.add("chr1", "ACGT");
        let cols = build_pileup(&recs, &refr, 0).unwrap();
        let text = to_pileup_text(&cols);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 4);
        // Column 1: both match A -> ".."
        assert!(lines[0].contains("\t..\t"));
        // Column 2: ref C, r1 has C (match -> .), r2 has G (mismatch -> G)
        assert!(lines[1].contains(".G") || lines[1].contains("G."));
    }

    #[test]
    fn strand_counts_split() {
        let mut fwd = mapped("r1", 1, "1M", "A", "I");
        fwd.flags = SamFlags(0);
        let mut rev = mapped("r2", 1, "1M", "A", "I");
        rev.flags = SamFlags(SamFlags::REVERSE);
        let cols = build_pileup(&[fwd, rev], &Reference::new(), 0).unwrap();
        assert_eq!(cols[0].strand_counts(b'A'), (1, 1));
    }

    #[test]
    fn build_pileup_rejects_record_past_max_span() {
        // Round-6 RED→GREEN: even if `Cigar::parse` is bypassed (a
        // hand-rolled SamRecord assembled directly), the pileup
        // builder must defensively reject any record whose CIGAR
        // total reference span exceeds MAX_PILEUP_SPAN — otherwise a
        // ~Mb-long `N` skip + a follow-up 1 M-op chain could still
        // amplify into hundreds of millions of BTreeMap inserts.
        //
        // Build a CIGAR with two `Match` ops totalling MAX + 1 bases.
        // (We construct CigarOp values directly, not via parse, so
        // the parser cap doesn't gate the test.)
        use crate::format::sam::{CigarKind, CigarOp};
        let oversized = MAX_PILEUP_SPAN + 1;
        let mut rec = mapped("r1", 1, "1M", "A", "I");
        // Two equal-ish ops summing past MAX_PILEUP_SPAN.
        let half = (oversized / 2) as u32;
        let rest = oversized as u32 - half;
        rec.cigar = crate::format::sam::Cigar {
            ops: vec![
                CigarOp {
                    len: half,
                    kind: CigarKind::Match,
                },
                CigarOp {
                    len: rest,
                    kind: CigarKind::Match,
                },
            ],
        };
        // The pileup builder must reject this BEFORE iterating the
        // CIGAR — no `Vec<f64; oversized>` allocations get made.
        let err = build_pileup(&[rec], &Reference::new(), 0).unwrap_err();
        match err {
            GenomicsError::InvalidRecord { kind, reason } => {
                assert_eq!(kind, "sam");
                assert!(reason.contains("pileup cap"), "msg: {reason}");
            }
            other => panic!("expected InvalidRecord, got {other:?}"),
        }
    }

    #[test]
    fn build_pileup_rejects_cigar_longer_than_seq() {
        // A hand-built record whose CIGAR consumes more query bases (10M)
        // than SEQ provides (4) must be rejected, not panic on `seq[4]`.
        // `SamFile::parse` catches this via `validate()`, but `build_pileup`
        // is a public entry point that accepts caller-built records directly.
        let rec = mapped("r1", 1, "10M", "ACGT", "IIII");
        let err = build_pileup(&[rec], &Reference::new(), 0).unwrap_err();
        match err {
            GenomicsError::InvalidRecord { kind, reason } => {
                assert_eq!(kind, "sam");
                assert!(reason.contains("CIGAR query length"), "msg: {reason}");
            }
            other => panic!("expected InvalidRecord, got {other:?}"),
        }
    }
}
