//! Active-region detection — find windows worth reassembling.
//!
//! The first stage of the GATK HaplotypeCaller is *active-region
//! detection*: the genome is scanned for windows that show **evidence
//! of variation** (mismatches, indels, soft clips, low base quality)
//! above a quality-weighted threshold. Quiet regions skip the
//! expensive reassembly stage and contribute zero calls.
//!
//! This module implements that scan on top of the existing
//! [`crate::format::pileup::PileupColumn`] stream. The
//! per-column **activity score** is a Phred-weighted sum of:
//!
//! - mismatch evidence — `Σ (1 − e_i)` over reads whose base differs
//!   from the reference base, with `e_i` the per-base error;
//! - indel evidence — `count(deletion_placeholder) + count(insertion)`
//!   each scaled by a small constant (indels are rarer than mismatches
//!   on the per-base scale so they get more weight per event);
//! - soft-clip / quality evidence — currently folded into the per-base
//!   mismatch term (a soft-clipped base never reaches the pileup, by
//!   construction of the pileup engine).
//!
//! Columns whose normalised activity exceeds a threshold are flagged
//! *active*; a contiguous run of active columns (with a small allowable
//! gap) plus a configurable left/right margin becomes one
//! [`ActiveRegion`].

use crate::format::pileup::PileupColumn;

/// One genomic window with variation evidence to reassemble.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ActiveRegion {
    /// The reference contig.
    pub chrom: String,
    /// Inclusive 1-based start position.
    pub start: i64,
    /// Inclusive 1-based end position.
    pub end: i64,
}

impl ActiveRegion {
    /// Inclusive span in base pairs.
    pub fn len(&self) -> usize {
        (self.end - self.start + 1).max(0) as usize
    }

    /// `true` when the region carries no positions.
    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }

    /// `true` when `pos` lies in `[start, end]` on `chrom`.
    pub fn contains(&self, chrom: &str, pos: i64) -> bool {
        chrom == self.chrom && self.start <= pos && pos <= self.end
    }
}

/// Tunable thresholds for [`detect_active_regions`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ActiveRegionParams {
    /// Minimum per-column activity score to flag the column active.
    /// A score around `2.0` corresponds to two high-quality mismatching
    /// reads.
    pub activity_threshold: f64,
    /// Number of contiguous non-active columns tolerated inside one
    /// region — small gaps inside otherwise-active windows would
    /// otherwise fragment the region.
    pub max_inner_gap: usize,
    /// Bases added to either side of the active span to give the local
    /// reassembler some flanking reference context to anchor on.
    pub flank: usize,
    /// Maximum region length. Larger candidate regions are split into
    /// chunks — bounded so the local assembler stays cheap.
    pub max_region_len: usize,
    /// Per-mismatch weight (each effective mismatch contributes
    /// `(1 − e_i) · mismatch_weight` to the column score).
    pub mismatch_weight: f64,
    /// Per-indel weight applied to deletion placeholders and to
    /// insertions attached to a base.
    pub indel_weight: f64,
}

impl Default for ActiveRegionParams {
    /// Reasonable defaults for ~30× germline data.
    fn default() -> Self {
        ActiveRegionParams {
            activity_threshold: 1.5,
            max_inner_gap: 5,
            flank: 25,
            max_region_len: 300,
            mismatch_weight: 1.0,
            indel_weight: 3.0,
        }
    }
}

/// Computes one column's activity score.
fn column_activity(col: &PileupColumn, params: &ActiveRegionParams) -> f64 {
    let ref_base = col.ref_base.to_ascii_uppercase();
    let mut score = 0.0f64;
    for b in &col.bases {
        // Deletion placeholder counts as indel evidence.
        if b.is_deletion() {
            score += params.indel_weight;
            continue;
        }
        // Insertion attached to this base counts as indel evidence.
        if !b.insertion.is_empty() {
            score += params.indel_weight;
        }
        let bb = b.base.to_ascii_uppercase();
        if matches!(ref_base, b'A' | b'C' | b'G' | b'T') && bb != ref_base {
            // Mismatch — weight by base quality.
            let q = b.quality.min(60) as f64;
            let e = 10f64.powf(-q / 10.0).clamp(1e-6, 0.75);
            score += params.mismatch_weight * (1.0 - e);
        }
    }
    score
}

/// Scans a sorted column stream and returns the active regions.
///
/// Columns must be sorted by `(chrom, pos)` — the output of
/// [`crate::format::pileup::build_pileup`] satisfies this. Columns from
/// different contigs are processed independently. Returned regions are
/// sorted by `(chrom, start)` and non-overlapping.
pub fn detect_active_regions(
    columns: &[PileupColumn],
    params: &ActiveRegionParams,
) -> Vec<ActiveRegion> {
    let mut regions: Vec<ActiveRegion> = Vec::new();
    // Split into contiguous per-chrom slices to keep gap-merging within
    // a single contig.
    let mut start = 0usize;
    while start < columns.len() {
        let chrom = &columns[start].chrom;
        let mut end = start + 1;
        while end < columns.len() && &columns[end].chrom == chrom {
            end += 1;
        }
        let slice = &columns[start..end];
        let chrom_regions = detect_one_contig(slice, params);
        regions.extend(chrom_regions);
        start = end;
    }
    regions
}

/// Per-contig active-region detection — runs on a slice already
/// restricted to one chromosome and sorted by `pos`.
fn detect_one_contig(columns: &[PileupColumn], params: &ActiveRegionParams) -> Vec<ActiveRegion> {
    if columns.is_empty() {
        return Vec::new();
    }
    let chrom = &columns[0].chrom;

    // Active markers per column index.
    let active: Vec<bool> = columns
        .iter()
        .map(|c| column_activity(c, params) >= params.activity_threshold)
        .collect();

    // Group active columns into runs, allowing up to max_inner_gap
    // non-active columns between two actives.
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut i = 0usize;
    while i < active.len() {
        if !active[i] {
            i += 1;
            continue;
        }
        let mut j = i;
        // Extend the run.
        while j + 1 < active.len() {
            if active[j + 1] {
                j += 1;
                continue;
            }
            // Look ahead for the next active within max_inner_gap.
            let mut k = j + 1;
            let mut gap = 0usize;
            while k < active.len() && gap < params.max_inner_gap && !active[k] {
                k += 1;
                gap += 1;
            }
            if k < active.len() && active[k] {
                j = k;
            } else {
                break;
            }
        }
        runs.push((i, j));
        i = j + 1;
    }

    // Translate each run to a 1-based [start, end] region, adding the
    // flank. Then merge runs whose flanked spans overlap; only after
    // that do we split overly-long merged regions into chunks (so the
    // split is the final step and chunks do not get re-merged).
    let mut regions: Vec<ActiveRegion> = Vec::new();
    for (s, e) in runs {
        let span_start = columns[s].pos - params.flank as i64;
        let span_end = columns[e].pos + params.flank as i64;
        let start = span_start.max(1);
        let end = span_end.max(start);
        regions.push(ActiveRegion {
            chrom: chrom.clone(),
            start,
            end,
        });
    }

    // Merge overlapping / touching flanked spans.
    regions.sort_by_key(|r| r.start);
    let mut merged: Vec<ActiveRegion> = Vec::new();
    for r in regions {
        if let Some(last) = merged.last_mut() {
            if r.start <= last.end + 1 && r.chrom == last.chrom {
                last.end = last.end.max(r.end);
                continue;
            }
        }
        merged.push(r);
    }

    // Split any region longer than max_region_len into chunks. This is
    // the final pass — no further merging happens after it.
    let mut out: Vec<ActiveRegion> = Vec::new();
    for r in merged {
        let mut cur = r.start;
        while cur <= r.end {
            let chunk_end = (cur + params.max_region_len as i64 - 1).min(r.end);
            out.push(ActiveRegion {
                chrom: r.chrom.clone(),
                start: cur,
                end: chunk_end,
            });
            cur = chunk_end + 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::pileup::{PileupBase, PileupColumn};

    fn pbase(base: u8, q: u8, rev: bool, ins: &[u8]) -> PileupBase {
        PileupBase {
            base,
            quality: q,
            reverse: rev,
            mapq: 60,
            read_pos: 0,
            insertion: ins.to_vec(),
        }
    }

    fn col(chrom: &str, pos: i64, refb: u8, bases: Vec<PileupBase>) -> PileupColumn {
        PileupColumn {
            chrom: chrom.to_string(),
            pos,
            ref_base: refb,
            bases,
        }
    }

    fn calm_column(pos: i64) -> PileupColumn {
        let bases: Vec<PileupBase> = (0..15).map(|_| pbase(b'A', 35, false, &[])).collect();
        col("chr1", pos, b'A', bases)
    }

    fn snv_column(pos: i64) -> PileupColumn {
        // 15 reads — 8 ref A, 7 alt G at Phred 35.
        let mut bases: Vec<PileupBase> = (0..8).map(|_| pbase(b'A', 35, false, &[])).collect();
        for _ in 0..7 {
            bases.push(pbase(b'G', 35, false, &[]));
        }
        col("chr1", pos, b'A', bases)
    }

    #[test]
    fn calm_region_yields_no_active_regions() {
        let cols: Vec<_> = (1..=20).map(calm_column).collect();
        let regs = detect_active_regions(&cols, &ActiveRegionParams::default());
        assert!(
            regs.is_empty(),
            "calm slice produced {} regions",
            regs.len()
        );
    }

    #[test]
    fn single_snv_column_makes_one_region() {
        let mut cols: Vec<_> = (1..=20).map(calm_column).collect();
        cols[10] = snv_column(11);
        let regs = detect_active_regions(&cols, &ActiveRegionParams::default());
        assert_eq!(regs.len(), 1, "expected one active region, got {regs:?}");
        let r = &regs[0];
        assert!(r.contains("chr1", 11));
        // Flank applied: start <= 11 - default_flank? At least the region
        // extends down to position 1 with the default 25 bp flank.
        assert!(r.start < 11);
        assert!(r.end > 11);
    }

    #[test]
    fn two_close_active_columns_become_one_region() {
        let mut cols: Vec<_> = (1..=40).map(calm_column).collect();
        cols[5] = snv_column(6);
        cols[8] = snv_column(9);
        let regs = detect_active_regions(&cols, &ActiveRegionParams::default());
        assert_eq!(regs.len(), 1);
    }

    #[test]
    fn two_far_active_columns_become_two_regions_if_far_enough() {
        // Use a small flank so distant active columns don't merge.
        let params = ActiveRegionParams {
            flank: 1,
            max_inner_gap: 1,
            ..ActiveRegionParams::default()
        };
        let mut cols: Vec<_> = (1..=200).map(calm_column).collect();
        cols[5] = snv_column(6);
        cols[100] = snv_column(101);
        let regs = detect_active_regions(&cols, &params);
        assert_eq!(regs.len(), 2);
    }

    #[test]
    fn indel_column_activates() {
        // 15 reads with deletion placeholder (`*`).
        let bases: Vec<_> = (0..15).map(|_| pbase(b'*', 0, false, &[])).collect();
        let mut cols: Vec<_> = (1..=20).map(calm_column).collect();
        cols[7] = col("chr1", 8, b'A', bases);
        let regs = detect_active_regions(&cols, &ActiveRegionParams::default());
        assert_eq!(regs.len(), 1);
        assert!(regs[0].contains("chr1", 8));
    }

    #[test]
    fn long_active_run_splits_at_max_region_len() {
        // 200 active columns in a row -> two regions when max=80.
        let params = ActiveRegionParams {
            max_region_len: 80,
            flank: 0,
            max_inner_gap: 1,
            ..ActiveRegionParams::default()
        };
        let cols: Vec<_> = (1..=200).map(snv_column).collect();
        let regs = detect_active_regions(&cols, &params);
        assert!(regs.len() >= 3, "expected splitting, got {regs:?}");
    }

    #[test]
    fn per_chrom_grouping() {
        let mut cols = Vec::new();
        for p in 1..=15 {
            cols.push(calm_column(p));
        }
        // chr1 SNV
        cols[5] = snv_column(6);
        // chr2 calm + SNV at pos 6
        for p in 1..=15 {
            let mut c = calm_column(p);
            c.chrom = "chr2".into();
            cols.push(c);
        }
        let mut sn = snv_column(6);
        sn.chrom = "chr2".into();
        cols[15 + 5] = sn;
        let regs = detect_active_regions(&cols, &ActiveRegionParams::default());
        assert_eq!(regs.len(), 2);
        assert_eq!(regs[0].chrom, "chr1");
        assert_eq!(regs[1].chrom, "chr2");
    }

    #[test]
    fn empty_columns_yield_empty_regions() {
        let regs = detect_active_regions(&[], &ActiveRegionParams::default());
        assert!(regs.is_empty());
    }

    #[test]
    fn region_contains_pos() {
        let r = ActiveRegion {
            chrom: "chr1".into(),
            start: 10,
            end: 20,
        };
        assert!(r.contains("chr1", 10));
        assert!(r.contains("chr1", 20));
        assert!(r.contains("chr1", 15));
        assert!(!r.contains("chr1", 9));
        assert!(!r.contains("chr2", 15));
    }
}
