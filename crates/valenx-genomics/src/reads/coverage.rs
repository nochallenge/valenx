//! Coverage and depth computation across a reference.
//!
//! Depth — the number of reads spanning each reference base — is the
//! single most-consulted QC metric for any resequencing experiment
//! (samtools `depth`, `bedtools genomecov`, `mosdepth`). This module
//! walks a set of [`crate::format::sam::SamRecord`]
//! alignments and produces a [`DepthProfile`]: the per-base depth
//! array of each contig plus summary statistics (mean, median,
//! breadth, the fraction of bases at or above a depth threshold).

use crate::format::sam::{CigarKind, SamRecord};
use std::collections::BTreeMap;

/// The per-base depth array of one contig.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContigDepth {
    /// The contig name.
    pub name: String,
    /// `depth[i]` is the read depth at the 0-based position `i`. The
    /// vector length is the highest covered position + 1 (or the
    /// declared contig length when one was supplied).
    pub depth: Vec<u32>,
}

impl ContigDepth {
    /// Mean depth across every position in the array.
    pub fn mean_depth(&self) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        self.depth.iter().map(|&d| d as f64).sum::<f64>() / self.depth.len() as f64
    }

    /// Median depth across every position.
    pub fn median_depth(&self) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        let mut v = self.depth.clone();
        v.sort_unstable();
        let mid = v.len() / 2;
        if v.len() % 2 == 1 {
            v[mid] as f64
        } else {
            (v[mid - 1] as f64 + v[mid] as f64) / 2.0
        }
    }

    /// Maximum depth seen anywhere on the contig.
    pub fn max_depth(&self) -> u32 {
        self.depth.iter().copied().max().unwrap_or(0)
    }

    /// **Breadth of coverage** — the fraction of positions with depth
    /// `>= 1`.
    pub fn breadth(&self) -> f64 {
        self.fraction_at_least(1)
    }

    /// Fraction of positions whose depth is at least `min_depth`.
    pub fn fraction_at_least(&self, min_depth: u32) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        let n = self.depth.iter().filter(|&&d| d >= min_depth).count();
        n as f64 / self.depth.len() as f64
    }

    /// Length of the contig as covered by the profile.
    pub fn len(&self) -> usize {
        self.depth.len()
    }

    /// `true` when the depth array is empty.
    pub fn is_empty(&self) -> bool {
        self.depth.is_empty()
    }
}

/// A genome-wide depth profile — one [`ContigDepth`] per contig.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DepthProfile {
    /// Per-contig depth arrays, keyed by contig name.
    pub contigs: BTreeMap<String, ContigDepth>,
}

impl DepthProfile {
    /// Looks up a contig's depth array.
    pub fn contig(&self, name: &str) -> Option<&ContigDepth> {
        self.contigs.get(name)
    }

    /// Mean depth averaged over every base of every contig.
    pub fn genome_mean_depth(&self) -> f64 {
        let mut total = 0.0f64;
        let mut bases = 0usize;
        for c in self.contigs.values() {
            total += c.depth.iter().map(|&d| d as f64).sum::<f64>();
            bases += c.depth.len();
        }
        if bases == 0 {
            0.0
        } else {
            total / bases as f64
        }
    }

    /// Genome-wide breadth — fraction of all bases with depth `>= 1`.
    pub fn genome_breadth(&self) -> f64 {
        self.genome_fraction_at_least(1)
    }

    /// Genome-wide fraction of bases with depth `>= min_depth`.
    pub fn genome_fraction_at_least(&self, min_depth: u32) -> f64 {
        let mut hit = 0usize;
        let mut bases = 0usize;
        for c in self.contigs.values() {
            hit += c.depth.iter().filter(|&&d| d >= min_depth).count();
            bases += c.depth.len();
        }
        if bases == 0 {
            0.0
        } else {
            hit as f64 / bases as f64
        }
    }

    /// Renders a samtools-`depth`-style text dump: `contig  pos
    /// depth`, 1-based positions, every position emitted.
    pub fn to_depth_text(&self) -> String {
        let mut s = String::new();
        for c in self.contigs.values() {
            for (i, &d) in c.depth.iter().enumerate() {
                s.push_str(&format!("{}\t{}\t{}\n", c.name, i + 1, d));
            }
        }
        s
    }
}

/// Optional declared contig lengths (so a depth array extends to the
/// full contig even past the last read).
pub type ContigLengths = BTreeMap<String, usize>;

/// Computes the [`DepthProfile`] of a set of alignments.
///
/// Only `M` / `=` / `X` / `D` CIGAR ops contribute to depth (a
/// reference base is "covered" by a read aligned across it, including
/// a deletion within the read — matching `samtools depth -a` default
/// behaviour for `D`); `N` (intron skips), `I`, `S` and `H` do not.
/// `min_mapq` excludes low-confidence reads; unmapped reads are
/// skipped. When `lengths` declares a contig length, the depth array
/// is padded to it.
pub fn compute_depth(records: &[SamRecord], min_mapq: u8, lengths: &ContigLengths) -> DepthProfile {
    // contig -> Vec<u32> grown lazily
    let mut depth: BTreeMap<String, Vec<u32>> = BTreeMap::new();

    // Seed declared lengths so empty contigs still appear.
    for (name, &len) in lengths {
        depth.entry(name.clone()).or_insert_with(|| vec![0; len]);
    }

    for rec in records {
        if rec.is_unmapped() || rec.pos <= 0 || rec.cigar.is_empty() {
            continue;
        }
        if rec.mapq < min_mapq {
            continue;
        }
        let arr = depth.entry(rec.rname.clone()).or_default();
        let mut ref_pos = (rec.pos - 1) as usize; // 0-based
        for op in &rec.cigar.ops {
            let n = op.len as usize;
            match op.kind {
                CigarKind::Match | CigarKind::Equal | CigarKind::Diff | CigarKind::Del => {
                    let end = ref_pos + n;
                    if arr.len() < end {
                        arr.resize(end, 0);
                    }
                    for d in arr.iter_mut().take(end).skip(ref_pos) {
                        *d = d.saturating_add(1);
                    }
                    ref_pos = end;
                }
                CigarKind::Skip => {
                    // Intron — advance the reference, add no depth.
                    ref_pos += n;
                }
                CigarKind::Ins | CigarKind::SoftClip | CigarKind::HardClip | CigarKind::Pad => {}
            }
        }
    }

    let contigs = depth
        .into_iter()
        .map(|(name, depth)| (name.clone(), ContigDepth { name, depth }))
        .collect();
    DepthProfile { contigs }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::sam::{Cigar, SamFlags};

    fn aln(pos: i64, cigar: &str, contig: &str) -> SamRecord {
        let mut r = SamRecord::unmapped("r");
        r.flags = SamFlags(0);
        r.rname = contig.to_string();
        r.pos = pos;
        r.mapq = 60;
        r.cigar = Cigar::parse(cigar).unwrap();
        // SEQ length must match the CIGAR query length.
        let qlen = r.cigar.query_len();
        r.seq = "A".repeat(qlen);
        r.qual = "I".repeat(qlen);
        r
    }

    #[test]
    fn single_read_depth() {
        let recs = vec![aln(1, "10M", "chr1")];
        let p = compute_depth(&recs, 0, &ContigLengths::new());
        let c = p.contig("chr1").unwrap();
        assert_eq!(c.depth.len(), 10);
        assert!(c.depth.iter().all(|&d| d == 1));
    }

    #[test]
    fn overlapping_reads_stack() {
        let recs = vec![aln(1, "10M", "chr1"), aln(5, "10M", "chr1")];
        let p = compute_depth(&recs, 0, &ContigLengths::new());
        let c = p.contig("chr1").unwrap();
        // Positions 5..10 (0-based 4..10) are covered twice.
        assert_eq!(c.depth[0], 1);
        assert_eq!(c.depth[4], 2);
        assert_eq!(c.depth[9], 2);
        assert_eq!(c.depth[10], 1);
        assert_eq!(c.depth.len(), 14);
    }

    #[test]
    fn intron_skip_leaves_a_gap() {
        // 5M100N5M: positions 6..105 get no depth.
        let recs = vec![aln(1, "5M100N5M", "chr1")];
        let p = compute_depth(&recs, 0, &ContigLengths::new());
        let c = p.contig("chr1").unwrap();
        assert_eq!(c.depth[0], 1);
        assert_eq!(c.depth[4], 1);
        assert_eq!(c.depth[5], 0); // inside the intron
        assert_eq!(c.depth[104], 0);
        assert_eq!(c.depth[105], 1);
    }

    #[test]
    fn deletion_counts_as_covered() {
        // 5M3D5M: the 3 deleted bases are still spanned by the read.
        let recs = vec![aln(1, "5M3D5M", "chr1")];
        let p = compute_depth(&recs, 0, &ContigLengths::new());
        let c = p.contig("chr1").unwrap();
        assert_eq!(c.depth.len(), 13);
        assert!(c.depth.iter().all(|&d| d == 1));
    }

    #[test]
    fn mapq_filter_excludes() {
        let mut low = aln(1, "10M", "chr1");
        low.mapq = 5;
        let p = compute_depth(&[low], 20, &ContigLengths::new());
        // Contig has no covered bases.
        assert!(p.contig("chr1").map(|c| c.is_empty()).unwrap_or(true));
    }

    #[test]
    fn declared_length_pads_array() {
        let mut lens = ContigLengths::new();
        lens.insert("chr1".to_string(), 50);
        let recs = vec![aln(1, "10M", "chr1")];
        let p = compute_depth(&recs, 0, &lens);
        let c = p.contig("chr1").unwrap();
        assert_eq!(c.depth.len(), 50);
        assert_eq!(c.breadth(), 10.0 / 50.0);
    }

    #[test]
    fn statistics() {
        let recs = vec![aln(1, "10M", "chr1"), aln(1, "10M", "chr1")];
        let p = compute_depth(&recs, 0, &ContigLengths::new());
        let c = p.contig("chr1").unwrap();
        assert_eq!(c.mean_depth(), 2.0);
        assert_eq!(c.median_depth(), 2.0);
        assert_eq!(c.max_depth(), 2);
        assert_eq!(c.fraction_at_least(2), 1.0);
        assert_eq!(c.fraction_at_least(3), 0.0);
        assert_eq!(p.genome_mean_depth(), 2.0);
        assert_eq!(p.genome_breadth(), 1.0);
    }
}
