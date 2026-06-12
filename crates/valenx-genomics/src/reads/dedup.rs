//! Duplicate marking — coordinate-based and sequence-based.
//!
//! PCR and optical duplicates inflate apparent depth and bias variant
//! calling. Picard `MarkDuplicates` and samtools `markdup` flag them
//! by alignment coordinate; sequence-identical de-duplication (the
//! `clumpify` / `seqkit rmdup` approach) flags them by read content.
//!
//! This module offers both:
//!
//! - [`mark_duplicates_coordinate`] — groups
//!   [`crate::format::sam::SamRecord`] alignments by `(rname, strand,
//!   unclipped 5′ position)`, keeps the highest-quality representative
//!   of each group and sets the SAM `DUPLICATE` flag on the rest.
//! - [`dedup_sequences`] — flags reads sharing an identical sequence
//!   string.
//!
//! ## v1 scope
//!
//! The coordinate model uses the **5′ unclipped** position (`POS`
//! shifted left by a leading soft-clip), which is what Picard keys on
//! for single-end reads. True paired-end duplicate detection — keying
//! on both mates' coordinates — and optical-duplicate detection from
//! flowcell tile coordinates are out of v1 scope; the module documents
//! the gap. Coordinate marking treats every read as single-end.

use crate::format::sam::{CigarKind, SamFlags, SamRecord};
use std::collections::HashMap;

/// A duplicate-group coordinate key: `(rname, is_reverse, unclipped 5′
/// position)`.
type GroupKey = (String, bool, i64);

/// The representative-quality key: `(summed base quality, MAPQ)`.
type QualityKey = (u64, u8);

/// Summary of a duplicate-marking pass.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct DedupStats {
    /// Total records inspected.
    pub total: usize,
    /// Records that are unique representatives (kept).
    pub unique: usize,
    /// Records flagged as duplicates.
    pub duplicates: usize,
    /// Unmapped records (never marked; passed through).
    pub unmapped: usize,
}

impl DedupStats {
    /// Duplicate fraction over the mapped records (`0.0` when none are
    /// mapped).
    pub fn duplicate_rate(&self) -> f64 {
        let mapped = self.unique + self.duplicates;
        if mapped == 0 {
            0.0
        } else {
            self.duplicates as f64 / mapped as f64
        }
    }
}

/// The 5′-most unclipped reference position of a read: `POS` minus any
/// leading soft / hard clip. Reverse-strand reads are keyed at their 3′
/// end shifted by trailing clips so a forward and reverse read of the
/// same fragment land on distinct keys.
fn unclipped_5prime(rec: &SamRecord) -> i64 {
    let reverse = rec.flags.is_reverse();
    if !reverse {
        let leading: i64 = rec
            .cigar
            .ops
            .iter()
            .take_while(|o| matches!(o.kind, CigarKind::SoftClip | CigarKind::HardClip))
            .map(|o| o.len as i64)
            .sum();
        rec.pos - leading
    } else {
        // Reverse strand: anchor at the unclipped 3′ end.
        let trailing: i64 = rec
            .cigar
            .ops
            .iter()
            .rev()
            .take_while(|o| matches!(o.kind, CigarKind::SoftClip | CigarKind::HardClip))
            .map(|o| o.len as i64)
            .sum();
        let span = rec.cigar.ref_len() as i64;
        rec.pos + span + trailing
    }
}

/// A "quality" key used to pick the representative of a duplicate
/// group: prefer higher summed base quality, then higher MAPQ.
fn quality_key(rec: &SamRecord) -> QualityKey {
    let sum: u64 = if rec.qual.is_empty() {
        0
    } else {
        rec.qual
            .as_bytes()
            .iter()
            .map(|&q| q.saturating_sub(33) as u64)
            .sum()
    };
    (sum, rec.mapq)
}

/// Marks PCR/optical duplicates by alignment coordinate, in place.
///
/// Records are grouped by `(rname, strand, unclipped 5′ position)`. In
/// each group the record with the highest quality key — summed base
/// quality, then MAPQ — keeps its `DUPLICATE` flag clear; every other
/// record in the group has the flag set. Unmapped records are left
/// untouched. The pre-existing `DUPLICATE` bit on every mapped record
/// is cleared first so the pass is idempotent.
pub fn mark_duplicates_coordinate(records: &mut [SamRecord]) -> DedupStats {
    let mut stats = DedupStats {
        total: records.len(),
        ..DedupStats::default()
    };

    // First pass: clear every mapped record's DUPLICATE flag (so the
    // routine is idempotent) and count the unmapped records.
    for rec in records.iter_mut() {
        if rec.is_unmapped() || rec.pos <= 0 {
            stats.unmapped += 1;
        } else {
            rec.flags.set(SamFlags::DUPLICATE, false);
        }
    }

    // Resolve the highest-quality representative of every coordinate
    // group: group key -> (winning index, its quality key).
    let mut group_best: HashMap<GroupKey, (usize, QualityKey)> = HashMap::new();
    for (idx, rec) in records.iter().enumerate() {
        if rec.is_unmapped() || rec.pos <= 0 {
            continue;
        }
        let key = (
            rec.rname.clone(),
            rec.flags.is_reverse(),
            unclipped_5prime(rec),
        );
        let qk = quality_key(rec);
        match group_best.get(&key) {
            Some(&(_, best_qk)) if best_qk >= qk => {}
            _ => {
                group_best.insert(key, (idx, qk));
            }
        }
    }

    // Second pass: flag everyone who is not the group winner.
    for (idx, rec) in records.iter_mut().enumerate() {
        if rec.is_unmapped() || rec.pos <= 0 {
            continue;
        }
        let key = (
            rec.rname.clone(),
            rec.flags.is_reverse(),
            unclipped_5prime(rec),
        );
        let winner = group_best.get(&key).map(|&(i, _)| i);
        if winner == Some(idx) {
            stats.unique += 1;
        } else {
            rec.flags.set(SamFlags::DUPLICATE, true);
            stats.duplicates += 1;
        }
    }
    stats
}

/// Flags reads sharing an identical sequence string.
///
/// Returns a parallel `Vec<bool>` — `true` at index `i` means read `i`
/// is a sequence duplicate of an earlier read. The **first**
/// occurrence of each distinct sequence is the kept representative.
/// Empty-sequence reads are never flagged.
pub fn dedup_sequences(reads: &[&str]) -> Vec<bool> {
    use std::collections::HashSet;
    let mut seen: HashSet<String> = HashSet::new();
    let mut is_dup = Vec::with_capacity(reads.len());
    for &r in reads {
        if r.is_empty() {
            is_dup.push(false);
            continue;
        }
        let upper = r.to_ascii_uppercase();
        is_dup.push(!seen.insert(upper));
    }
    is_dup
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::sam::{Cigar, SamFlags};

    fn aln(name: &str, pos: i64, cigar: &str, seq: &str, qual: &str, rev: bool) -> SamRecord {
        let mut r = SamRecord::unmapped(name);
        r.flags = SamFlags(if rev { SamFlags::REVERSE } else { 0 });
        r.rname = "chr1".to_string();
        r.pos = pos;
        r.mapq = 60;
        r.cigar = Cigar::parse(cigar).unwrap();
        r.seq = seq.to_string();
        r.qual = qual.to_string();
        r
    }

    #[test]
    fn coordinate_marks_same_position_reads() {
        // Three reads at the same forward 5' coordinate.
        let mut recs = vec![
            aln("r1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
            aln("r2", 100, "10M", "ACGTACGTAC", "##########", false),
            aln("r3", 100, "10M", "ACGTACGTAC", "++++++++++", false),
        ];
        let stats = mark_duplicates_coordinate(&mut recs);
        assert_eq!(stats.unique, 1);
        assert_eq!(stats.duplicates, 2);
        // The highest-quality read (r1, all 'I') is the representative.
        assert!(!recs[0].flags.is_duplicate());
        assert!(recs[1].flags.is_duplicate());
        assert!(recs[2].flags.is_duplicate());
    }

    #[test]
    fn different_positions_not_duplicates() {
        let mut recs = vec![
            aln("r1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
            aln("r2", 200, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
        ];
        let stats = mark_duplicates_coordinate(&mut recs);
        assert_eq!(stats.unique, 2);
        assert_eq!(stats.duplicates, 0);
    }

    #[test]
    fn strand_separates_groups() {
        // Same POS but opposite strands -> different fragments.
        let mut recs = vec![
            aln("r1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
            aln("r2", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", true),
        ];
        let stats = mark_duplicates_coordinate(&mut recs);
        assert_eq!(stats.unique, 2);
    }

    #[test]
    fn soft_clip_normalizes_position() {
        // r2 has a 5-base leading soft clip; its unclipped 5' equals
        // r1's POS.
        let mut recs = vec![
            aln("r1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
            aln(
                "r2",
                105,
                "5S10M",
                "AAAAAACGTACGTAC",
                "IIIIIIIIIIIIIII",
                false,
            ),
        ];
        let stats = mark_duplicates_coordinate(&mut recs);
        assert_eq!(stats.duplicates, 1);
    }

    #[test]
    fn unmapped_reads_passed_through() {
        let mut recs = vec![SamRecord::unmapped("r1"), SamRecord::unmapped("r2")];
        let stats = mark_duplicates_coordinate(&mut recs);
        assert_eq!(stats.unmapped, 2);
        assert_eq!(stats.duplicates, 0);
    }

    #[test]
    fn idempotent() {
        let mut recs = vec![
            aln("r1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
            aln("r2", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", false),
        ];
        let s1 = mark_duplicates_coordinate(&mut recs);
        let s2 = mark_duplicates_coordinate(&mut recs);
        assert_eq!(s1, s2);
    }

    #[test]
    fn sequence_dedup_flags_repeats() {
        let reads = vec!["ACGT", "TTTT", "acgt", "GGGG", "ACGT"];
        let dup = dedup_sequences(&reads);
        // index 2 ("acgt" == "ACGT") and index 4 are duplicates.
        assert_eq!(dup, vec![false, false, true, false, true]);
    }

    #[test]
    fn duplicate_rate() {
        let mut recs = vec![
            aln("r1", 100, "5M", "ACGTA", "IIIII", false),
            aln("r2", 100, "5M", "ACGTA", "IIIII", false),
        ];
        let stats = mark_duplicates_coordinate(&mut recs);
        assert!((stats.duplicate_rate() - 0.5).abs() < 1e-9);
    }
}
