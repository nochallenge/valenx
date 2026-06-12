//! CRISPR edit-outcome analysis (CRISPResso-class).
//!
//! After a CRISPR experiment, amplicon reads spanning the cut site are
//! sequenced and compared against the unedited reference amplicon.
//! CRISPResso aligns each read to the reference, classifies the
//! outcome (unmodified, insertion, deletion, substitution), quantifies
//! the indel-size spectrum, and reports the **editing efficiency** —
//! the fraction of reads bearing an indel.
//!
//! This module implements that pipeline on top of `valenx-align`'s
//! global pairwise alignment ([`valenx_align::needleman_wunsch`]):
//!
//! - [`analyze_amplicon`] — aligns every edited read to the reference
//!   amplicon and produces an [`EditOutcomeReport`];
//! - the report carries the per-class read counts, the indel-size
//!   histogram, the editing efficiency and an optional
//!   quantification-window restriction (count indels only if they
//!   overlap a window around the predicted cut site — the CRISPResso
//!   `--quantification_window_size` behaviour).
//!
//! ## v1 scope
//!
//! A read is classified by the *net* indels in its global alignment to
//! the reference. The analysis is single-amplicon and does not do
//! CRISPResso's allele-frequency table, prime-editing / base-editing
//! sub-modes, or HDR-template quantification — it is the core
//! indel-spectrum quantifier. Substitutions are counted but, like
//! CRISPResso's default, do not by themselves mark a read as "edited".

use crate::error::{GenomicsError, Result};
use std::collections::BTreeMap;
use valenx_align::pairwise::global::gotoh;
use valenx_align::ScoringScheme;

/// The edit class of a single read.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum EditClass {
    /// The read matches the reference (no indel; substitutions allowed).
    Unmodified,
    /// The read carries a net insertion relative to the reference.
    Insertion,
    /// The read carries a net deletion relative to the reference.
    Deletion,
    /// The read carries both insertion and deletion events.
    InsertionDeletion,
}

impl EditClass {
    /// `true` for any class bearing an indel — i.e. an "edited" read.
    pub fn is_edited(self) -> bool {
        !matches!(self, EditClass::Unmodified)
    }
}

/// The per-read alignment outcome.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadOutcome {
    /// The edit class.
    pub class: EditClass,
    /// Total inserted bases (relative to the reference).
    pub inserted: usize,
    /// Total deleted bases (relative to the reference).
    pub deleted: usize,
    /// Number of mismatched (substituted) aligned columns.
    pub substitutions: usize,
    /// Net length change `inserted - deleted` (signed).
    pub net_indel: i64,
}

/// Parameters controlling the edit-outcome analysis.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct EditAnalysisParams {
    /// 0-based reference position of the predicted cut site (typically
    /// 3 bp 5′ of the PAM). Used to centre the quantification window.
    pub cut_site: usize,
    /// Half-width of the quantification window in bases. An indel is
    /// only counted as an edit when it overlaps `[cut_site -
    /// window, cut_site + window]`. Set to `usize::MAX` to count indels
    /// anywhere in the amplicon.
    pub window: usize,
}

impl Default for EditAnalysisParams {
    /// Count indels anywhere (no quantification-window restriction).
    fn default() -> Self {
        EditAnalysisParams {
            cut_site: 0,
            window: usize::MAX,
        }
    }
}

/// The full CRISPResso-style report for an amplicon experiment.
#[derive(Clone, Debug, PartialEq)]
pub struct EditOutcomeReport {
    /// Total reads analysed.
    pub total_reads: usize,
    /// Reads classified as unmodified.
    pub unmodified: usize,
    /// Reads with an insertion only.
    pub insertions: usize,
    /// Reads with a deletion only.
    pub deletions: usize,
    /// Reads with both an insertion and a deletion.
    pub indels: usize,
    /// Indel-size spectrum: `net length change -> read count`. A
    /// negative key is a net deletion, a positive key a net insertion.
    pub size_spectrum: BTreeMap<i64, usize>,
    /// Per-read outcomes, in input order.
    pub outcomes: Vec<ReadOutcome>,
}

impl EditOutcomeReport {
    /// Editing efficiency — the fraction of reads that carry an indel.
    pub fn editing_efficiency(&self) -> f64 {
        if self.total_reads == 0 {
            return 0.0;
        }
        let edited = self.insertions + self.deletions + self.indels;
        edited as f64 / self.total_reads as f64
    }

    /// The most common net indel size (`0` means most reads are
    /// unmodified). Returns `None` for an empty report.
    pub fn modal_indel_size(&self) -> Option<i64> {
        self.size_spectrum
            .iter()
            .max_by_key(|(_, &count)| count)
            .map(|(&size, _)| size)
    }

    /// A compact text summary.
    pub fn summary_text(&self) -> String {
        format!(
            "reads={} unmodified={} ins={} del={} indel={} efficiency={:.1}%",
            self.total_reads,
            self.unmodified,
            self.insertions,
            self.deletions,
            self.indels,
            self.editing_efficiency() * 100.0,
        )
    }
}

/// Classifies one aligned read against the reference.
///
/// `ref_row` and `read_row` are the gapped rows of a global alignment
/// (same length; `-` for gaps). The routine walks the columns counting
/// insertions (gap in the reference), deletions (gap in the read) and
/// substitutions, and — when a quantification window is set — only
/// counts an indel column that falls inside the window.
fn classify_aligned(ref_row: &[u8], read_row: &[u8], params: &EditAnalysisParams) -> ReadOutcome {
    let mut inserted = 0usize;
    let mut deleted = 0usize;
    let mut substitutions = 0usize;
    let mut ref_pos = 0usize; // 0-based position in the ungapped reference

    let window_lo = params.cut_site.saturating_sub(params.window);
    let window_hi = params.cut_site.saturating_add(params.window);
    let in_window =
        |p: usize| -> bool { params.window == usize::MAX || (p >= window_lo && p <= window_hi) };

    for (&r, &q) in ref_row.iter().zip(read_row) {
        match (r, q) {
            (b'-', _) => {
                // Insertion relative to the reference. Anchored at the
                // current reference position.
                if in_window(ref_pos) {
                    inserted += 1;
                }
            }
            (_, b'-') => {
                // Deletion from the read.
                if in_window(ref_pos) {
                    deleted += 1;
                }
                ref_pos += 1;
            }
            (rr, qq) => {
                if !rr.eq_ignore_ascii_case(&qq) {
                    substitutions += 1;
                }
                ref_pos += 1;
            }
        }
    }

    let class = match (inserted > 0, deleted > 0) {
        (false, false) => EditClass::Unmodified,
        (true, false) => EditClass::Insertion,
        (false, true) => EditClass::Deletion,
        (true, true) => EditClass::InsertionDeletion,
    };
    ReadOutcome {
        class,
        inserted,
        deleted,
        substitutions,
        net_indel: inserted as i64 - deleted as i64,
    }
}

/// Analyses an amplicon CRISPR experiment.
///
/// `reference` is the unedited amplicon; `reads` are the sequenced
/// edited-amplicon reads. Each read is globally aligned to the
/// reference with a DNA scoring scheme; the outcomes are tallied into
/// an [`EditOutcomeReport`].
///
/// Returns [`GenomicsError::Invalid`] for an empty reference.
pub fn analyze_amplicon(
    reference: &[u8],
    reads: &[&[u8]],
    params: &EditAnalysisParams,
) -> Result<EditOutcomeReport> {
    if reference.is_empty() {
        return Err(GenomicsError::invalid("reference", "amplicon is empty"));
    }
    let scheme = ScoringScheme::dna_default();
    let mut report = EditOutcomeReport {
        total_reads: reads.len(),
        unmodified: 0,
        insertions: 0,
        deletions: 0,
        indels: 0,
        size_spectrum: BTreeMap::new(),
        outcomes: Vec::with_capacity(reads.len()),
    };

    for &read in reads {
        if read.is_empty() {
            // An empty read is a total deletion of the amplicon.
            let outcome = ReadOutcome {
                class: EditClass::Deletion,
                inserted: 0,
                deleted: reference.len(),
                substitutions: 0,
                net_indel: -(reference.len() as i64),
            };
            tally(&mut report, &outcome);
            continue;
        }
        // Align read (row1) to reference (row2) with an *affine* gap
        // model. A linear gap penalty makes two cheap single-base gaps
        // (a delete + an insert) outscore one mismatch, so a plain
        // substitution would be mis-represented as an indel and
        // wrongly counted as an edit. Affine gaps (open + extend) keep
        // a real NHEJ indel as one event and leave substitutions as
        // substitutions.
        let aln = gotoh(read, reference, &scheme)
            .map_err(|e| GenomicsError::invalid("alignment", e.to_string()))?;
        // row1 = read, row2 = reference.
        let outcome = classify_aligned(&aln.row2, &aln.row1, params);
        tally(&mut report, &outcome);
    }
    Ok(report)
}

/// Folds one read outcome into the running report.
fn tally(report: &mut EditOutcomeReport, outcome: &ReadOutcome) {
    match outcome.class {
        EditClass::Unmodified => report.unmodified += 1,
        EditClass::Insertion => report.insertions += 1,
        EditClass::Deletion => report.deletions += 1,
        EditClass::InsertionDeletion => report.indels += 1,
    }
    *report.size_spectrum.entry(outcome.net_indel).or_insert(0) += 1;
    report.outcomes.push(outcome.clone());
}

#[cfg(test)]
mod tests {
    use super::*;

    const AMPLICON: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGT";

    #[test]
    fn unmodified_read_classified() {
        let reads: Vec<&[u8]> = vec![AMPLICON, AMPLICON];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.unmodified, 2);
        assert_eq!(r.editing_efficiency(), 0.0);
    }

    #[test]
    fn deletion_read_classified() {
        // Drop 4 bases from the middle of the amplicon.
        let mut edited = AMPLICON[..14].to_vec();
        edited.extend_from_slice(&AMPLICON[18..]);
        let reads: Vec<&[u8]> = vec![&edited];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.deletions, 1);
        assert_eq!(r.outcomes[0].deleted, 4);
        assert_eq!(r.outcomes[0].net_indel, -4);
        assert!((r.editing_efficiency() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn insertion_read_classified() {
        // Insert "TTTT" into the amplicon.
        let mut edited = AMPLICON[..16].to_vec();
        edited.extend_from_slice(b"TTTT");
        edited.extend_from_slice(&AMPLICON[16..]);
        let reads: Vec<&[u8]> = vec![&edited];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.insertions, 1);
        assert_eq!(r.outcomes[0].inserted, 4);
        assert_eq!(r.outcomes[0].net_indel, 4);
    }

    #[test]
    fn mixed_pool_efficiency() {
        let mut del = AMPLICON[..14].to_vec();
        del.extend_from_slice(&AMPLICON[17..]);
        // 2 unmodified, 1 deletion, 1 insertion-laden.
        let mut ins = AMPLICON[..16].to_vec();
        ins.extend_from_slice(b"GG");
        ins.extend_from_slice(&AMPLICON[16..]);
        let reads: Vec<&[u8]> = vec![AMPLICON, AMPLICON, &del, &ins];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.total_reads, 4);
        assert_eq!(r.unmodified, 2);
        // Editing efficiency = 2/4 = 0.5.
        assert!((r.editing_efficiency() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn size_spectrum_built() {
        let mut del = AMPLICON[..14].to_vec();
        del.extend_from_slice(&AMPLICON[17..]); // -3
        let reads: Vec<&[u8]> = vec![AMPLICON, &del, &del];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.size_spectrum.get(&0), Some(&1));
        assert_eq!(r.size_spectrum.get(&-3), Some(&2));
        assert_eq!(r.modal_indel_size(), Some(-3));
    }

    #[test]
    fn substitution_does_not_count_as_edit() {
        // One substituted base, no indel.
        let mut sub = AMPLICON.to_vec();
        sub[15] = if sub[15] == b'A' { b'C' } else { b'A' };
        let reads: Vec<&[u8]> = vec![&sub];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.unmodified, 1);
        assert_eq!(r.outcomes[0].substitutions, 1);
        assert_eq!(r.editing_efficiency(), 0.0);
    }

    #[test]
    fn quantification_window_restricts_counting() {
        // A deletion far from the cut site is ignored when the window
        // is tight.
        let mut del = AMPLICON[..2].to_vec(); // delete bases 2..5 (far 5')
        del.extend_from_slice(&AMPLICON[5..]);
        let reads: Vec<&[u8]> = vec![&del];
        // Cut site at position 28, window 3 -> the 5' deletion is out.
        let params = EditAnalysisParams {
            cut_site: 28,
            window: 3,
        };
        let r = analyze_amplicon(AMPLICON, &reads, &params).unwrap();
        assert_eq!(r.unmodified, 1, "5' deletion should be outside the window");
    }

    #[test]
    fn empty_read_is_total_deletion() {
        let empty: &[u8] = b"";
        let reads: Vec<&[u8]> = vec![empty];
        let r = analyze_amplicon(AMPLICON, &reads, &EditAnalysisParams::default()).unwrap();
        assert_eq!(r.deletions, 1);
        assert_eq!(r.outcomes[0].deleted, AMPLICON.len());
    }

    #[test]
    fn rejects_empty_reference() {
        let reads: Vec<&[u8]> = vec![b"ACGT"];
        assert!(analyze_amplicon(b"", &reads, &EditAnalysisParams::default()).is_err());
    }
}
