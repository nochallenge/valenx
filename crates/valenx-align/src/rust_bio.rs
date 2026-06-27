//! `rust-bio` adapters — battle-tested foundational algorithms.
//!
//! This module wraps the upstream [`bio`](https://crates.io/crates/bio)
//! crate (rust-bio) so its widely-used reference implementations are
//! available alongside the in-house dynamic-programming routines in
//! [`crate::pairwise`] and [`crate::search`].
//!
//! Two capability families are exposed:
//!
//! - **Pairwise alignment** — Needleman-Wunsch (global), Smith-Waterman
//!   (local) and semi-global alignment with an affine gap penalty, via
//!   [`global_align`] / [`local_align`] / [`semiglobal_align`]. Each
//!   returns an [`RbAlignment`] carrying the score, the aligned spans
//!   and a CIGAR string.
//! - **Full-text indexing** — a Burrows-Wheeler-transform / FM-index
//!   built over a suffix array ([`FmIndex`]) for exact substring
//!   search returning every match position, plus direct access to the
//!   underlying suffix array.
//!
//! ## Why both this and the in-house code?
//!
//! The in-house [`crate::pairwise`] and [`crate::search::fmindex`]
//! modules remain the primary, dependency-free implementations. This
//! adapter complements them: rust-bio is a long-standing community
//! reference, so these wrappers double as an independent oracle for
//! cross-validating the native routines (see the tests below, which
//! check the two FM-indexes agree), and they give callers a familiar
//! API surface. Per the crate's "real working v1" philosophy, anything
//! that would conflict with the in-house data model is kept in-house;
//! here we adopt rust-bio cleanly for the algorithms it provides well.
//!
//! ```
//! use valenx_align::rust_bio::{global_align, FmIndex};
//!
//! // Global (Needleman-Wunsch) alignment. The reference has one extra
//! // base, so relative to the query it reads as a deletion.
//! let aln = global_align(b"ACGT", b"ACGGT", 1, -1, -5, -1).unwrap();
//! assert_eq!(aln.cigar(), "2M1D2M");
//!
//! // FM-index exact search.
//! let fm = FmIndex::new(b"ACGTACGTACGT").unwrap();
//! assert_eq!(fm.search(b"ACGT"), vec![0, 4, 8]);
//! ```

use crate::error::{AlignError, Result};

// Disambiguate the external `bio` crate from this crate's own
// `crate::bio` (valenx-bioseq interop) module via the absolute path.
use ::bio::alignment::pairwise::Aligner;
use ::bio::alignment::AlignmentOperation;
use ::bio::alphabets::dna;
use ::bio::data_structures::bwt::{bwt, less, Occ};
use ::bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
use ::bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};

/// A pairwise alignment produced by one of the rust-bio aligners.
///
/// Coordinates follow rust-bio's convention: `x` is the query, `y` is
/// the reference, and the `*start` / `*end` fields are half-open spans
/// into each sequence. For global alignment the spans always cover the
/// full inputs; for local alignment they bound the optimal sub-alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RbAlignment {
    /// Optimal alignment score under the supplied scoring.
    pub score: i32,
    /// Start offset of the alignment in the query (`x`).
    pub xstart: usize,
    /// End offset (exclusive) of the alignment in the query (`x`).
    pub xend: usize,
    /// Start offset of the alignment in the reference (`y`).
    pub ystart: usize,
    /// End offset (exclusive) of the alignment in the reference (`y`).
    pub yend: usize,
    /// The column-by-column edit operations of the alignment.
    pub operations: Vec<AlignmentOperation>,
}

impl RbAlignment {
    /// Renders the alignment as a run-length CIGAR string.
    ///
    /// `Match` and `Subst` both consume one query and one reference
    /// base and collapse into `M` (the standard, non-extended CIGAR
    /// alphabet); `Ins` (a base present in the query but not the
    /// reference) becomes `I`; `Del` becomes `D`; soft clips
    /// (`Xclip` / `Yclip`) become `S`. Consecutive identical
    /// operations are merged, e.g. `4M1I`.
    #[must_use]
    pub fn cigar(&self) -> String {
        let mut out = String::new();
        let mut run: Option<(char, usize)> = None;
        let mut flush = |run: &mut Option<(char, usize)>| {
            if let Some((op, n)) = run.take() {
                out.push_str(&n.to_string());
                out.push(op);
            }
        };
        for op in &self.operations {
            let c = match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => 'M',
                AlignmentOperation::Ins => 'I',
                AlignmentOperation::Del => 'D',
                AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => 'S',
            };
            match &mut run {
                Some((cur, n)) if *cur == c => *n += 1,
                _ => {
                    flush(&mut run);
                    run = Some((c, 1));
                }
            }
        }
        flush(&mut run);
        out
    }

    /// Number of `Match` (identical) columns in the alignment.
    #[must_use]
    pub fn matches(&self) -> usize {
        self.operations
            .iter()
            .filter(|o| matches!(o, AlignmentOperation::Match))
            .count()
    }
}

/// Builds an [`Aligner`] and runs the supplied alignment mode.
///
/// The scoring is the classic match/mismatch model with an affine gap
/// penalty: a base pair scores `match_score` if equal else
/// `mismatch_score`; opening a gap costs `gap_open` and each further
/// gap base costs `gap_extend`. By rust-bio's convention the penalties
/// are passed as non-positive integers (e.g. `gap_open = -5`).
fn run_align(
    x: &[u8],
    y: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    mode: Mode,
) -> Result<RbAlignment> {
    if x.is_empty() || y.is_empty() {
        return Err(AlignError::invalid(
            "seq",
            "both sequences must be non-empty",
        ));
    }
    let score = |a: u8, b: u8| if a == b { match_score } else { mismatch_score };
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score);
    let a = match mode {
        Mode::Global => aligner.global(x, y),
        Mode::Local => aligner.local(x, y),
        Mode::Semiglobal => aligner.semiglobal(x, y),
    };
    Ok(RbAlignment {
        score: a.score,
        xstart: a.xstart,
        xend: a.xend,
        ystart: a.ystart,
        yend: a.yend,
        operations: a.operations,
    })
}

#[derive(Clone, Copy)]
enum Mode {
    Global,
    Local,
    Semiglobal,
}

/// Global (Needleman-Wunsch) alignment of query `x` against reference
/// `y` with an affine gap penalty.
///
/// Returns [`AlignError::Invalid`] if either sequence is empty.
pub fn global_align(
    x: &[u8],
    y: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<RbAlignment> {
    run_align(
        x,
        y,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
        Mode::Global,
    )
}

/// Local (Smith-Waterman) alignment: the highest-scoring sub-alignment
/// of `x` and `y` under the given affine scoring.
///
/// Returns [`AlignError::Invalid`] if either sequence is empty.
pub fn local_align(
    x: &[u8],
    y: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<RbAlignment> {
    run_align(
        x,
        y,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
        Mode::Local,
    )
}

/// Semi-global alignment: global in the query `x`, free end gaps in the
/// reference `y` (useful for fitting a short read into a longer
/// reference).
///
/// Returns [`AlignError::Invalid`] if either sequence is empty.
pub fn semiglobal_align(
    x: &[u8],
    y: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<RbAlignment> {
    run_align(
        x,
        y,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
        Mode::Semiglobal,
    )
}

/// An FM-index (Burrows-Wheeler transform over a suffix array) for
/// exact substring search, backed by rust-bio.
///
/// Construction runs in `O(n)` time and stores the suffix array, the
/// BWT, the `less` table and an `Occ` (rank) table; queries then run in
/// time proportional to the pattern length plus the number of hits.
///
/// A `$` sentinel (which must not occur in the input) is appended
/// internally as rust-bio's BWT requires; it is transparent to callers
/// and never reported as a match.
pub struct FmIndex {
    fm: FMIndex<Vec<u8>, Vec<usize>, Occ>,
    sa: RawSuffixArray,
    text_len: usize,
}

impl FmIndex {
    /// Builds an FM-index over `text`.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::Invalid`] if `text` is empty or already
    /// contains the reserved `$` sentinel byte.
    pub fn new(text: &[u8]) -> Result<Self> {
        if text.is_empty() {
            return Err(AlignError::invalid("text", "text must be non-empty"));
        }
        if text.contains(&b'$') {
            return Err(AlignError::invalid(
                "text",
                "text must not contain the reserved '$' sentinel byte",
            ));
        }
        let mut buf = Vec::with_capacity(text.len() + 1);
        buf.extend_from_slice(text);
        buf.push(b'$');

        // An IUPAC DNA alphabet plus the sentinel covers nucleotide
        // input; for arbitrary bytes we widen to every symbol present.
        let mut alphabet = dna::iupac_alphabet();
        alphabet.insert(b'$');
        for &b in text {
            alphabet.insert(b);
        }

        let sa = suffix_array(&buf);
        let bwt_vec = bwt(&buf, &sa);
        let less_vec = less(&bwt_vec, &alphabet);
        let occ = Occ::new(&bwt_vec, 3, &alphabet);
        let fm = FMIndex::new(bwt_vec, less_vec, occ);
        Ok(Self {
            fm,
            sa,
            text_len: text.len(),
        })
    }

    /// Returns every start position in the original text where
    /// `pattern` occurs exactly, in ascending order.
    ///
    /// An empty pattern, or one not present, yields an empty vector.
    #[must_use]
    pub fn search(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() {
            return Vec::new();
        }
        match self.fm.backward_search(pattern.iter()) {
            BackwardSearchResult::Complete(interval) => {
                let mut hits: Vec<usize> = interval
                    .occ(&self.sa)
                    .into_iter()
                    // Drop the sentinel position if it ever surfaces.
                    .filter(|&p| p < self.text_len)
                    .collect();
                hits.sort_unstable();
                hits
            }
            BackwardSearchResult::Partial(..) | BackwardSearchResult::Absent => Vec::new(),
        }
    }

    /// Number of exact occurrences of `pattern` in the text.
    #[must_use]
    pub fn count(&self, pattern: &[u8]) -> usize {
        self.search(pattern).len()
    }

    /// `true` if `pattern` occurs at least once in the text.
    #[must_use]
    pub fn contains(&self, pattern: &[u8]) -> bool {
        !self.search(pattern).is_empty()
    }

    /// Length of the indexed text (excluding the internal sentinel).
    #[must_use]
    pub fn text_len(&self) -> usize {
        self.text_len
    }

    /// The suffix array of the sentinel-terminated text.
    ///
    /// Index `0` is always the sentinel suffix; the remaining entries
    /// are the suffix start positions of the original text in
    /// lexicographic order.
    #[must_use]
    pub fn suffix_array(&self) -> &[usize] {
        &self.sa
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Pairwise alignment ------------------------------------------

    #[test]
    fn global_acgt_vs_acggt_scores_and_cigar() {
        // ACGT (query) vs ACGGT (reference): the reference carries one
        // extra G, which relative to the query is a single-base
        // deletion. 4 matches @ +1 each, one gap-open @ -5 ->
        // score 4 - 5 = -1. rust-bio places the gap mid-sequence,
        // giving the CIGAR 2M1D2M.
        let aln = global_align(b"ACGT", b"ACGGT", 1, -1, -5, -1).unwrap();
        assert_eq!(aln.score, -1);
        assert_eq!(aln.cigar(), "2M1D2M");
        assert_eq!(aln.matches(), 4);
        // Global alignment spans the whole query and reference.
        assert_eq!((aln.xstart, aln.xend), (0, 4));
        assert_eq!((aln.ystart, aln.yend), (0, 5));
    }

    #[test]
    fn global_identical_sequences_full_match() {
        let aln = global_align(b"ACGTACGT", b"ACGTACGT", 1, -1, -5, -1).unwrap();
        assert_eq!(aln.score, 8);
        assert_eq!(aln.cigar(), "8M");
        assert_eq!(aln.matches(), 8);
    }

    #[test]
    fn local_finds_embedded_exact_substring() {
        // The query AAAA sits inside the reference; local alignment
        // should recover it with a perfect score of 4.
        let aln = local_align(b"AAAA", b"GGGGAAAAGGGG", 1, -1, -5, -1).unwrap();
        assert_eq!(aln.score, 4);
        assert_eq!(aln.cigar(), "4M");
        assert_eq!((aln.ystart, aln.yend), (4, 8));
    }

    #[test]
    fn semiglobal_fits_query_into_reference() {
        let aln = semiglobal_align(b"ACGT", b"TTACGTTT", 1, -1, -5, -1).unwrap();
        assert_eq!(aln.score, 4);
        assert_eq!(aln.matches(), 4);
    }

    #[test]
    fn empty_input_is_rejected() {
        assert!(global_align(b"", b"ACGT", 1, -1, -5, -1).is_err());
        assert!(local_align(b"ACGT", b"", 1, -1, -5, -1).is_err());
    }

    // --- FM-index search ---------------------------------------------

    #[test]
    fn fmindex_search_returns_all_positions() {
        let fm = FmIndex::new(b"ACGTACGTACGT").unwrap();
        // ACGT occurs at offsets 0, 4 and 8.
        assert_eq!(fm.search(b"ACGT"), vec![0, 4, 8]);
        assert_eq!(fm.count(b"ACGT"), 3);
        assert!(fm.contains(b"GTAC"));
        assert_eq!(fm.search(b"GTAC"), vec![2, 6]);
    }

    #[test]
    fn fmindex_single_and_full_and_absent() {
        let fm = FmIndex::new(b"MISSISSIPPI").unwrap();
        // Classic suffix-structure example: 'ISS' occurs twice.
        assert_eq!(fm.search(b"ISS"), vec![1, 4]);
        // 'SSI' occurs twice.
        assert_eq!(fm.search(b"SSI"), vec![2, 5]);
        // Whole-text match is a single hit at 0.
        assert_eq!(fm.search(b"MISSISSIPPI"), vec![0]);
        // Absent pattern.
        assert!(fm.search(b"XYZ").is_empty());
        // Empty pattern.
        assert!(fm.search(b"").is_empty());
        assert_eq!(fm.text_len(), 11);
    }

    #[test]
    fn fmindex_rejects_empty_and_sentinel_text() {
        assert!(FmIndex::new(b"").is_err());
        assert!(FmIndex::new(b"AC$GT").is_err());
    }

    /// Cross-validate the rust-bio FM-index against the in-house one:
    /// for a battery of patterns they must report identical position
    /// sets, which independently exercises both implementations.
    #[test]
    fn fmindex_agrees_with_in_house() {
        let text = b"ACGTACGTTTACGTACACGT";
        let rb = FmIndex::new(text).unwrap();
        let native = crate::search::fmindex::FmIndex::build(text).unwrap();
        for pat in [
            &b"A"[..],
            b"AC",
            b"ACGT",
            b"CGT",
            b"TTT",
            b"ACAC",
            b"GG", // absent
            text,
        ] {
            let mut want = native.locate(pat);
            want.sort_unstable();
            assert_eq!(rb.search(pat), want, "mismatch for {pat:?}");
        }
    }
}
