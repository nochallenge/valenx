//! Global pairwise alignment — Needleman-Wunsch and Gotoh.
//!
//! Two routines, both end-to-end (every residue of both sequences is
//! placed):
//!
//! - [`needleman_wunsch`] — the classic O(nm) global DP with a *linear*
//!   gap penalty (every gap residue costs the same).
//! - [`gotoh`] — global DP with an *affine* gap penalty (open + extend)
//!   via Gotoh's three-matrix recurrence, still O(nm) time.
//!
//! Both return an [`Alignment`] with `span1 = (0, len1)` and
//! `span2 = (0, len2)` — a global alignment spans both inputs fully.

use super::hirschberg::hirschberg;
use super::result::Alignment;
use crate::error::Result;
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};
use crate::matrix::{GapCost, ScoringScheme};

/// A very negative sentinel used as `-inf` in the DP without risking
/// `i32` overflow when a finite penalty is added to it.
const NEG_INF: i32 = i32::MIN / 4;

/// Needleman-Wunsch global alignment with a **linear** gap penalty.
///
/// The gap penalty used is `scheme.gap.extend` per gap residue; the
/// `open` field is ignored (use [`gotoh`] for affine gaps). Runs in
/// O(`a.len()` · `b.len()`) time and space.
///
/// When the `(n+1)·(m+1)` score matrix would exceed
/// [`MAX_DP_CELLS`], this transparently
/// falls back to the linear-space
/// [`hirschberg`] routine — it
/// uses the same linear gap penalty and yields the identical optimal
/// global alignment in O(min(n, m)) space, so the caller never has to
/// handle an error for a merely-large input.
pub fn needleman_wunsch(a: &[u8], b: &[u8], scheme: &ScoringScheme) -> Result<Alignment> {
    needleman_wunsch_capped(a, b, scheme, MAX_DP_CELLS)
}

/// [`needleman_wunsch`] with an explicit cell cap (so tests can force
/// the linear-space fallback without allocating a production-size
/// matrix). Over the cap it dispatches to [`hirschberg`], which
/// produces the same linear-gap optimum.
fn needleman_wunsch_capped(
    a: &[u8],
    b: &[u8],
    scheme: &ScoringScheme,
    max_cells: usize,
) -> Result<Alignment> {
    let n = a.len();
    let m = b.len();
    let g = scheme.gap.extend;

    // Bound the quadratic allocation: over the cap, the linear-space
    // Hirschberg routine gives the same result without the O(nm) matrix.
    if check_dp_size_with(n + 1, m + 1, max_cells).is_err() {
        return hirschberg(a, b, scheme);
    }

    // Score matrix (n+1) x (m+1) and a traceback matrix.
    let w = m + 1;
    let mut s = vec![0i32; (n + 1) * w];
    // 0 = diag, 1 = up (gap in b), 2 = left (gap in a)
    let mut tb = vec![0u8; (n + 1) * w];

    for j in 1..=m {
        s[j] = s[j - 1] - g;
        tb[j] = 2;
    }
    for i in 1..=n {
        s[i * w] = s[(i - 1) * w] - g;
        tb[i * w] = 1;
    }

    for i in 1..=n {
        for j in 1..=m {
            let diag = s[(i - 1) * w + j - 1] + scheme.sub(a[i - 1], b[j - 1]);
            let up = s[(i - 1) * w + j] - g;
            let left = s[i * w + j - 1] - g;
            let (best, op) = max3(diag, up, left);
            s[i * w + j] = best;
            tb[i * w + j] = op;
        }
    }

    let score = s[n * w + m];
    let (row1, row2) = traceback_linear(a, b, &tb, w, n, m);
    Alignment::new(row1, row2, score, (0, n), (0, m))
}

/// Gotoh global alignment with an **affine** gap penalty.
///
/// Uses the standard three-matrix recurrence: `m` (ending in a
/// match/mismatch), `ix` (ending in a gap in `b`), `iy` (ending in a
/// gap in `a`). A gap of length `L` costs `open + extend * L`. Runs in
/// O(`a.len()` · `b.len()`) time and space.
///
/// Returns [`AlignError::TooLarge`](crate::error::AlignError::TooLarge)
/// when the `(n+1)·(m+1)` matrix would exceed
/// [`MAX_DP_CELLS`]. Unlike
/// [`needleman_wunsch`] there is no linear-space fallback: Hirschberg
/// models a *linear* gap and would silently change this routine's
/// affine scoring, so an oversized affine alignment is rejected rather
/// than answered with a different model. Use a banded affine routine
/// ([`crate::pairwise::banded::banded_affine`]) for long near-diagonal
/// inputs.
pub fn gotoh(a: &[u8], b: &[u8], scheme: &ScoringScheme) -> Result<Alignment> {
    gotoh_capped(a, b, scheme, MAX_DP_CELLS)
}

/// [`gotoh`] with an explicit cell cap so tests can exercise the
/// over-cap rejection without a production-size allocation.
fn gotoh_capped(a: &[u8], b: &[u8], scheme: &ScoringScheme, max_cells: usize) -> Result<Alignment> {
    let n = a.len();
    let m = b.len();
    let GapCost { open, extend } = scheme.gap;
    let w = m + 1;

    // Six i32 / u8 matrices of (n+1)·(m+1) — bound the allocation.
    check_dp_size_with(n + 1, m + 1, max_cells)?;

    let mut mm = vec![NEG_INF; (n + 1) * w]; // ends in (mis)match
    let mut ix = vec![NEG_INF; (n + 1) * w]; // ends in gap in b (consume a)
    let mut iy = vec![NEG_INF; (n + 1) * w]; // ends in gap in a (consume b)

    // Traceback: which matrix each cell came from, per layer.
    // 0 = from mm, 1 = from ix, 2 = from iy.
    let mut tm = vec![0u8; (n + 1) * w];
    let mut tx = vec![0u8; (n + 1) * w];
    let mut ty = vec![0u8; (n + 1) * w];

    mm[0] = 0;
    for j in 1..=m {
        iy[j] = super::affine_init(open, extend, j, NEG_INF);
        ty[j] = 2; // extend within iy
    }
    for i in 1..=n {
        ix[i * w] = super::affine_init(open, extend, i, NEG_INF);
        tx[i * w] = 1;
    }

    for i in 1..=n {
        for j in 1..=m {
            let idx = i * w + j;
            let sub = scheme.sub(a[i - 1], b[j - 1]);

            // mm[i][j] = best previous cell + substitution score.
            let pdiag = i.saturating_sub(1) * w + j - 1;
            let (best_prev, src) = max3(mm[pdiag], ix[pdiag], iy[pdiag]);
            mm[idx] = best_prev + sub;
            tm[idx] = src;

            // ix: gap in b => consume a, came from row i-1.
            let up = (i - 1) * w + j;
            let open_x = mm[up] - open - extend;
            let ext_x = ix[up] - extend;
            if open_x >= ext_x {
                ix[idx] = open_x;
                tx[idx] = 0;
            } else {
                ix[idx] = ext_x;
                tx[idx] = 1;
            }

            // iy: gap in a => consume b, came from col j-1.
            let left = i * w + j - 1;
            let open_y = mm[left] - open - extend;
            let ext_y = iy[left] - extend;
            if open_y >= ext_y {
                iy[idx] = open_y;
                ty[idx] = 0;
            } else {
                iy[idx] = ext_y;
                ty[idx] = 2;
            }
        }
    }

    let last = n * w + m;
    let (score, start_layer) = max3_score(mm[last], ix[last], iy[last]);

    // Traceback across the three layers.
    let mut row1 = Vec::new();
    let mut row2 = Vec::new();
    let (mut i, mut j) = (n, m);
    let mut layer = start_layer;
    while i > 0 || j > 0 {
        let idx = i * w + j;
        match layer {
            0 => {
                // came from a (mis)match column
                row1.push(a[i - 1]);
                row2.push(b[j - 1]);
                let next = tm[idx];
                i -= 1;
                j -= 1;
                layer = next;
            }
            1 => {
                // gap in b
                row1.push(a[i - 1]);
                row2.push(b'-');
                let next = tx[idx];
                i -= 1;
                layer = next;
            }
            _ => {
                // gap in a
                row1.push(b'-');
                row2.push(b[j - 1]);
                let next = ty[idx];
                j -= 1;
                layer = next;
            }
        }
    }
    row1.reverse();
    row2.reverse();
    Alignment::new(row1, row2, score, (0, n), (0, m))
}

/// Returns `(max, op)` where `op` is 0 for `diag`, 1 for `up`, 2 for
/// `left` — biased toward diagonal then up on ties (stable traceback).
fn max3(diag: i32, up: i32, left: i32) -> (i32, u8) {
    let mut best = diag;
    let mut op = 0u8;
    if up > best {
        best = up;
        op = 1;
    }
    if left > best {
        best = left;
        op = 2;
    }
    (best, op)
}

/// Like [`max3`] but the index meaning is "layer" (0/1/2) for Gotoh.
fn max3_score(a: i32, b: i32, c: i32) -> (i32, u8) {
    max3(a, b, c)
}

/// Walks the linear-gap traceback matrix back to the origin.
fn traceback_linear(
    a: &[u8],
    b: &[u8],
    tb: &[u8],
    w: usize,
    mut i: usize,
    mut j: usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut row1 = Vec::new();
    let mut row2 = Vec::new();
    while i > 0 || j > 0 {
        match tb[i * w + j] {
            0 => {
                row1.push(a[i - 1]);
                row2.push(b[j - 1]);
                i -= 1;
                j -= 1;
            }
            1 => {
                row1.push(a[i - 1]);
                row2.push(b'-');
                i -= 1;
            }
            _ => {
                row1.push(b'-');
                row2.push(b[j - 1]);
                j -= 1;
            }
        }
    }
    row1.reverse();
    row2.reverse();
    (row1, row2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn dna_scheme(open: i32, extend: i32) -> ScoringScheme {
        ScoringScheme::new(
            SubstitutionMatrix::dna_simple(1, -1),
            GapCost::new(open, extend),
        )
    }

    #[test]
    fn nw_identical_sequences() {
        let s = dna_scheme(0, 1);
        let al = needleman_wunsch(b"ACGTACGT", b"ACGTACGT", &s).unwrap();
        assert_eq!(al.score, 8);
        assert_eq!(al.row1_str(), "ACGTACGT");
        assert_eq!(al.row2_str(), "ACGTACGT");
        assert!((al.percent_identity() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn nw_single_gap() {
        // ACGT vs ACT — one deletion of G
        let s = dna_scheme(0, 1);
        let al = needleman_wunsch(b"ACGT", b"ACT", &s).unwrap();
        // 3 matches (+3), 1 gap (-1) => 2
        assert_eq!(al.score, 2);
        assert_eq!(al.len(), 4);
        assert!(al.row2.contains(&b'-'));
    }

    #[test]
    fn nw_span_covers_both() {
        let s = dna_scheme(0, 1);
        let al = needleman_wunsch(b"AAAA", b"AAAAAA", &s).unwrap();
        assert_eq!(al.span1, (0, 4));
        assert_eq!(al.span2, (0, 6));
    }

    #[test]
    fn gotoh_affine_prefers_one_long_gap() {
        // Reference: A A A G G G G A A A  vs query missing the GGGG.
        // With a steep open penalty, affine should keep the gap
        // contiguous (one open) rather than splitting it.
        let a = b"AAAGGGGAAA";
        let b = b"AAAAAA";
        let s = dna_scheme(10, 1);
        let al = gotoh(a, b, &s).unwrap();
        // One contiguous 4-gap: 6 matches (+6) - (open10 + 4*1) = -8
        assert_eq!(al.score, -8);
        // Exactly one run of gaps in row2.
        let gap_runs = al
            .row2
            .split(|&c| c != b'-')
            .filter(|r| !r.is_empty())
            .count();
        assert_eq!(gap_runs, 1, "affine gap should stay contiguous");
    }

    #[test]
    fn gotoh_identical_is_full_score() {
        let s = ScoringScheme::new(SubstitutionMatrix::blosum62(), GapCost::new(11, 1));
        let al = gotoh(b"MKVLAAG", b"MKVLAAG", &s).unwrap();
        let expect: i32 = b"MKVLAAG".iter().map(|&r| s.sub(r, r)).sum();
        assert_eq!(al.score, expect);
    }

    #[test]
    fn gotoh_matches_nw_when_linear() {
        // With open=0 the affine model degenerates to linear; the two
        // routines must agree on the score.
        let a = b"ACGTACGTTT";
        let b = b"ACGTCGTT";
        let lin = dna_scheme(0, 1);
        let nw = needleman_wunsch(a, b, &lin).unwrap();
        let go = gotoh(a, b, &lin).unwrap();
        assert_eq!(nw.score, go.score);
    }

    #[test]
    fn nw_equals_negative_levenshtein_under_unit_costs() {
        // Pins the NW <-> edit-distance contract: with match 0,
        // mismatch -1, gap -1 every alignment column scores 0 (match)
        // or -1 (substitution / indel), so the optimal NW score is
        // exactly the negated minimum edit-operation count, i.e.
        // -levenshtein(a, b). The IDENTITY matrix scores identical
        // residues 0 and any differing pair -1 over the full A..Z
        // alphabet (dna_simple would mis-score letters outside ACGT).
        use crate::util::editdist::levenshtein;
        let s = ScoringScheme::new(SubstitutionMatrix::identity(0, -1), GapCost::new(0, 1));

        // Classic substitutions-plus-insertion example: edit dist 3.
        let kit = needleman_wunsch(b"kitten", b"sitting", &s).unwrap();
        assert_eq!(kit.score, -3, "kitten/sitting edit distance is 3");
        assert_eq!(kit.score, -(levenshtein(b"kitten", b"sitting") as i32));

        // A pair dominated by an indel.
        let ins = needleman_wunsch(b"GATTACA", b"GATTTACA", &s).unwrap();
        assert_eq!(ins.score, -(levenshtein(b"GATTACA", b"GATTTACA") as i32));

        // A mixed substitution + indel pair.
        let mix = needleman_wunsch(b"flaw", b"lawn", &s).unwrap();
        assert_eq!(mix.score, -(levenshtein(b"flaw", b"lawn") as i32));
    }

    #[test]
    fn gotoh_empty_sequence() {
        let s = dna_scheme(2, 1);
        let al = gotoh(b"", b"ACGT", &s).unwrap();
        // All-gap alignment of 4: -(open + 4*extend) = -(2+4) = -6
        assert_eq!(al.score, -6);
        assert_eq!(al.len(), 4);
        assert!(al.row1.iter().all(|&c| c == b'-'));
    }

    #[test]
    fn nw_over_cap_falls_back_to_hirschberg() {
        // With a tiny injected cap, NW must NOT allocate the O(nm)
        // matrix; it dispatches to linear-space Hirschberg, which under
        // a linear gap gives the *identical* optimal score. We assert
        // the capped path equals both the full NW (computed with a
        // generous cap) and Hirschberg.
        use crate::pairwise::hirschberg::hirschberg;
        let a = b"ACGTACGTACGT";
        let b = b"ACGTCGTACGT";
        let s = dna_scheme(0, 1);
        // 13*13 = 169 cells; a cap of 8 forces the fallback.
        let capped = needleman_wunsch_capped(a, b, &s, 8).unwrap();
        let full = needleman_wunsch_capped(a, b, &s, usize::MAX).unwrap();
        let hb = hirschberg(a, b, &s).unwrap();
        assert_eq!(capped.score, full.score, "fallback score == full NW score");
        assert_eq!(capped.score, hb.score, "fallback delegates to Hirschberg");
        assert_eq!(capped.span1, (0, a.len()));
        assert_eq!(capped.span2, (0, b.len()));
    }

    #[test]
    fn gotoh_over_cap_errors_not_silent_fallback() {
        // Gotoh is affine; Hirschberg models a linear gap, so a fallback
        // would silently change the scoring model. The oversized affine
        // alignment must be rejected with TooLarge instead.
        use crate::error::AlignError;
        let a = b"ACGTACGTACGT";
        let b = b"ACGTCGTACGT";
        let s = dna_scheme(10, 1);
        let err = gotoh_capped(a, b, &s, 8).unwrap_err();
        assert!(
            matches!(err, AlignError::TooLarge { .. }),
            "expected TooLarge, got {err:?}"
        );
        // A generous cap still computes normally.
        assert!(gotoh_capped(a, b, &s, usize::MAX).is_ok());
    }

    #[test]
    fn gotoh_long_leading_gap_does_not_overflow_init() {
        // The init row computes -(open + extend*len); the i64-clamped
        // helper keeps this correct. Exercise a real (small) leading-gap
        // alignment to confirm the score is sane and negative.
        let s = dna_scheme(5, 2);
        let al = gotoh(b"", b"ACGTACGT", &s).unwrap();
        // -(open + 8*extend) = -(5 + 16) = -21
        assert_eq!(al.score, -21);
    }
}
