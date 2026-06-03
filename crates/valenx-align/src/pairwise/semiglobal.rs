//! Semi-global ("overlap" / "glocal") alignment with free end gaps.
//!
//! Semi-global alignment places *one* sequence entirely but lets the
//! *other*'s leading and trailing gaps be free — the model for read
//! mapping (a short read against a long reference) and for detecting
//! the overlap between two assembly contigs.
//!
//! [`semi_global`] takes an [`EndGapPolicy`] that independently frees
//! the start / end gaps of each sequence; the common presets are
//! exposed as [`EndGapPolicy::overlap`] (free all four ends — true
//! overlap detection) and [`EndGapPolicy::fit`] (free both ends of the
//! reference only — fit the whole query into the reference).

use super::result::Alignment;
use crate::error::Result;
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};
use crate::matrix::ScoringScheme;

/// A very negative floor for the free-end-gap init rows — clamps the
/// leading-gap cost in `i64` space so a long sequence cannot overflow
/// `i32`, while leaving headroom for the inner-loop additions.
const NEG_FLOOR: i32 = i32::MIN / 4;

/// Which end gaps are free (unpenalised) in a semi-global alignment.
///
/// "Sequence 1" is `a`, "sequence 2" is `b` in the
/// [`semi_global`] call. A freed end gap costs nothing, so the
/// alignment can slide that end without penalty.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct EndGapPolicy {
    /// Free leading gaps of sequence 1 (`a`).
    pub free_start1: bool,
    /// Free trailing gaps of sequence 1 (`a`).
    pub free_end1: bool,
    /// Free leading gaps of sequence 2 (`b`).
    pub free_start2: bool,
    /// Free trailing gaps of sequence 2 (`b`).
    pub free_end2: bool,
}

impl EndGapPolicy {
    /// Free every end gap — classic overlap / dovetail detection.
    pub fn overlap() -> Self {
        EndGapPolicy {
            free_start1: true,
            free_end1: true,
            free_start2: true,
            free_end2: true,
        }
    }

    /// Fit the whole of sequence 2 (`b`, the read) somewhere inside
    /// sequence 1 (`a`, the reference): the read is placed in full while
    /// the reference may overhang at both ends for free.
    ///
    /// The reference overhang appears in the alignment as reference
    /// residues opposite a *gap in the read* — i.e. as end gaps of
    /// sequence 2 — so it is `free_start2` / `free_end2` that must be
    /// set, while the read's own ends stay penalised so it cannot
    /// dangle off the reference.
    pub fn fit() -> Self {
        EndGapPolicy {
            free_start1: false,
            free_end1: false,
            free_start2: true,
            free_end2: true,
        }
    }

    /// Penalise every end gap — equivalent to a true global alignment.
    pub fn global() -> Self {
        EndGapPolicy {
            free_start1: false,
            free_end1: false,
            free_start2: false,
            free_end2: false,
        }
    }
}

/// Semi-global alignment with a **linear** gap penalty and a
/// configurable end-gap policy.
///
/// Uses `scheme.gap.extend` as the per-residue gap penalty (this is the
/// linear-gap routine; the `open` field is ignored). The returned
/// `span1` / `span2` reflect which portion of each sequence is
/// non-end-gap aligned. Runs in O(`a.len()` · `b.len()`).
pub fn semi_global(
    a: &[u8],
    b: &[u8],
    scheme: &ScoringScheme,
    policy: EndGapPolicy,
) -> Result<Alignment> {
    semi_global_capped(a, b, scheme, policy, MAX_DP_CELLS)
}

/// [`semi_global`] with an explicit cell cap (test seam). Over the cap
/// it returns [`AlignError::TooLarge`](crate::error::AlignError::TooLarge);
/// semi-global has no linear-space variant here.
fn semi_global_capped(
    a: &[u8],
    b: &[u8],
    scheme: &ScoringScheme,
    policy: EndGapPolicy,
    max_cells: usize,
) -> Result<Alignment> {
    let n = a.len();
    let m = b.len();
    let g = scheme.gap.extend;
    let w = m + 1;

    // Two matrices of (n+1)·(m+1) — bound the allocation.
    check_dp_size_with(n + 1, m + 1, max_cells)?;

    let mut s = vec![0i32; (n + 1) * w];
    let mut tb = vec![0u8; (n + 1) * w]; // 0 diag, 1 up, 2 left

    // First column: `b` has not started, `a`'s leading residues are
    // consumed as gaps IN `b` — so this is a leading gap of sequence 2,
    // freed by `free_start2` (not `free_start1`). The penalised cost is
    // computed in i64 and clamped so a long sequence cannot overflow.
    for i in 1..=n {
        s[i * w] = if policy.free_start2 {
            0
        } else {
            super::affine_init(0, g, i, NEG_FLOOR)
        };
        tb[i * w] = 1;
    }
    // First row: `a` has not started, `b`'s leading residues are
    // consumed as gaps IN `a` — a leading gap of sequence 1, freed by
    // `free_start1`.
    for j in 1..=m {
        s[j] = if policy.free_start1 {
            0
        } else {
            super::affine_init(0, g, j, NEG_FLOOR)
        };
        tb[j] = 2;
    }

    for i in 1..=n {
        for j in 1..=m {
            let diag = s[(i - 1) * w + j - 1] + scheme.sub(a[i - 1], b[j - 1]);
            // Up = gap in b. Free only on the last column when
            // free_end2 (b is exhausted, trailing gap of b).
            let up_free = j == m && policy.free_end2;
            let up = s[(i - 1) * w + j] - if up_free { 0 } else { g };
            // Left = gap in a. Free only on the last row when
            // free_end1.
            let left_free = i == n && policy.free_end1;
            let left = s[i * w + j - 1] - if left_free { 0 } else { g };

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
            s[i * w + j] = best;
            tb[i * w + j] = op;
        }
    }

    // Find the traceback start: the best cell allowed to be a free end.
    let (mut i, mut j, score) = find_end(&s, w, n, m, policy);

    // Walk the free trailing gaps (cells from (i,j) to (n,m)) without
    // emitting penalised columns — but we still emit them as end gaps.
    let mut tail1 = Vec::new();
    let mut tail2 = Vec::new();
    {
        let (mut ti, mut tj) = (n, m);
        while ti > i {
            tail1.push(a[ti - 1]);
            tail2.push(b'-');
            ti -= 1;
        }
        while tj > j {
            tail1.push(b'-');
            tail2.push(b[tj - 1]);
            tj -= 1;
        }
    }

    // Core traceback.
    let mut row1 = Vec::new();
    let mut row2 = Vec::new();
    let (end_i, end_j) = (i, j);
    while i > 0 && j > 0 {
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
    let core_start_i = i;
    let core_start_j = j;
    // Leading end gaps.
    while i > 0 {
        row1.push(a[i - 1]);
        row2.push(b'-');
        i -= 1;
    }
    while j > 0 {
        row1.push(b'-');
        row2.push(b[j - 1]);
        j -= 1;
    }

    row1.reverse();
    row2.reverse();
    // Append the trailing tail (tail vectors are in reverse order).
    tail1.reverse();
    tail2.reverse();
    row1.extend(tail1);
    row2.extend(tail2);

    // span = the non-leading-end-gap, non-trailing-end-gap core region.
    let span1 = (core_start_i, end_i.max(core_start_i));
    let span2 = (core_start_j, end_j.max(core_start_j));
    Alignment::new(row1, row2, score, span1, span2)
}

/// Finds the cell to start traceback from: the best score over the
/// cells reachable by free trailing gaps. Returns `(i, j, score)`.
///
/// When a whole free-end edge ties at the maximum score (e.g. the entire
/// last column is equal because the trailing gaps of `b` are free), the
/// cell *closest to the body of the alignment* is chosen — the smallest
/// `i` on the column scan, the smallest `j` on the row scan — so the
/// tied stretch becomes a genuine free end gap and the reported `span`
/// excludes it. (A strict `>` starting from `s[n][m]` would instead
/// pin the start at the far corner and swallow the free region into the
/// span.)
fn find_end(
    s: &[i32],
    w: usize,
    n: usize,
    m: usize,
    policy: EndGapPolicy,
) -> (usize, usize, i32) {
    let mut best = s[n * w + m];
    let mut bi = n;
    let mut bj = m;
    // Free trailing gap of b => last column free: scan column m, taking
    // the smallest row index that achieves the maximum.
    if policy.free_end2 {
        for i in 0..=n {
            let v = s[i * w + m];
            if v > best || (v == best && i < bi && bj == m) {
                best = v;
                bi = i;
                bj = m;
            }
        }
    }
    // Free trailing gap of a => last row free: scan row n, taking the
    // smallest column index that achieves the maximum.
    if policy.free_end1 {
        for j in 0..=m {
            let v = s[n * w + j];
            if v > best || (v == best && j < bj && bi == n) {
                best = v;
                bi = n;
                bj = j;
            }
        }
    }
    (bi, bj, best)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn dna(open: i32, ext: i32) -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(1, -1), GapCost::new(open, ext))
    }

    #[test]
    fn overlap_dovetail_no_end_penalty() {
        // a ends with the same suffix that b starts with.
        let a = b"AAAAGATTACA";
        let b = b"GATTACATTTT";
        let s = dna(0, 2);
        let al = semi_global(a, b, &s, EndGapPolicy::overlap()).unwrap();
        // The 7-residue overlap GATTACA matches: score 7, end gaps free.
        assert_eq!(al.score, 7);
    }

    #[test]
    fn fit_short_query_into_reference() {
        // b fits entirely inside a; a's ends are free.
        let a = b"TTTTTACGTACGTTTTTT";
        let b = b"ACGTACGT";
        let s = dna(0, 2);
        let al = semi_global(a, b, &s, EndGapPolicy::fit()).unwrap();
        assert_eq!(al.score, 8); // query placed perfectly
        // The query span inside the reference.
        assert_eq!(al.span1, (5, 13));
    }

    #[test]
    fn global_policy_penalises_ends() {
        // With the global policy this is just NW; long end gaps hurt.
        let a = b"ACGT";
        let b = b"ACGTACGTACGT";
        let s = dna(0, 1);
        let glob = semi_global(a, b, &s, EndGapPolicy::global()).unwrap();
        let over = semi_global(a, b, &s, EndGapPolicy::overlap()).unwrap();
        assert!(over.score > glob.score, "free ends must score higher");
    }

    #[test]
    fn rows_stay_equal_length() {
        let s = dna(0, 1);
        let al = semi_global(b"AAGATTACA", b"GATTACAGG", &s, EndGapPolicy::overlap()).unwrap();
        assert_eq!(al.row1.len(), al.row2.len());
    }

    #[test]
    fn semiglobal_over_cap_errors() {
        use crate::error::AlignError;
        let s = dna(0, 1);
        let err =
            semi_global_capped(b"ACGT", b"ACGTACGTACGT", &s, EndGapPolicy::overlap(), 8)
                .unwrap_err();
        assert!(matches!(err, AlignError::TooLarge { .. }), "got {err:?}");
        assert!(
            semi_global_capped(b"ACGT", b"ACGTACGTACGT", &s, EndGapPolicy::overlap(), usize::MAX)
                .is_ok()
        );
    }

    #[test]
    fn semiglobal_global_policy_long_leading_gap_no_overflow() {
        // Global policy penalises the leading gap; the i64-clamped init
        // keeps it correct. Small case: query ACGT vs 12-mer, global.
        let s = dna(0, 1);
        let glob = semi_global(b"ACGT", b"ACGTACGTACGT", &s, EndGapPolicy::global()).unwrap();
        // Score must be a sane negative-or-small number, never a wrapped
        // positive overflow artefact.
        assert!(glob.score <= 4, "global score {} sane", glob.score);
    }
}
