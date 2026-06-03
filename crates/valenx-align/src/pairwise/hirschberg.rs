//! Hirschberg linear-space global alignment.
//!
//! Needleman-Wunsch needs O(nm) memory for its traceback matrix —
//! prohibitive for megabase sequences. Hirschberg's algorithm computes
//! the *same optimal global alignment* in O(min(n, m)) space by
//! divide-and-conquer: it splits sequence `a` in half, uses the
//! linear-space "score-only" DP forward over the first half and
//! backward over the second half to find the column where the optimal
//! path crosses the split row, and recurses on the two sub-problems.
//!
//! [`hirschberg`] uses a **linear** gap penalty (`scheme.gap.extend`)
//! and returns a result identical in score to
//! [`crate::pairwise::global::needleman_wunsch`].

use super::result::Alignment;
use crate::error::Result;
use crate::matrix::ScoringScheme;

/// Hirschberg global alignment in O(min(n, m)) space, linear gap cost.
///
/// The result is the optimal global alignment — the same score
/// `needleman_wunsch` would produce — computed without an O(nm)
/// traceback matrix. Time stays O(nm).
pub fn hirschberg(a: &[u8], b: &[u8], scheme: &ScoringScheme) -> Result<Alignment> {
    let g = scheme.gap.extend;
    let (row1, row2) = solve(a, b, scheme, g);
    let score = score_rows(&row1, &row2, scheme, g);
    Alignment::new(row1, row2, score, (0, a.len()), (0, b.len()))
}

/// Recursive divide-and-conquer. Returns the two gapped rows.
fn solve(a: &[u8], b: &[u8], scheme: &ScoringScheme, g: i32) -> (Vec<u8>, Vec<u8>) {
    let n = a.len();
    let m = b.len();

    // Base cases.
    if n == 0 {
        return (vec![b'-'; m], b.to_vec());
    }
    if m == 0 {
        return (a.to_vec(), vec![b'-'; n]);
    }
    if n == 1 {
        return align_single_row(a[0], b, scheme, g);
    }
    if m == 1 {
        // Symmetric: align b[0] against a, then swap rows.
        let (r_b, r_a) = align_single_row(b[0], a, scheme, g);
        return (r_a, r_b);
    }

    // Split a at its midpoint.
    let mid = n / 2;
    let score_l = nw_score_row(&a[..mid], b, scheme, g);
    let score_r_rev = nw_score_row_rev(&a[mid..], b, scheme, g);

    // Find the column j minimising score_l[j] + score_r_rev[m - j].
    let mut best = i32::MIN;
    let mut split = 0;
    for j in 0..=m {
        let v = score_l[j].saturating_add(score_r_rev[m - j]);
        if v > best {
            best = v;
            split = j;
        }
    }

    let (l1, l2) = solve(&a[..mid], &b[..split], scheme, g);
    let (r1, r2) = solve(&a[mid..], &b[split..], scheme, g);
    let mut row1 = l1;
    row1.extend(r1);
    let mut row2 = l2;
    row2.extend(r2);
    (row1, row2)
}

/// Aligns a single residue `x` against the whole of `rest`, choosing
/// the cheapest column to place the match. Returns `(row_x, row_rest)`.
fn align_single_row(x: u8, rest: &[u8], scheme: &ScoringScheme, g: i32) -> (Vec<u8>, Vec<u8>) {
    let m = rest.len();
    // Option A: x matches some rest[j], the others are gaps in row_x.
    let mut best_score = i32::MIN;
    let mut best_j: Option<usize> = None;
    for (j, &rj) in rest.iter().enumerate() {
        // (j gaps before) + match + (m-1-j gaps after); each gap costs g.
        let s = scheme.sub(x, rj) - g * (m as i32 - 1);
        if s > best_score {
            best_score = s;
            best_j = Some(j);
        }
    }
    // Option B: x is a gap-insert, everything in rest gapped on row_x
    // and x placed as a column with a gap in row_rest.
    let all_gap = -g * (m as i32 + 1);
    if all_gap > best_score {
        best_j = None;
    }

    match best_j {
        Some(j) => {
            let mut row_x = vec![b'-'; m];
            row_x[j] = x;
            (row_x, rest.to_vec())
        }
        None => {
            // x in its own column, then all of rest gapped against x.
            let mut row_x = vec![x];
            row_x.extend(std::iter::repeat_n(b'-', m));
            let mut row_rest = vec![b'-'];
            row_rest.extend_from_slice(rest);
            (row_x, row_rest)
        }
    }
}

/// Score-only forward NW: returns the last DP row (length `m + 1`).
fn nw_score_row(a: &[u8], b: &[u8], scheme: &ScoringScheme, g: i32) -> Vec<i32> {
    let m = b.len();
    let mut prev: Vec<i32> = (0..=m).map(|j| -(g * j as i32)).collect();
    let mut cur = vec![0i32; m + 1];
    for &ca in a {
        cur[0] = prev[0] - g;
        for j in 1..=m {
            let diag = prev[j - 1] + scheme.sub(ca, b[j - 1]);
            let up = prev[j] - g;
            let left = cur[j - 1] - g;
            cur[j] = diag.max(up).max(left);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev
}

/// Score-only NW over the *reversed* sequences — the last row of the
/// DP run right-to-left. Used for the backward half of the split.
fn nw_score_row_rev(a: &[u8], b: &[u8], scheme: &ScoringScheme, g: i32) -> Vec<i32> {
    let ar: Vec<u8> = a.iter().rev().copied().collect();
    let br: Vec<u8> = b.iter().rev().copied().collect();
    nw_score_row(&ar, &br, scheme, g)
}

/// Scores a pair of already-gapped rows under the linear-gap scheme.
fn score_rows(row1: &[u8], row2: &[u8], scheme: &ScoringScheme, g: i32) -> i32 {
    let mut total = 0;
    for (&a, &b) in row1.iter().zip(row2) {
        if a == b'-' || b == b'-' {
            total -= g;
        } else {
            total += scheme.sub(a, b);
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};
    use crate::pairwise::global::needleman_wunsch;

    fn dna(g: i32) -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(1, -1), GapCost::new(0, g))
    }

    #[test]
    fn hirschberg_matches_nw_score() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGTACGT", b"ACGTACGT"),
            (b"ACGTACGT", b"ACGACGT"),
            (b"AAAAAAAA", b"AAAA"),
            (b"GATTACA", b"GCATGCU".as_slice()),
            (b"ACGTACGTACGTACGT", b"ACGTTCGTACATACGT"),
        ];
        for &(a, b) in cases {
            let s = dna(1);
            let nw = needleman_wunsch(a, b, &s).unwrap();
            let hb = hirschberg(a, b, &s).unwrap();
            assert_eq!(
                hb.score,
                nw.score,
                "score mismatch on {:?}/{:?}",
                std::str::from_utf8(a),
                std::str::from_utf8(b),
            );
        }
    }

    #[test]
    fn hirschberg_rows_are_consistent() {
        let s = dna(1);
        let hb = hirschberg(b"ACGTACGTAC", b"ACGTCGTAC", &s).unwrap();
        assert_eq!(hb.row1.len(), hb.row2.len());
        // Stripping gaps must recover the originals.
        let r1: Vec<u8> = hb.row1.iter().copied().filter(|&c| c != b'-').collect();
        let r2: Vec<u8> = hb.row2.iter().copied().filter(|&c| c != b'-').collect();
        assert_eq!(r1, b"ACGTACGTAC");
        assert_eq!(r2, b"ACGTCGTAC");
    }

    #[test]
    fn hirschberg_empty_inputs() {
        let s = dna(2);
        let hb = hirschberg(b"", b"ACGT", &s).unwrap();
        assert_eq!(hb.score, -8);
        let hb2 = hirschberg(b"ACG", b"", &s).unwrap();
        assert_eq!(hb2.score, -6);
    }

    #[test]
    fn hirschberg_protein() {
        let s = ScoringScheme::new(SubstitutionMatrix::blosum62(), GapCost::new(0, 4));
        let nw = needleman_wunsch(b"MKVLAAGGK", b"MKVLAGGK", &s).unwrap();
        let hb = hirschberg(b"MKVLAAGGK", b"MKVLAGGK", &s).unwrap();
        assert_eq!(hb.score, nw.score);
    }
}
