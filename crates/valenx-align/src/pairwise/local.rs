//! Local alignment — Smith-Waterman.
//!
//! [`smith_waterman`] finds the single highest-scoring *local* segment
//! pair: the DP is the same O(nm) recurrence as Needleman-Wunsch but
//! every cell is floored at `0`, and the traceback starts at the
//! global maximum and stops at the first `0`. The returned
//! [`Alignment`]'s `span1` / `span2` are the half-open coordinate
//! ranges of the aligned segment within each input.
//!
//! Affine gaps are supported through the same three-matrix Gotoh
//! recurrence used by [`crate::pairwise::global::gotoh`], floored at
//! zero.

use super::result::Alignment;
use crate::error::Result;
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};
use crate::matrix::{GapCost, ScoringScheme};

const NEG_INF: i32 = i32::MIN / 4;

/// Smith-Waterman local alignment with an affine gap penalty.
///
/// Returns the best-scoring local segment pair. If no positively
/// scoring segment exists the result is an empty alignment with score
/// `0` and zero-length spans. Runs in O(`a.len()` · `b.len()`) time
/// and space.
///
/// Returns [`AlignError::TooLarge`](crate::error::AlignError::TooLarge)
/// when the `(n+1)·(m+1)` matrices would exceed
/// [`MAX_DP_CELLS`]; local alignment has no
/// linear-space variant here, so an oversized input is rejected.
pub fn smith_waterman(a: &[u8], b: &[u8], scheme: &ScoringScheme) -> Result<Alignment> {
    smith_waterman_capped(a, b, scheme, MAX_DP_CELLS)
}

/// [`smith_waterman`] with an explicit cell cap (test seam).
fn smith_waterman_capped(
    a: &[u8],
    b: &[u8],
    scheme: &ScoringScheme,
    max_cells: usize,
) -> Result<Alignment> {
    let n = a.len();
    let m = b.len();
    let GapCost { open, extend } = scheme.gap;
    let w = m + 1;

    // Four matrices of (n+1)·(m+1) — bound the allocation.
    check_dp_size_with(n + 1, m + 1, max_cells)?;

    let mut mm = vec![0i32; (n + 1) * w];
    let mut ix = vec![NEG_INF; (n + 1) * w];
    let mut iy = vec![NEG_INF; (n + 1) * w];
    // Traceback: 0 stop, 1 diag(from mm), 2 up(ix), 3 left(iy).
    let mut tb = vec![0u8; (n + 1) * w];

    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=n {
        for j in 1..=m {
            let idx = i * w + j;
            let sub = scheme.sub(a[i - 1], b[j - 1]);

            // gap in b (consume a)
            let up = (i - 1) * w + j;
            ix[idx] = (mm[up] - open - extend).max(ix[up] - extend);
            // gap in a (consume b)
            let left = i * w + j - 1;
            iy[idx] = (mm[left] - open - extend).max(iy[left] - extend);

            let diag = mm[(i - 1) * w + j - 1] + sub;
            // Local: floor at 0.
            let mut cell = 0;
            let mut op = 0u8;
            if diag > cell {
                cell = diag;
                op = 1;
            }
            if ix[idx] > cell {
                cell = ix[idx];
                op = 2;
            }
            if iy[idx] > cell {
                cell = iy[idx];
                op = 3;
            }
            mm[idx] = cell;
            tb[idx] = op;

            if cell > best {
                best = cell;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best == 0 {
        // No positively scoring segment.
        return Alignment::new(Vec::new(), Vec::new(), 0, (0, 0), (0, 0));
    }

    // Traceback from the maximum cell to the first zero.
    let mut row1 = Vec::new();
    let mut row2 = Vec::new();
    let (mut i, mut j) = (best_i, best_j);
    while i > 0 && j > 0 {
        let idx = i * w + j;
        match tb[idx] {
            0 => break,
            1 => {
                row1.push(a[i - 1]);
                row2.push(b[j - 1]);
                i -= 1;
                j -= 1;
            }
            2 => {
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
    let span1 = (i, best_i);
    let span2 = (j, best_j);
    Alignment::new(row1, row2, best, span1, span2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn dna_scheme(open: i32, extend: i32) -> ScoringScheme {
        ScoringScheme::new(
            SubstitutionMatrix::dna_simple(2, -1),
            GapCost::new(open, extend),
        )
    }

    #[test]
    fn sw_finds_embedded_segment() {
        // The shared core "GATTACA" is embedded in flanking noise.
        let a = b"TTTTGATTACATTTT";
        let b = b"CCCGATTACACCC";
        let s = dna_scheme(4, 1);
        let al = smith_waterman(a, b, &s).unwrap();
        assert_eq!(al.row1_str(), "GATTACA");
        assert_eq!(al.row2_str(), "GATTACA");
        assert_eq!(al.score, 14); // 7 matches * 2
    }

    #[test]
    fn sw_spans_point_at_segment() {
        let a = b"TTTTGATTACATTTT";
        let b = b"CCCGATTACACCC";
        let s = dna_scheme(4, 1);
        let al = smith_waterman(a, b, &s).unwrap();
        // GATTACA starts at index 4 of a, index 3 of b.
        assert_eq!(al.span1, (4, 11));
        assert_eq!(al.span2, (3, 10));
    }

    #[test]
    fn sw_no_similarity_is_empty() {
        let s = dna_scheme(4, 1);
        let al = smith_waterman(b"AAAAAA", b"TTTTTT", &s).unwrap();
        assert_eq!(al.score, 0);
        assert!(al.is_empty());
    }

    #[test]
    fn sw_local_beats_global_with_noise() {
        // A perfect core flanked by mismatches: local score should be
        // exactly the core score, ignoring the flanks.
        let a = b"XXXACGTACGTXXX";
        let b = b"ACGTACGT";
        let s = ScoringScheme::new(
            SubstitutionMatrix::dna_simple(2, -1),
            GapCost::new(4, 1),
        );
        let al = smith_waterman(a, b, &s).unwrap();
        assert_eq!(al.score, 16); // 8 matches * 2, flanks excluded
    }

    #[test]
    fn sw_protein_local() {
        let s = ScoringScheme::new(SubstitutionMatrix::blosum62(), GapCost::new(11, 1));
        let al = smith_waterman(b"AAMKVLAA", b"CCMKVLCC", &s).unwrap();
        assert_eq!(al.row1_str(), "MKVL");
        assert!(al.score > 0);
    }

    #[test]
    fn sw_over_cap_errors() {
        use crate::error::AlignError;
        let a = b"TTTTGATTACATTTT";
        let b = b"CCCGATTACACCC";
        let s = dna_scheme(4, 1);
        // A tiny cap rejects without allocating the four matrices.
        let err = smith_waterman_capped(a, b, &s, 8).unwrap_err();
        assert!(matches!(err, AlignError::TooLarge { .. }), "got {err:?}");
        // A generous cap computes normally — same as the public fn.
        let ok = smith_waterman_capped(a, b, &s, usize::MAX).unwrap();
        assert_eq!(ok.score, 14);
    }
}
