//! Banded global alignment — diagonal-restricted dynamic programming.
//!
//! When two sequences are known to be similar (so the optimal path
//! stays close to the main diagonal), full O(nm) DP wastes work.
//! [`banded_global`] restricts the DP to a band of `±k` diagonals
//! around the main diagonal, giving O(n·k) time. This is the inner
//! loop of seed extension and the standard speed-up for read mapping.
//!
//! The band is *anti-diagonal-indexed*: cell `(i, j)` lies on diagonal
//! `d = j - i`, and only `|d| <= k` (plus the length-difference offset)
//! cells are computed. If the band is too narrow to contain a valid
//! path to the far corner the routine returns
//! [`AlignError::Dimension`] — the caller should widen `k` and retry.
//!
//! [`banded_affine`] is the Gotoh three-matrix version, banded — the
//! BWA-MEM / minimap2 inner-loop primitive for chained-anchor
//! base-level extension. It uses the same band as [`banded_global`]
//! but with `mm` / `ix` / `iy` layers so affine open/extend costs are
//! honoured exactly.
//!
//! # Optimality is conditional on the band
//!
//! Banded DP is **not** an exact global aligner. The returned alignment
//! is the optimal one *only if the true optimal path stays within the
//! `±k` diagonal band*. If the optimal alignment would require more
//! than `k` net indels at some prefix (i.e. it strays past `k`
//! diagonals from the main diagonal), these routines return the best
//! path that *does* stay in the band — a possibly globally-suboptimal
//! alignment — rather than an error. An error
//! ([`AlignError::Dimension`]) is returned only in the degenerate case
//! where *no* in-band path can reach the far corner at all (e.g. the
//! length difference already exceeds `k`). To guarantee a globally
//! optimal result, choose `k` ≥ the maximum number of indels the
//! optimal alignment could contain, or fall back to full
//! [`crate::pairwise::global`] DP.

use super::result::Alignment;
use crate::error::{AlignError, Result};
use crate::limits::check_dp_size;
use crate::matrix::ScoringScheme;

const NEG_INF: i32 = i32::MIN / 4;

/// Banded global alignment with a **linear** gap penalty.
///
/// `k` is the half-bandwidth: a path may stray at most `k` diagonals
/// from the line joining the two corners. Choose `k` ≥ the expected
/// number of indels. Uses `scheme.gap.extend` per gap residue.
///
/// Returns [`AlignError::Invalid`] for `k == 0` with sequences of
/// different lengths (no valid path can exist), and
/// [`AlignError::Dimension`] if the band as specified cannot reach the
/// (`n`, `m`) corner.
///
/// **Optimality is conditional on `k`.** The returned alignment is
/// globally optimal *only if the true optimal path stays within the
/// `±k` band*. If the optimal alignment needs more than `k` net indels
/// it falls outside the band, and this function returns the best
/// *in-band* (possibly globally-suboptimal) alignment — **not** an
/// error. Widen `k` to recover the exact optimum.
pub fn banded_global(a: &[u8], b: &[u8], scheme: &ScoringScheme, k: usize) -> Result<Alignment> {
    let n = a.len();
    let m = b.len();
    let g = scheme.gap.extend;

    let len_diff = (n as isize - m as isize).unsigned_abs();
    if k < len_diff {
        return Err(AlignError::dimension(format!(
            "band half-width {k} too small for length difference {len_diff}"
        )));
    }

    // Diagonal d = j - i ranges over [-k_lo, k_hi]. To always reach the
    // (n,m) corner (diagonal m-n) the band must include it.
    let d_lo = -(k as isize);
    let d_hi = k as isize;
    let band_w = (d_hi - d_lo + 1) as usize;

    // Bound the banded allocation: a caller passing a large `k` would
    // otherwise request `(n+1)·band_w` cells unchecked. Same chokepoint
    // the full-matrix DP sites use.
    check_dp_size(n + 1, band_w)?;

    // s[i][col] where col = (j - i) - d_lo. Out-of-band cells are NEG_INF.
    let mut s = vec![NEG_INF; (n + 1) * band_w];
    let mut tb = vec![0u8; (n + 1) * band_w]; // 0 diag, 1 up, 2 left

    let col_of = |i: usize, j: usize| -> Option<usize> {
        let d = j as isize - i as isize;
        if d < d_lo || d > d_hi {
            None
        } else {
            Some((d - d_lo) as usize)
        }
    };

    // Origin.
    if let Some(c) = col_of(0, 0) {
        s[c] = 0;
    }
    // First row (i=0): gaps in a.
    for j in 1..=m {
        if let Some(c) = col_of(0, j) {
            s[c] = -(g * j as i32);
            tb[c] = 2;
        }
    }
    // First column (j=0): gaps in b.
    for i in 1..=n {
        if let Some(c) = col_of(i, 0) {
            s[i * band_w + c] = -(g * i as i32);
            tb[i * band_w + c] = 1;
        }
    }

    for i in 1..=n {
        // j range that stays in band for this row.
        let j_min = ((i as isize + d_lo).max(1)) as usize;
        let j_max = ((i as isize + d_hi).min(m as isize)).max(0) as usize;
        for j in j_min..=j_max {
            let c = match col_of(i, j) {
                Some(c) => c,
                None => continue,
            };
            // diag (i-1, j-1)
            let diag = match col_of(i - 1, j - 1) {
                Some(pc) => {
                    let v = s[(i - 1) * band_w + pc];
                    if v <= NEG_INF {
                        NEG_INF
                    } else {
                        v + scheme.sub(a[i - 1], b[j - 1])
                    }
                }
                None => NEG_INF,
            };
            // up (i-1, j)
            let up = match col_of(i - 1, j) {
                Some(pc) => {
                    let v = s[(i - 1) * band_w + pc];
                    if v <= NEG_INF {
                        NEG_INF
                    } else {
                        v - g
                    }
                }
                None => NEG_INF,
            };
            // left (i, j-1)
            let left = if j >= 1 {
                match col_of(i, j - 1) {
                    Some(pc) => {
                        let v = s[i * band_w + pc];
                        if v <= NEG_INF {
                            NEG_INF
                        } else {
                            v - g
                        }
                    }
                    None => NEG_INF,
                }
            } else {
                NEG_INF
            };

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
            s[i * band_w + c] = best;
            tb[i * band_w + c] = op;
        }
    }

    let final_col = col_of(n, m)
        .ok_or_else(|| AlignError::dimension("band does not reach the (n,m) corner"))?;
    let score = s[n * band_w + final_col];
    if score <= NEG_INF {
        return Err(AlignError::dimension(
            "no in-band alignment path exists; widen the band",
        ));
    }

    // Traceback.
    let mut row1 = Vec::new();
    let mut row2 = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let c = col_of(i, j).ok_or_else(|| AlignError::dimension("traceback left the band"))?;
        match tb[i * band_w + c] {
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
    Alignment::new(row1, row2, score, (0, n), (0, m))
}

/// Banded **affine-gap** global alignment — the Gotoh three-matrix DP
/// restricted to a band of `±k` diagonals.
///
/// Time is `O((n + m) · k)`, space `O((n + m) · k)`. The result is a
/// proper [`Alignment`] (gapped rows + score + spans) — drop-in for the
/// caller wherever it would have run [`crate::pairwise::global::gotoh`]
/// on a slice it already knows is near-diagonal (e.g. the chained
/// region of a read mapper).
///
/// **Optimality is conditional on `k`.** As with [`banded_global`], the
/// returned alignment is globally optimal *only if the optimal path
/// stays within the `±k` band*. If the optimum needs more than `k` net
/// indels, the best *in-band* (possibly globally-suboptimal) alignment
/// is returned rather than an error; an error is raised only when no
/// in-band path reaches the corner. Widen `k` for an exact result.
pub fn banded_affine(
    a: &[u8],
    b: &[u8],
    scheme: &ScoringScheme,
    k: usize,
) -> Result<Alignment> {
    let n = a.len();
    let m = b.len();
    let open = scheme.gap.open;
    let extend = scheme.gap.extend;

    let len_diff = (n as isize - m as isize).unsigned_abs();
    if k < len_diff {
        return Err(AlignError::dimension(format!(
            "band half-width {k} too small for length difference {len_diff}"
        )));
    }

    let d_lo = -(k as isize);
    let d_hi = k as isize;
    let band_w = (d_hi - d_lo + 1) as usize;

    // Bound the banded allocation before laying out the six DP layers; a
    // large `k` would otherwise request `(n+1)·band_w` cells unchecked.
    check_dp_size(n + 1, band_w)?;

    let col_of = |i: usize, j: usize| -> Option<usize> {
        let d = j as isize - i as isize;
        if d < d_lo || d > d_hi {
            None
        } else {
            Some((d - d_lo) as usize)
        }
    };

    // Three DP layers — mm, ix (gap-in-b ending), iy (gap-in-a ending).
    let mut mm = vec![NEG_INF; (n + 1) * band_w];
    let mut ix = vec![NEG_INF; (n + 1) * band_w];
    let mut iy = vec![NEG_INF; (n + 1) * band_w];
    // Traceback: which layer the cell came from. 0=mm, 1=ix, 2=iy.
    let mut tm = vec![0u8; (n + 1) * band_w];
    let mut tx = vec![0u8; (n + 1) * band_w];
    let mut ty = vec![0u8; (n + 1) * band_w];

    if let Some(c) = col_of(0, 0) {
        mm[c] = 0;
    }
    // First row (i = 0): only iy gets a value (gaps in a).
    for j in 1..=m {
        if let Some(c) = col_of(0, j) {
            iy[c] = -(open + extend * j as i32);
            ty[c] = 2;
        }
    }
    // First column (j = 0): only ix.
    for i in 1..=n {
        if let Some(c) = col_of(i, 0) {
            ix[i * band_w + c] = -(open + extend * i as i32);
            tx[i * band_w + c] = 1;
        }
    }

    for i in 1..=n {
        let j_min = ((i as isize + d_lo).max(1)) as usize;
        let j_max = ((i as isize + d_hi).min(m as isize)).max(0) as usize;
        for j in j_min..=j_max {
            let c = match col_of(i, j) {
                Some(c) => c,
                None => continue,
            };
            let idx = i * band_w + c;
            let sub = scheme.sub(a[i - 1], b[j - 1]);

            // mm: from diag.
            if let Some(pc) = col_of(i - 1, j - 1) {
                let pidx = (i - 1) * band_w + pc;
                let cm = mm[pidx];
                let cx = ix[pidx];
                let cy = iy[pidx];
                let (best, src) = max3(cm, cx, cy);
                if best > NEG_INF / 2 {
                    mm[idx] = best + sub;
                    tm[idx] = src;
                }
            }
            // ix: gap in b — consume a, came from (i-1, j).
            if let Some(pc) = col_of(i - 1, j) {
                let pidx = (i - 1) * band_w + pc;
                let open_x = mm[pidx] - open - extend;
                let ext_x = ix[pidx] - extend;
                if open_x >= ext_x {
                    if open_x > NEG_INF / 2 {
                        ix[idx] = open_x;
                        tx[idx] = 0;
                    }
                } else if ext_x > NEG_INF / 2 {
                    ix[idx] = ext_x;
                    tx[idx] = 1;
                }
            }
            // iy: gap in a — consume b, came from (i, j-1).
            if j >= 1 {
                if let Some(pc) = col_of(i, j - 1) {
                    let pidx = i * band_w + pc;
                    let open_y = mm[pidx] - open - extend;
                    let ext_y = iy[pidx] - extend;
                    if open_y >= ext_y {
                        if open_y > NEG_INF / 2 {
                            iy[idx] = open_y;
                            ty[idx] = 0;
                        }
                    } else if ext_y > NEG_INF / 2 {
                        iy[idx] = ext_y;
                        ty[idx] = 2;
                    }
                }
            }
        }
    }

    let final_col = col_of(n, m)
        .ok_or_else(|| AlignError::dimension("band does not reach the (n,m) corner"))?;
    let cm = mm[n * band_w + final_col];
    let cx = ix[n * band_w + final_col];
    let cy = iy[n * band_w + final_col];
    let (score, start_layer) = max3(cm, cx, cy);
    if score <= NEG_INF / 2 {
        return Err(AlignError::dimension(
            "no in-band affine alignment exists; widen the band",
        ));
    }

    // Traceback.
    let mut row1 = Vec::new();
    let mut row2 = Vec::new();
    let (mut i, mut j) = (n, m);
    let mut layer = start_layer;
    while i > 0 || j > 0 {
        let c = col_of(i, j).ok_or_else(|| AlignError::dimension("traceback left the band"))?;
        let idx = i * band_w + c;
        match layer {
            0 => {
                row1.push(a[i - 1]);
                row2.push(b[j - 1]);
                let next = tm[idx];
                i -= 1;
                j -= 1;
                layer = next;
            }
            1 => {
                row1.push(a[i - 1]);
                row2.push(b'-');
                let next = tx[idx];
                i -= 1;
                layer = next;
            }
            _ => {
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

/// Returns `(max, index)` where index is 0/1/2 for the first/second/
/// third argument. Stable: ties resolve toward the first argument.
fn max3(a: i32, b: i32, c: i32) -> (i32, u8) {
    let mut best = a;
    let mut idx = 0u8;
    if b > best {
        best = b;
        idx = 1;
    }
    if c > best {
        best = c;
        idx = 2;
    }
    (best, idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};
    use crate::pairwise::global::{gotoh, needleman_wunsch};

    fn dna(g: i32) -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(1, -1), GapCost::new(0, g))
    }

    #[test]
    fn banded_matches_full_dp_when_band_wide() {
        // A wide enough band must give exactly the NW score.
        let a = b"ACGTACGTACGTAC";
        let b = b"ACGTACTTACGTAC";
        let s = dna(1);
        let full = needleman_wunsch(a, b, &s).unwrap();
        let band = banded_global(a, b, &s, 5).unwrap();
        assert_eq!(band.score, full.score);
        assert_eq!(band.row1.len(), band.row2.len());
    }

    #[test]
    fn banded_identical_zero_band() {
        // Identical sequences: even a zero-width band suffices.
        let a = b"ACGTACGT";
        let s = dna(1);
        let band = banded_global(a, a, &s, 0).unwrap();
        assert_eq!(band.score, 8);
    }

    #[test]
    fn banded_rejects_too_narrow_band() {
        // Length difference 4 cannot fit in a band of half-width 1.
        let a = b"ACGTACGTACGT";
        let b = b"ACGTACGT";
        let s = dna(1);
        assert!(banded_global(a, b, &s, 1).is_err());
        // Band wide enough for the 4-residue length gap works.
        assert!(banded_global(a, b, &s, 4).is_ok());
    }

    #[test]
    fn banded_handles_indel_inside_band() {
        // One deletion; a band of half-width 2 contains the path.
        let a = b"ACGTACGT";
        let b = b"ACGACGT";
        let s = dna(1);
        let full = needleman_wunsch(a, b, &s).unwrap();
        let band = banded_global(a, b, &s, 2).unwrap();
        assert_eq!(band.score, full.score);
    }

    fn affine_scheme(open: i32, extend: i32) -> ScoringScheme {
        ScoringScheme::new(
            SubstitutionMatrix::dna_simple(2, -3),
            GapCost::new(open, extend),
        )
    }

    #[test]
    fn banded_affine_matches_gotoh_when_band_wide() {
        let a = b"ACGTACGTACGTAC";
        let b = b"ACGTACTTACGTAC";
        let s = affine_scheme(5, 2);
        let full = gotoh(a, b, &s).unwrap();
        let band = banded_affine(a, b, &s, 6).unwrap();
        assert_eq!(band.score, full.score);
        assert_eq!(band.row1.len(), band.row2.len());
    }

    #[test]
    fn banded_affine_handles_affine_gap() {
        // Reference has a 4-base run not in the query; affine should
        // keep it as one contiguous gap.
        let a = b"AAAGGGGAAA";
        let b = b"AAAAAA";
        let s = affine_scheme(10, 1);
        let full = gotoh(a, b, &s).unwrap();
        let band = banded_affine(a, b, &s, 4).unwrap();
        assert_eq!(band.score, full.score);
        // Exactly one gap run in row2 (proper affine behaviour).
        let runs = band
            .row2
            .split(|&c| c != b'-')
            .filter(|r| !r.is_empty())
            .count();
        assert_eq!(runs, 1);
    }

    #[test]
    fn banded_affine_rejects_narrow_band() {
        let a = b"ACGTACGTACGT";
        let b = b"ACGTACGT";
        let s = affine_scheme(5, 2);
        assert!(banded_affine(a, b, &s, 1).is_err());
        assert!(banded_affine(a, b, &s, 4).is_ok());
    }

    #[test]
    fn banded_global_over_cap_errors() {
        // Equal-length inputs (so len_diff == 0 <= k clears the
        // narrow-band check) with a huge `k`: band_w = 2k+1 = 16385 and
        // (n+1)·band_w = 8193 · 16385 ≈ 134M cells, past the 64 Mi
        // MAX_DP_CELLS cap. Pre-guard this attempted a ~134M-cell alloc;
        // the guard must reject it as TooLarge before allocating.
        let a = vec![b'A'; 8192];
        let b = vec![b'A'; 8192];
        let s = dna(1);
        let err = banded_global(&a, &b, &s, 8192).unwrap_err();
        assert!(
            matches!(err, AlignError::TooLarge { .. }),
            "expected TooLarge, got {err:?}"
        );
    }

    #[test]
    fn banded_affine_over_cap_errors() {
        // Same oversized band as the linear case; the affine path lays
        // out six layers, so the unchecked allocation is even larger.
        let a = vec![b'A'; 8192];
        let b = vec![b'A'; 8192];
        let s = affine_scheme(5, 2);
        let err = banded_affine(&a, &b, &s, 8192).unwrap_err();
        assert!(
            matches!(err, AlignError::TooLarge { .. }),
            "expected TooLarge, got {err:?}"
        );
    }
}
