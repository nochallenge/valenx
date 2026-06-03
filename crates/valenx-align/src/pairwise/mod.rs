//! Pairwise sequence alignment.
//!
//! Every routine takes two byte slices and a
//! [`crate::matrix::ScoringScheme`] and returns an
//! [`result::Alignment`]. The family:
//!
//! | Function | Mode | Gaps | Space |
//! |---|---|---|---|
//! | [`global::needleman_wunsch`] | global | linear | O(nm) |
//! | [`global::gotoh`] | global | affine | O(nm) |
//! | [`local::smith_waterman`] | local | affine | O(nm) |
//! | [`semiglobal::semi_global`] | semi-global | linear | O(nm) |
//! | [`banded::banded_global`] | global | linear | O(nk) |
//! | [`hirschberg::hirschberg`] | global | linear | O(min(n,m)) |
//!
//! The [`result`] module defines the shared [`result::Alignment`],
//! the [`result::Cigar`] type and the [`result::AlignStats`]
//! identity / similarity statistics.

pub mod banded;
pub mod global;
pub mod hirschberg;
pub mod local;
pub mod result;
pub mod semiglobal;

pub use result::{Alignment, AlignStats, Cigar, CigarOp};

use crate::matrix::SubstitutionMatrix;

/// The initial-row score of an affine-gap DP: the cost of a leading run
/// of `len` gap residues, `-(open + extend · len)`, computed in `i64`
/// and clamped to `floor` so a very long `len` cannot overflow `i32`
/// and wrap to a spuriously positive value. `floor` is the routine's
/// `-inf` sentinel (e.g. `i32::MIN / 4`), which keeps headroom for the
/// inner-loop additions while never overflowing.
///
/// (In practice the DP-size cap rejects inputs long enough to overflow,
/// but the init rows are computed before any per-cell arithmetic, so we
/// keep them correct regardless.)
#[inline]
pub(crate) fn affine_init(open: i32, extend: i32, len: usize, floor: i32) -> i32 {
    // `len as i64` would wrap for len > i64::MAX (e.g. usize::MAX -> -1);
    // saturate the conversion so a huge length stays huge.
    let len_i64 = i64::try_from(len).unwrap_or(i64::MAX);
    let cost = (open as i64).saturating_add((extend as i64).saturating_mul(len_i64));
    (-cost).max(floor as i64) as i32
}

/// Convenience: percent identity of two *already-aligned* equal-length
/// gapped rows, in `[0, 1]`. Returns `0.0` if there are no aligned
/// (non-gap) columns. This is the "pairwise identity" feature exposed
/// at module scope; the richer breakdown is
/// [`Alignment::stats`](result::Alignment::stats).
pub fn percent_identity(row1: &[u8], row2: &[u8]) -> f64 {
    let mut ident = 0usize;
    let mut aligned = 0usize;
    for (&a, &b) in row1.iter().zip(row2) {
        if a == b'-' || b == b'-' {
            continue;
        }
        aligned += 1;
        if a.eq_ignore_ascii_case(&b) {
            ident += 1;
        }
    }
    if aligned == 0 {
        0.0
    } else {
        ident as f64 / aligned as f64
    }
}

/// The number of substitutions scoring positively under `matrix`
/// between two aligned rows — BLAST's "positives" count.
pub fn positive_count(row1: &[u8], row2: &[u8], matrix: &SubstitutionMatrix) -> usize {
    row1.iter()
        .zip(row2)
        .filter(|(&a, &b)| a != b'-' && b != b'-' && matrix.score(a, b) > 0)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn module_percent_identity() {
        assert!((percent_identity(b"ACGT", b"ACGT") - 1.0).abs() < 1e-9);
        assert!((percent_identity(b"ACGT", b"ACTT") - 0.75).abs() < 1e-9);
        assert!((percent_identity(b"A-GT", b"ACGT") - 1.0).abs() < 1e-9);
        assert_eq!(percent_identity(b"--", b"AC"), 0.0);
    }

    #[test]
    fn module_positive_count() {
        let m = SubstitutionMatrix::blosum62();
        // L/I is a positive-scoring substitution; A/W is not.
        assert_eq!(positive_count(b"LA", b"IW", &m), 1);
    }

    #[test]
    fn affine_init_matches_naive_for_small_lengths() {
        let floor = i32::MIN / 4;
        for len in 0..50usize {
            let l = len as i32;
            assert_eq!(affine_init(10, 1, len, floor), -(10 + l));
            assert_eq!(affine_init(0, 5, len, floor), -(5 * l));
        }
    }

    #[test]
    fn affine_init_does_not_overflow_for_huge_len() {
        // A naive `-(open + extend * len as i32)` overflows i32 here and
        // could wrap to a positive score; the i64 path clamps to floor.
        let floor = i32::MIN / 4;
        let v = affine_init(11, 2, usize::MAX, floor);
        assert_eq!(v, floor, "huge gap run must clamp to the -inf floor");
        assert!(v < 0, "must stay negative, never wrap positive");
        // Just under the floor boundary stays exact.
        let small = affine_init(0, 1, 1000, floor);
        assert_eq!(small, -1000);
    }
}
