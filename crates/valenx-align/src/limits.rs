//! Resource limits for the dynamic-programming routines.
//!
//! Every full-matrix DP in this crate allocates `(n+1)·(m+1)` cells —
//! and several routines allocate that matrix *more than once* (Gotoh
//! and the profile-profile aligner keep six `i32` layers). Two long
//! sequences therefore cost quadratic memory: a pair of ~50 kb
//! sequences is ~2.5 × 10⁹ cells, and at six 4-byte layers that is
//! ~60 GB — an easy out-of-memory / denial-of-service.
//!
//! [`check_dp_size`] is the single chokepoint every DP allocation site
//! calls *before* it allocates. It computes the cell count with a
//! checked multiply (so the size computation itself cannot overflow and
//! wrap to a small number) and rejects anything past
//! [`MAX_DP_CELLS`].

use crate::error::{AlignError, Result};

/// Maximum number of cells a single dynamic-programming matrix may
/// have, i.e. the largest `(n+1)·(m+1)` this crate will allocate.
///
/// Set to **64 Mi cells** — an 8192 × 8192 matrix. The rationale is the
/// *worst-case layer count* at any one site: Gotoh affine alignment and
/// the profile-profile aligner each hold six `i32` matrices, so the
/// peak resident memory at the cap is roughly
/// `64 Mi × 6 × 4 bytes ≈ 1.5 GiB` — large enough to align sequences
/// well past the point where a quadratic-space algorithm is the right
/// tool, yet bounded so a pathological input cannot exhaust RAM.
///
/// Sequences over this size should use a linear-space global aligner
/// ([`crate::pairwise::hirschberg::hirschberg`], O(min(n, m)) space) or
/// a banded routine; the global entry points fall back to Hirschberg
/// automatically, and the others return [`AlignError::TooLarge`].
pub const MAX_DP_CELLS: usize = 64 * 1024 * 1024;

/// Verifies that a `(rows × cols)` DP matrix is within
/// [`MAX_DP_CELLS`], returning the cell count on success.
///
/// `rows` and `cols` are the full matrix dimensions the caller is about
/// to allocate (typically `n + 1` and `m + 1`). The product is computed
/// with [`usize::checked_mul`] so the size check cannot itself overflow
/// and wrap to a deceptively small value; an overflow is reported as
/// [`AlignError::TooLarge`] with `cells = usize::MAX`.
///
/// Returns [`AlignError::TooLarge`] when the matrix would exceed the
/// cap, otherwise `Ok(rows * cols)`.
#[inline]
pub fn check_dp_size(rows: usize, cols: usize) -> Result<usize> {
    match rows.checked_mul(cols) {
        Some(cells) if cells <= MAX_DP_CELLS => Ok(cells),
        Some(cells) => Err(AlignError::too_large(cells, MAX_DP_CELLS)),
        None => Err(AlignError::too_large(usize::MAX, MAX_DP_CELLS)),
    }
}

/// Like [`check_dp_size`] but checks against a caller-supplied cap —
/// used by tests to exercise the guard without allocating a matrix at
/// the production ceiling.
#[inline]
pub fn check_dp_size_with(rows: usize, cols: usize, max: usize) -> Result<usize> {
    match rows.checked_mul(cols) {
        Some(cells) if cells <= max => Ok(cells),
        Some(cells) => Err(AlignError::too_large(cells, max)),
        None => Err(AlignError::too_large(usize::MAX, max)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_matrices_pass() {
        assert_eq!(check_dp_size(10, 10).unwrap(), 100);
        assert_eq!(check_dp_size(1, 1).unwrap(), 1);
        // Exactly at the ceiling is allowed.
        assert_eq!(check_dp_size(MAX_DP_CELLS, 1).unwrap(), MAX_DP_CELLS);
    }

    #[test]
    fn oversized_matrix_is_rejected() {
        // One cell past the cap.
        let err = check_dp_size(MAX_DP_CELLS + 1, 1).unwrap_err();
        match err {
            AlignError::TooLarge { cells, max } => {
                assert_eq!(cells, MAX_DP_CELLS + 1);
                assert_eq!(max, MAX_DP_CELLS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    #[test]
    fn two_large_sequences_would_oom() {
        // ~50 kb × ~50 kb ≈ 2.5e9 cells — comfortably over the cap, and
        // the whole point of the guard. No allocation happens.
        let n = 50_000usize;
        let err = check_dp_size(n + 1, n + 1).unwrap_err();
        assert!(matches!(err, AlignError::TooLarge { .. }));
    }

    #[test]
    fn dimension_product_overflow_is_caught() {
        // A multiply that overflows usize must NOT wrap to a small
        // value and slip past the cap; it is reported as TooLarge with
        // cells == usize::MAX.
        let err = check_dp_size(usize::MAX, 2).unwrap_err();
        match err {
            AlignError::TooLarge { cells, max } => {
                assert_eq!(cells, usize::MAX);
                assert_eq!(max, MAX_DP_CELLS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    #[test]
    fn custom_cap_helper_works() {
        assert_eq!(check_dp_size_with(5, 5, 100).unwrap(), 25);
        assert!(check_dp_size_with(11, 11, 100).is_err());
        assert!(check_dp_size_with(usize::MAX, 2, 100).is_err());
    }
}
