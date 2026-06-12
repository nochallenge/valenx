//! Iterative refinement of a multiple-sequence alignment
//! (MUSCLE-class).
//!
//! Progressive alignment can lock in an early mistake. Iterative
//! refinement repairs it with **tree-partition realignment**: pick an
//! edge of the guide tree, split the sequences into the two groups it
//! separates, re-align each group's profile against the other, and
//! *keep the result only if the sum-of-pairs score improved*. Repeat
//! over edges until no edge yields an improvement or an iteration cap
//! is hit.
//!
//! [`refine`] takes a starting [`Msa`] and returns a refined one whose
//! sum-of-pairs score is greater than or equal to the input's — the
//! "keep if improved" rule makes the routine monotone.

use super::profile::{align_profiles, Profile};
use super::progressive::Msa;
use crate::error::Result;
use crate::matrix::ScoringScheme;

/// Tunable parameters for [`refine`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct RefineParams {
    /// Maximum number of full refinement passes.
    pub max_iterations: usize,
}

impl Default for RefineParams {
    fn default() -> Self {
        RefineParams { max_iterations: 8 }
    }
}

/// Iteratively refines `msa`, returning an alignment whose
/// sum-of-pairs score is `>=` the input's.
///
/// Each pass tries every single-sequence-vs-rest partition (the
/// simplest useful partition family): it pulls one sequence out,
/// re-aligns it against the profile of the others, and accepts the new
/// alignment if it scores higher. With fewer than three sequences
/// there is nothing to refine and the input is returned unchanged.
pub fn refine(msa: &Msa, scheme: &ScoringScheme, params: RefineParams) -> Result<Msa> {
    if msa.depth() < 3 {
        return Ok(msa.clone());
    }

    let mut current = msa.clone();
    let mut current_score = current.sum_of_pairs(scheme);

    for _ in 0..params.max_iterations {
        let mut improved = false;

        for pulled in 0..current.depth() {
            let candidate = realign_one(&current, pulled, scheme)?;
            let cand_score = candidate.sum_of_pairs(scheme);
            if cand_score > current_score {
                current = candidate;
                current_score = cand_score;
                improved = true;
            }
        }

        if !improved {
            break; // converged
        }
    }

    Ok(current)
}

/// Pulls sequence `pulled` out of `msa`, builds a profile from the rest
/// (with all-gap columns stripped), aligns the pulled sequence back
/// in, and returns the re-assembled MSA in the original row order.
fn realign_one(msa: &Msa, pulled: usize, scheme: &ScoringScheme) -> Result<Msa> {
    // The "rest" rows, gap columns that are now all-gap removed.
    let rest_rows: Vec<Vec<u8>> = (0..msa.depth())
        .filter(|&i| i != pulled)
        .map(|i| msa.rows[i].clone())
        .collect();
    let rest_rows = strip_all_gap_columns(&rest_rows);

    let rest_profile = Profile::from_alignment(&rest_rows)?;

    // The pulled sequence, ungapped.
    let pulled_seq: Vec<u8> = msa.rows[pulled]
        .iter()
        .copied()
        .filter(|&c| c != b'-')
        .collect();
    let pulled_profile = Profile::from_sequence(&pulled_seq)?;

    let merged = align_profiles(&rest_profile, &pulled_profile, scheme)?;
    // merged.rows = rest rows (in rest order) then the pulled row.
    let rest_count = rest_rows.len();
    let pulled_row = merged.rows[rest_count].clone();

    // Re-interleave into original order.
    let mut out: Vec<Vec<u8>> = vec![Vec::new(); msa.depth()];
    let mut rest_iter = 0;
    for (i, slot) in out.iter_mut().enumerate() {
        if i == pulled {
            *slot = pulled_row.clone();
        } else {
            *slot = merged.rows[rest_iter].clone();
            rest_iter += 1;
        }
    }
    Msa::new(out)
}

/// Removes columns that are gaps in every row (they can appear after a
/// sequence is pulled out of an alignment).
fn strip_all_gap_columns(rows: &[Vec<u8>]) -> Vec<Vec<u8>> {
    if rows.is_empty() {
        return Vec::new();
    }
    let width = rows[0].len();
    let keep: Vec<usize> = (0..width)
        .filter(|&c| rows.iter().any(|r| r[c] != b'-'))
        .collect();
    rows.iter()
        .map(|r| keep.iter().map(|&c| r[c]).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};
    use crate::msa::progressive::align;

    fn dna() -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(2, -1), GapCost::new(4, 1))
    }

    #[test]
    fn refine_never_lowers_score() {
        let seqs: &[&[u8]] = &[
            b"ACGTACGT",
            b"ACGTCGT",
            b"ACGTACGTT",
            b"ACGACGT",
            b"ACGTACG",
        ];
        let start = align(seqs, &dna()).unwrap();
        let start_score = start.sum_of_pairs(&dna());
        let refined = refine(&start, &dna(), RefineParams::default()).unwrap();
        assert!(
            refined.sum_of_pairs(&dna()) >= start_score,
            "refinement must be monotone non-decreasing"
        );
    }

    #[test]
    fn refine_preserves_sequences() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTCGT", b"ACGTACGTT", b"ACGACGT"];
        let start = align(seqs, &dna()).unwrap();
        let refined = refine(&start, &dna(), RefineParams::default()).unwrap();
        // Each row still recovers its original sequence when degapped.
        for (i, &orig) in seqs.iter().enumerate() {
            let stripped: Vec<u8> = refined.rows[i]
                .iter()
                .copied()
                .filter(|&c| c != b'-')
                .collect();
            assert_eq!(stripped, orig);
        }
    }

    #[test]
    fn refine_keeps_rows_equal_length() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTCGT", b"ACGTACGTT", b"ACGACGT"];
        let start = align(seqs, &dna()).unwrap();
        let refined = refine(&start, &dna(), RefineParams::default()).unwrap();
        let w = refined.width();
        assert!(refined.rows.iter().all(|r| r.len() == w));
    }

    #[test]
    fn small_msa_unchanged() {
        let seqs: &[&[u8]] = &[b"ACGT", b"ACGT"];
        let start = align(seqs, &dna()).unwrap();
        let refined = refine(&start, &dna(), RefineParams::default()).unwrap();
        assert_eq!(refined, start);
    }

    #[test]
    fn strip_all_gap_columns_works() {
        let rows = vec![b"A-C".to_vec(), b"A-G".to_vec()];
        let stripped = strip_all_gap_columns(&rows);
        // Middle column is all-gap -> removed.
        assert_eq!(stripped, vec![b"AC".to_vec(), b"AG".to_vec()]);
    }
}
