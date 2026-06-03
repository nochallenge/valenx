//! Consensus structure from aligned sequences (RNAalifold-class).
//!
//! When several homologous RNAs are aligned, their *common*
//! secondary structure can be inferred far more reliably than from
//! any single sequence — because evolution leaves a fingerprint:
//! base-paired columns *co-vary* (a mutation on one side of a helix is
//! compensated by a mutation on the other, preserving the pair).
//! RNAalifold (Hofacker *et al.* 2002) folds the *alignment*, scoring
//! each candidate pair of columns by an average folding energy plus a
//! covariation bonus.
//!
//! ## Method
//!
//! For a pair of alignment columns `(i, j)` the score combines:
//!
//! - the **average pairability** — the fraction of sequences in which
//!   columns `i` and `j` hold a canonical pair (a stand-in for the
//!   averaged Turner energy);
//! - a **covariation bonus** — rewards columns that pair *and* differ
//!   between sequences (compensatory mutations), penalises
//!   inconsistent columns where a pair is impossible.
//!
//! A Nussinov-style maximisation over this column-pair score yields
//! the consensus structure. This is a faithful v1 of the RNAalifold
//! objective; the energy term is averaged-pairability rather than the
//! full averaged nearest-neighbor model.

use crate::error::{Result, RnaStructError};
use crate::fold::energy::encode_base;
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::structure::Structure;

/// Weight of the covariation term relative to the pairability term.
pub const COVARIATION_WEIGHT: f64 = 1.0;

/// The result of consensus folding.
#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusResult {
    /// The consensus secondary structure over the alignment columns.
    pub structure: Structure,
    /// The total consensus score the structure achieves.
    pub score: f64,
    /// The number of alignment columns (consensus length).
    pub columns: usize,
    /// The number of sequences in the input alignment.
    pub sequences: usize,
}

/// Computes the consensus structure of a set of aligned RNA
/// sequences (RNAalifold-class).
///
/// Every sequence in `alignment` must have the same length (it is a
/// gapped alignment — `-` marks a gap). The consensus structure spans
/// the alignment columns.
///
/// # Errors
/// - [`RnaStructError::Invalid`] if the alignment is empty or the
///   rows differ in length.
pub fn consensus_structure(alignment: &[&str]) -> Result<ConsensusResult> {
    if alignment.is_empty() {
        return Err(RnaStructError::invalid(
            "alignment",
            "need at least one aligned sequence",
        ));
    }
    let cols = alignment[0].len();
    for (r, row) in alignment.iter().enumerate() {
        if row.len() != cols {
            return Err(RnaStructError::invalid(
                "alignment",
                format!("row {r} has length {} (expected {cols})", row.len()),
            ));
        }
    }
    if cols == 0 {
        return Ok(ConsensusResult {
            structure: Structure::empty(0),
            score: 0.0,
            columns: 0,
            sequences: alignment.len(),
        });
    }

    // Pre-encode the alignment as a column-major matrix of optional
    // base codes (None for a gap / non-RNA char).
    let nseq = alignment.len();
    let mut matrix: Vec<Vec<Option<u8>>> = vec![vec![None; cols]; nseq];
    for (s, row) in alignment.iter().enumerate() {
        for (c, ch) in row.bytes().enumerate() {
            matrix[s][c] = if ch == b'-' {
                None
            } else {
                encode_base(ch)
            };
        }
    }

    // Column-pair score function.
    let pair_score = |i: usize, j: usize| -> f64 {
        column_pair_score(&matrix, nseq, i, j)
    };

    // Nussinov-style maximisation over the column-pair score.
    let mut m = vec![0.0_f64; cols * cols];
    let at = |i: usize, j: usize| i * cols + j;
    for span in 1..cols {
        for i in 0..(cols - span) {
            let j = i + span;
            let mut best = m[at(i + 1, j)].max(m[at(i, j - 1)]);
            if j - i > MIN_HAIRPIN {
                let s = pair_score(i, j);
                if s > 0.0 {
                    let inner = if i + 2 <= j {
                        m[at(i + 1, j - 1)]
                    } else {
                        0.0
                    };
                    best = best.max(inner + s);
                }
            }
            for k in (i + 1)..j {
                let cand = m[at(i, k)] + m[at(k + 1, j)];
                if cand > best {
                    best = cand;
                }
            }
            m[at(i, j)] = best;
        }
    }

    // Traceback.
    let mut partner: Vec<Option<usize>> = vec![None; cols];
    let mut stack = vec![(0usize, cols - 1)];
    let feq = |a: f64, b: f64| (a - b).abs() < 1e-9;
    while let Some((i, j)) = stack.pop() {
        if i >= j {
            continue;
        }
        let here = m[at(i, j)];
        if feq(here, m[at(i + 1, j)]) {
            stack.push((i + 1, j));
            continue;
        }
        if feq(here, m[at(i, j - 1)]) {
            stack.push((i, j - 1));
            continue;
        }
        if j - i > MIN_HAIRPIN {
            let s = pair_score(i, j);
            let inner = if i + 2 <= j {
                m[at(i + 1, j - 1)]
            } else {
                0.0
            };
            if s > 0.0 && feq(here, inner + s) {
                partner[i] = Some(j);
                partner[j] = Some(i);
                if i + 2 <= j {
                    stack.push((i + 1, j - 1));
                }
                continue;
            }
        }
        for k in (i + 1)..j {
            if feq(here, m[at(i, k)] + m[at(k + 1, j)]) {
                stack.push((i, k));
                stack.push((k + 1, j));
                break;
            }
        }
    }

    let structure = Structure::from_partner(partner)?;
    Ok(ConsensusResult {
        structure,
        score: m[at(0, cols - 1)],
        columns: cols,
        sequences: nseq,
    })
}

/// The RNAalifold-style score of pairing alignment columns `i` and
/// `j`: average pairability plus a covariation bonus.
fn column_pair_score(
    matrix: &[Vec<Option<u8>>],
    nseq: usize,
    i: usize,
    j: usize,
) -> f64 {
    let mut canonical = 0usize; // sequences with a canonical pair here
    let mut inconsistent = 0usize; // sequences where the pair is impossible
    let mut distinct_pairs = std::collections::HashSet::new();
    for row in matrix {
        match (row[i], row[j]) {
            (Some(a), Some(b)) => {
                if crate::fold::energy::can_pair_codes(a, b) {
                    canonical += 1;
                    distinct_pairs.insert((a, b));
                } else {
                    inconsistent += 1;
                }
            }
            _ => {
                // a gap on either side: that sequence is neutral
            }
        }
    }
    if canonical == 0 {
        return 0.0;
    }
    // Average pairability in [0, 1].
    let pairability = canonical as f64 / nseq as f64;
    // Covariation bonus: more distinct canonical pair types observed
    // across the alignment => stronger compensatory-mutation signal.
    let covariation = (distinct_pairs.len().saturating_sub(1)) as f64;
    // Inconsistency penalty: a column pair some sequences cannot form.
    let penalty = inconsistent as f64 / nseq as f64;
    pairability + COVARIATION_WEIGHT * covariation - penalty
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_sequence_consensus_is_a_plain_fold() {
        // one sequence -> the consensus is just its maximum-pairing
        // structure
        let aln = ["GGGGAAAACCCC"];
        let r = consensus_structure(&aln).unwrap();
        assert_eq!(r.sequences, 1);
        assert_eq!(r.columns, 12);
        assert!(r.structure.n_pairs() > 0);
    }

    #[test]
    fn covarying_alignment_recovers_the_helix() {
        // three sequences, all forming the same hairpin but with
        // compensatory mutations in the stem
        let aln = [
            "GGGGAAAACCCC",
            "GCGCAAAAGCGC",
            "AUAUAAAAAUAU",
        ];
        let r = consensus_structure(&aln).unwrap();
        // the outer columns should pair in the consensus
        assert!(r.structure.n_pairs() >= 3, "expected a consensus stem");
        assert_eq!(r.structure.partner(0), Some(11));
    }

    #[test]
    fn inconsistent_columns_are_not_paired() {
        // a column pair that no sequence can form should be rejected
        let aln = [
            "AAAAAAAAAAAA",
            "AAAAAAAAAAAA",
        ];
        let r = consensus_structure(&aln).unwrap();
        assert_eq!(r.structure.n_pairs(), 0);
    }

    #[test]
    fn rejects_ragged_alignment() {
        let aln = ["GGGGCCCC", "GGGCCC"];
        assert!(consensus_structure(&aln).is_err());
    }

    #[test]
    fn rejects_empty_alignment() {
        let aln: [&str; 0] = [];
        assert!(consensus_structure(&aln).is_err());
    }

    #[test]
    fn handles_gaps() {
        // gapped alignment: '-' columns are neutral
        let aln = [
            "GGGG-AAAA-CCCC",
            "GGGGAAAAAACCCC",
        ];
        let r = consensus_structure(&aln).unwrap();
        assert_eq!(r.columns, 14);
        assert!(r.structure.is_nested());
    }

    #[test]
    fn consensus_structure_is_always_nested() {
        let aln = [
            "GGGGCCCCAAAGGGGCCCC",
            "GCGCGCGCAAAGCGCGCGC",
        ];
        let r = consensus_structure(&aln).unwrap();
        assert!(r.structure.is_nested());
    }
}
