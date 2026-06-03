//! Feature 17 — consensus scoring.
//!
//! Every docking scoring function has its own systematic bias. A
//! ligand that several *different* functions all rank highly is a
//! safer bet than one a single function loves. *Consensus scoring*
//! formalises that: each ligand is scored by several functions, and
//! the per-function results are combined into one consensus ranking.
//!
//! Because the scoring functions live on different energy scales (a
//! Vina kcal/mol is not an AutoDock4 kcal/mol), the robust way to
//! combine them is *rank aggregation* — convert each function's scores
//! to ranks, then combine the ranks. This module offers three
//! aggregation rules via [`ConsensusMethod`]:
//!
//! - [`ConsensusMethod::RankSum`] — sum the per-function ranks
//!   (Borda-style); the ligand with the lowest total rank wins.
//! - [`ConsensusMethod::RankMean`] — the mean rank (identical
//!   ordering to RankSum, expressed per-function).
//! - [`ConsensusMethod::BestRank`] — a ligand's consensus rank is its
//!   *best* rank across functions (optimistic; surfaces compounds at
//!   least one function loves).
//!
//! All three operate on the *scores* — there is no docking here, just
//! aggregation of scores the caller already computed.

use crate::error::{DockScreenError, Result};

/// How per-function ranks are aggregated into a consensus rank.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConsensusMethod {
    /// Sum of per-function ranks (Borda count). Lower total = better.
    RankSum,
    /// Mean of per-function ranks. Same ordering as [`RankSum`](Self::RankSum).
    RankMean,
    /// The best (lowest) rank across functions — optimistic.
    BestRank,
}

/// A ligand's consensus result.
#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusEntry {
    /// The ligand's index in the input score table.
    pub ligand_index: usize,
    /// The ligand's per-function ranks (0 = best), one per scoring
    /// function, in input column order.
    pub per_function_rank: Vec<usize>,
    /// The aggregated consensus value (rank sum / mean / best — its
    /// meaning depends on the [`ConsensusMethod`]).
    pub consensus_value: f64,
}

/// Combine the scores of several scoring functions into a consensus
/// ranking.
///
/// `scores` is a table: `scores[f][l]` is the score function `f`
/// assigned to ligand `l`. Every row must have the same length (one
/// entry per ligand). Lower scores are better (the docking
/// convention).
///
/// Returns the ligands ranked best-consensus-first.
///
/// Returns [`DockScreenError::Invalid`] if `scores` is empty, has an
/// empty function row, or has ragged rows.
pub fn consensus_rank(
    scores: &[Vec<f64>],
    method: ConsensusMethod,
) -> Result<Vec<ConsensusEntry>> {
    if scores.is_empty() {
        return Err(DockScreenError::invalid(
            "scores",
            "consensus needs at least one scoring function",
        ));
    }
    let n_ligands = scores[0].len();
    if n_ligands == 0 {
        return Err(DockScreenError::invalid(
            "scores",
            "scoring functions have no ligand scores",
        ));
    }
    for (f, row) in scores.iter().enumerate() {
        if row.len() != n_ligands {
            return Err(DockScreenError::invalid(
                "scores",
                format!(
                    "scoring-function row {f} has {} entries, expected {n_ligands}",
                    row.len()
                ),
            ));
        }
    }

    // Rank each scoring function's column independently. Rank 0 = the
    // ligand with the lowest (best) score for that function. Ties get
    // the same ("dense") rank in score order.
    let n_functions = scores.len();
    let mut ranks: Vec<Vec<usize>> = vec![vec![0; n_functions]; n_ligands];
    for (f, row) in scores.iter().enumerate() {
        let mut order: Vec<usize> = (0..n_ligands).collect();
        order.sort_by(|&a, &b| {
            row[a]
                .partial_cmp(&row[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (rank, &ligand) in order.iter().enumerate() {
            ranks[ligand][f] = rank;
        }
    }

    // Aggregate.
    let mut entries: Vec<ConsensusEntry> = (0..n_ligands)
        .map(|l| {
            let per = &ranks[l];
            let value = match method {
                ConsensusMethod::RankSum => per.iter().sum::<usize>() as f64,
                ConsensusMethod::RankMean => {
                    per.iter().sum::<usize>() as f64 / n_functions as f64
                }
                ConsensusMethod::BestRank => {
                    per.iter().copied().min().unwrap_or(0) as f64
                }
            };
            ConsensusEntry {
                ligand_index: l,
                per_function_rank: per.clone(),
                consensus_value: value,
            }
        })
        .collect();

    // Lower consensus value = better; ties broken by ligand index for
    // a deterministic order.
    entries.sort_by(|a, b| {
        a.consensus_value
            .partial_cmp(&b.consensus_value)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.ligand_index.cmp(&b.ligand_index))
    });
    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_empty_or_ragged_tables() {
        assert!(consensus_rank(&[], ConsensusMethod::RankSum).is_err());
        assert!(consensus_rank(&[vec![]], ConsensusMethod::RankSum).is_err());
        // Ragged: row 1 is shorter.
        let ragged = vec![vec![1.0, 2.0, 3.0], vec![1.0, 2.0]];
        assert!(consensus_rank(&ragged, ConsensusMethod::RankSum).is_err());
    }

    #[test]
    fn unanimous_winner_ranks_first() {
        // Two functions, three ligands. Ligand 1 is best for both.
        let scores = vec![
            vec![-5.0, -9.0, -3.0], // function A
            vec![-6.0, -8.0, -4.0], // function B
        ];
        let ranked = consensus_rank(&scores, ConsensusMethod::RankSum).unwrap();
        assert_eq!(ranked[0].ligand_index, 1);
        // Ligand 1 has rank 0 in both → rank sum 0.
        assert_eq!(ranked[0].consensus_value, 0.0);
    }

    #[test]
    fn rank_sum_resolves_a_disagreement() {
        // Function A loves ligand 0; function B loves ligand 1.
        // Ligand 2 is mediocre for both. RankSum: 0→0+2=2, 1→2+0=2,
        // 2→1+1=2 — a three-way tie broken by index.
        let scores = vec![
            vec![-9.0, -3.0, -5.0],
            vec![-3.0, -9.0, -5.0],
        ];
        let ranked = consensus_rank(&scores, ConsensusMethod::RankSum).unwrap();
        assert_eq!(ranked.len(), 3);
        // All three have the same rank sum.
        assert!(ranked.iter().all(|e| (e.consensus_value - 2.0).abs() < 1e-9));
        // Deterministic tie-break by ligand index.
        assert_eq!(ranked[0].ligand_index, 0);
    }

    #[test]
    fn best_rank_surfaces_a_polarising_compound() {
        // Ligand 0 is the favourite of function A (rank 0) but the
        // worst for function B. BestRank gives it consensus 0.
        let scores = vec![
            vec![-9.0, -5.0, -4.0],
            vec![-3.0, -8.0, -9.0],
        ];
        let ranked = consensus_rank(&scores, ConsensusMethod::BestRank).unwrap();
        // Ligand 0 has best rank 0 (from function A).
        let l0 = ranked.iter().find(|e| e.ligand_index == 0).unwrap();
        assert_eq!(l0.consensus_value, 0.0);
    }

    #[test]
    fn rank_mean_orders_identically_to_rank_sum() {
        let scores = vec![
            vec![-5.0, -9.0, -3.0],
            vec![-6.0, -8.0, -4.0],
        ];
        let by_sum = consensus_rank(&scores, ConsensusMethod::RankSum).unwrap();
        let by_mean = consensus_rank(&scores, ConsensusMethod::RankMean).unwrap();
        let sum_order: Vec<usize> = by_sum.iter().map(|e| e.ligand_index).collect();
        let mean_order: Vec<usize> = by_mean.iter().map(|e| e.ligand_index).collect();
        assert_eq!(sum_order, mean_order);
    }

    #[test]
    fn per_function_ranks_are_recorded() {
        let scores = vec![vec![-5.0, -9.0], vec![-8.0, -6.0]];
        let ranked = consensus_rank(&scores, ConsensusMethod::RankSum).unwrap();
        // Every entry carries one rank per scoring function.
        for e in &ranked {
            assert_eq!(e.per_function_rank.len(), 2);
        }
    }

    #[test]
    fn single_function_consensus_is_just_that_function() {
        // With one function, consensus order is the score order.
        let scores = vec![vec![-3.0, -9.0, -5.0]];
        let ranked = consensus_rank(&scores, ConsensusMethod::RankSum).unwrap();
        assert_eq!(ranked[0].ligand_index, 1); // -9.0 is best
        assert_eq!(ranked[2].ligand_index, 0); // -3.0 is worst
    }
}
