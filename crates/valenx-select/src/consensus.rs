//! Consensus ranking across orthogonal scoring methods.

use serde::{Deserialize, Serialize};

use crate::error::SelectError;

/// Per-candidate consensus and disagreement across the input methods.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConsensusResult {
    /// Mean fractional rank across methods, in `[0, 1]` (1 = best by consensus).
    pub consensus: Vec<f64>,
    /// Standard deviation of fractional ranks across methods — high means the
    /// methods disagree about this candidate (treat as lower confidence).
    pub disagreement: Vec<f64>,
}

/// Fractional ranks in `[0, 1]` (1 = highest score), ties averaged.
fn fractional_ranks(scores: &[f64]) -> Vec<f64> {
    let n = scores.len();
    if n == 1 {
        return vec![1.0];
    }
    let denom = (n - 1) as f64;
    scores
        .iter()
        .map(|&s| {
            let mut less = 0usize;
            let mut ties = 0usize;
            for &o in scores {
                match o.partial_cmp(&s) {
                    Some(std::cmp::Ordering::Less) => less += 1,
                    Some(std::cmp::Ordering::Equal) => ties += 1, // includes self
                    _ => {}
                }
            }
            // exclude self from the tie group
            (less as f64 + 0.5 * (ties as f64 - 1.0)) / denom
        })
        .collect()
}

/// Aggregate per-method scores (each `method_scores[m][i]` = method `m`'s score
/// for candidate `i`, higher is better) into a consensus rank and a
/// disagreement signal.
///
/// All methods must score the same number of candidates, with finite scores.
pub fn consensus_borda(method_scores: &[Vec<f64>]) -> Result<ConsensusResult, SelectError> {
    if method_scores.is_empty() {
        return Err(SelectError::Empty { what: "methods" });
    }
    let n = method_scores[0].len();
    if n == 0 {
        return Err(SelectError::Empty { what: "candidates" });
    }
    for m in method_scores {
        if m.len() != n {
            return Err(SelectError::Inconsistent {
                what: "method length",
            });
        }
        for &s in m {
            if !s.is_finite() {
                return Err(SelectError::NonFinite { what: "score" });
            }
        }
    }

    let fr: Vec<Vec<f64>> = method_scores.iter().map(|m| fractional_ranks(m)).collect();
    let nm = fr.len() as f64;
    let mut consensus = vec![0.0; n];
    let mut disagreement = vec![0.0; n];
    for i in 0..n {
        let mean = fr.iter().map(|f| f[i]).sum::<f64>() / nm;
        let var = fr.iter().map(|f| (f[i] - mean).powi(2)).sum::<f64>() / nm;
        consensus[i] = mean;
        disagreement[i] = var.sqrt();
    }
    Ok(ConsensusResult {
        consensus,
        disagreement,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fractional_ranks_basic() {
        // 3 is best (1.0), 1 worst (0.0), 2 middle (0.5)
        assert_eq!(fractional_ranks(&[3.0, 1.0, 2.0]), vec![1.0, 0.0, 0.5]);
    }

    #[test]
    fn fractional_ranks_handles_ties() {
        // both tied -> both at the average rank 0.5
        assert_eq!(fractional_ranks(&[5.0, 5.0]), vec![0.5, 0.5]);
    }

    #[test]
    fn agreeing_methods_have_zero_disagreement() {
        let r = consensus_borda(&[vec![3.0, 1.0, 2.0], vec![9.0, 1.0, 5.0]]).unwrap();
        // both rank candidate 0 best, 1 worst, 2 middle
        assert!((r.consensus[0] - 1.0).abs() < 1e-12);
        assert!((r.consensus[2] - 0.5).abs() < 1e-12);
        assert!(r.disagreement.iter().all(|&d| d.abs() < 1e-12));
    }

    #[test]
    fn conflicting_methods_flag_disagreement() {
        // method A: 0 best, 1 worst; method B: 0 worst, 1 best
        let r = consensus_borda(&[vec![3.0, 1.0, 2.0], vec![1.0, 3.0, 2.0]]).unwrap();
        // candidate 0 (ranks 1.0 and 0.0) disagrees; candidate 2 (0.5, 0.5) agrees
        assert!(r.disagreement[0] > r.disagreement[2]);
        assert!((r.disagreement[0] - 0.5).abs() < 1e-12);
        assert!(r.disagreement[2].abs() < 1e-12);
        // consensus pulls the conflicted candidate to the middle
        assert!((r.consensus[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(consensus_borda(&[]).unwrap_err().code(), "empty");
        assert_eq!(
            consensus_borda(&[vec![1.0, 2.0], vec![1.0]])
                .unwrap_err()
                .code(),
            "inconsistent"
        );
        assert_eq!(
            consensus_borda(&[vec![f64::NAN]]).unwrap_err().code(),
            "non_finite"
        );
    }
}
