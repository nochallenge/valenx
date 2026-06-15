//! The end-to-end selection funnel: consensus → diversify → top-`N`.

use serde::{Deserialize, Serialize};

use crate::consensus::consensus_borda;
use crate::diversity::sphere_exclusion_select;
use crate::error::SelectError;

/// One candidate entering the funnel. Scores from each orthogonal method, a
/// feature vector for diversity, and the safety/confidence metadata to carry
/// through to the shortlist.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Candidate {
    /// Stable candidate identifier.
    pub id: String,
    /// One score per scoring method (higher is better), same order for every
    /// candidate.
    pub method_scores: Vec<f64>,
    /// Feature vector used for diversity selection.
    pub features: Vec<f64>,
    /// Optional calibrated confidence (e.g. from `valenx-calibrate`).
    pub calibrated_confidence: Option<f64>,
    /// Safety flags raised upstream (e.g. off-target / immunogenicity hits).
    /// An empty list is **not** an assertion of safety.
    pub safety_flags: Vec<String>,
}

/// One row of the selected shortlist.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ShortlistEntry {
    /// The candidate's identifier.
    pub id: String,
    /// Consensus rank in `[0, 1]` (1 = best by agreement across methods).
    pub consensus_score: f64,
    /// Cross-method disagreement (higher = treat as lower confidence).
    pub disagreement: f64,
    /// Carried-through calibrated confidence, if any.
    pub calibrated_confidence: Option<f64>,
    /// Carried-through safety flags.
    pub safety_flags: Vec<String>,
    /// Position in the diverse selection (0 = picked first).
    pub diversity_rank: usize,
}

/// The funnel's output: a ranked, diverse shortlist plus the size it was cut to.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Shortlist {
    /// The selected entries, in diversity-selection order.
    pub entries: Vec<ShortlistEntry>,
    /// The requested cut size `n` (the shortlist may be shorter if diversity is
    /// exhausted).
    pub requested_n: usize,
}

/// Run the funnel: consensus-rank the candidates across their methods, then take
/// a diverse top-`n` by sphere exclusion at `radius`, carrying confidence and
/// safety flags through.
///
/// `n` is just where the operator cuts (5, 20, 100); the underlying ranking is
/// the same. Every candidate must carry the same number of method scores and
/// the same feature dimension.
pub fn select_shortlist(
    candidates: &[Candidate],
    n: usize,
    radius: f64,
) -> Result<Shortlist, SelectError> {
    if candidates.is_empty() {
        return Err(SelectError::Empty { what: "candidates" });
    }
    if n == 0 {
        return Err(SelectError::ZeroN);
    }
    let n_methods = candidates[0].method_scores.len();
    if n_methods == 0 {
        return Err(SelectError::Empty { what: "methods" });
    }
    for c in candidates {
        if c.method_scores.len() != n_methods {
            return Err(SelectError::Inconsistent {
                what: "method length",
            });
        }
    }

    // Transpose to per-method score columns, then consensus-rank.
    let per_method: Vec<Vec<f64>> = (0..n_methods)
        .map(|m| candidates.iter().map(|c| c.method_scores[m]).collect())
        .collect();
    let consensus = consensus_borda(&per_method)?;

    // Diverse top-n, ordered by consensus.
    let features: Vec<Vec<f64>> = candidates.iter().map(|c| c.features.clone()).collect();
    let selected = sphere_exclusion_select(&features, &consensus.consensus, n, radius)?;

    let entries = selected
        .iter()
        .enumerate()
        .map(|(rank, &i)| ShortlistEntry {
            id: candidates[i].id.clone(),
            consensus_score: consensus.consensus[i],
            disagreement: consensus.disagreement[i],
            calibrated_confidence: candidates[i].calibrated_confidence,
            safety_flags: candidates[i].safety_flags.clone(),
            diversity_rank: rank,
        })
        .collect();

    Ok(Shortlist {
        entries,
        requested_n: n,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cand(id: &str, scores: Vec<f64>, feat: f64, flags: Vec<&str>) -> Candidate {
        Candidate {
            id: id.to_string(),
            method_scores: scores,
            features: vec![feat],
            calibrated_confidence: Some(0.5),
            safety_flags: flags.into_iter().map(|s| s.to_string()).collect(),
        }
    }

    #[test]
    fn funnel_selects_diverse_top_n_and_carries_flags() {
        // Two near-duplicate strong candidates near feature 0, one near 10.
        let cands = vec![
            cand("A", vec![0.95, 0.9], 0.0, vec!["GDF11_crossreactivity"]),
            cand("A2", vec![0.94, 0.9], 0.1, vec![]),
            cand("B", vec![0.7, 0.7], 10.0, vec![]),
        ];
        let sl = select_shortlist(&cands, 2, 1.0).unwrap();
        assert_eq!(sl.requested_n, 2);
        // Picks the best of the near-duplicate cluster (A) then the distant B —
        // not A and its near-twin A2.
        let ids: Vec<&str> = sl.entries.iter().map(|e| e.id.as_str()).collect();
        assert_eq!(ids, vec!["A", "B"]);
        // The safety flag rides along to the shortlist.
        assert_eq!(sl.entries[0].safety_flags, vec!["GDF11_crossreactivity"]);
        assert_eq!(sl.entries[0].diversity_rank, 0);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(select_shortlist(&[], 2, 1.0).unwrap_err().code(), "empty");
        let c = vec![cand("A", vec![1.0], 0.0, vec![])];
        assert_eq!(select_shortlist(&c, 0, 1.0).unwrap_err().code(), "zero_n");
        let mixed = vec![
            cand("A", vec![1.0, 2.0], 0.0, vec![]),
            cand("B", vec![1.0], 1.0, vec![]),
        ];
        assert_eq!(
            select_shortlist(&mixed, 2, 1.0).unwrap_err().code(),
            "inconsistent"
        );
    }

    #[test]
    fn shortlist_serde_round_trips() {
        let c = vec![cand("A", vec![1.0, 1.0], 0.0, vec!["flag"])];
        let sl = select_shortlist(&c, 1, 1.0).unwrap();
        let j = serde_json::to_string(&sl).unwrap();
        let back: Shortlist = serde_json::from_str(&j).unwrap();
        assert_eq!(sl, back);
    }
}
