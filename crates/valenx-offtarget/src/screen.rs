//! Screening a candidate against a reference set and summarising the risk.

use serde::{Deserialize, Serialize};

use crate::aa::first_invalid;
use crate::error::OffTargetError;
use crate::similarity::{best_ungapped_identity, kmer_jaccard};

/// One reference's similarity to the candidate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OffTargetHit {
    /// The reference's identifier.
    pub reference_id: String,
    /// Best ungapped fractional identity to the candidate, in `[0, 1]`.
    pub identity: f64,
    /// k-mer Jaccard overlap with the candidate, in `[0, 1]`.
    pub jaccard: f64,
}

/// The result of screening a candidate against a reference set.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OffTargetReport {
    /// Number of references screened.
    pub n_references: usize,
    /// References whose identity met or exceeded the threshold, sorted by
    /// descending identity (ties broken by reference id).
    pub flagged: Vec<OffTargetHit>,
    /// The maximum identity seen over all references, in `[0, 1]`.
    pub max_identity: f64,
}

impl OffTargetReport {
    /// Whether any reference was flagged (a potential off-target liability).
    pub fn has_risk(&self) -> bool {
        !self.flagged.is_empty()
    }
}

/// Screen `candidate` against `references` (`(id, sequence)` pairs), flagging
/// every reference whose best ungapped identity meets `identity_threshold`.
///
/// `k` is the k-mer length for the Jaccard measure (reported alongside each
/// hit). All sequences must be non-empty and all-standard; `references` must be
/// non-empty; `k >= 1`; `identity_threshold` must be finite.
pub fn screen(
    candidate: &str,
    references: &[(&str, &str)],
    k: usize,
    identity_threshold: f64,
) -> Result<OffTargetReport, OffTargetError> {
    if candidate.is_empty() {
        return Err(OffTargetError::Empty { what: "candidate" });
    }
    if let Some((pos, residue)) = first_invalid(candidate) {
        return Err(OffTargetError::InvalidResidue {
            which: "candidate".to_string(),
            pos,
            residue,
        });
    }
    if references.is_empty() {
        return Err(OffTargetError::Empty {
            what: "reference set",
        });
    }
    if k == 0 {
        return Err(OffTargetError::ZeroK);
    }
    if !identity_threshold.is_finite() {
        return Err(OffTargetError::NonFiniteThreshold(identity_threshold));
    }

    let mut hits = Vec::with_capacity(references.len());
    let mut max_identity = 0.0_f64;
    for (id, seq) in references {
        if seq.is_empty() {
            return Err(OffTargetError::Empty { what: "reference" });
        }
        if let Some((pos, residue)) = first_invalid(seq) {
            return Err(OffTargetError::InvalidResidue {
                which: (*id).to_string(),
                pos,
                residue,
            });
        }
        let identity = best_ungapped_identity(candidate, seq)?;
        let jaccard = kmer_jaccard(candidate, seq, k)?;
        if identity > max_identity {
            max_identity = identity;
        }
        hits.push(OffTargetHit {
            reference_id: (*id).to_string(),
            identity,
            jaccard,
        });
    }

    let mut flagged: Vec<OffTargetHit> = hits
        .into_iter()
        .filter(|h| h.identity >= identity_threshold)
        .collect();
    flagged.sort_by(|a, b| {
        b.identity
            .partial_cmp(&a.identity)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.reference_id.cmp(&b.reference_id))
    });

    Ok(OffTargetReport {
        n_references: references.len(),
        flagged,
        max_identity,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs() -> Vec<(&'static str, &'static str)> {
        vec![
            ("self_1", "MKTAYIAKQR"), // identical to the candidate below
            ("near", "MKTAYIAKQQ"),   // one mismatch at the end
            ("unrelated", "WWWWWWWWWW"),
        ]
    }

    #[test]
    fn flags_self_and_near_not_unrelated() {
        let report = screen("MKTAYIAKQR", &refs(), 3, 0.8).unwrap();
        assert_eq!(report.n_references, 3);
        // self_1 (1.0) and near (0.9) are >= 0.8; unrelated (0.0) is not.
        assert_eq!(report.flagged.len(), 2);
        assert_eq!(report.flagged[0].reference_id, "self_1");
        assert!((report.flagged[0].identity - 1.0).abs() < 1e-12);
        assert_eq!(report.flagged[1].reference_id, "near");
        assert!((report.flagged[1].identity - 0.9).abs() < 1e-12);
        assert!((report.max_identity - 1.0).abs() < 1e-12);
        assert!(report.has_risk());
    }

    #[test]
    fn no_risk_when_threshold_high() {
        let report = screen("MKTAYIAKQR", &refs(), 3, 0.95).unwrap();
        // only the identical self_1 clears 0.95
        assert_eq!(report.flagged.len(), 1);
        assert_eq!(report.flagged[0].reference_id, "self_1");
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(screen("", &refs(), 3, 0.8).unwrap_err().code(), "empty");
        assert_eq!(
            screen("MKTAYIAKQR", &[], 3, 0.8).unwrap_err().code(),
            "empty"
        );
        assert_eq!(
            screen("MKTAYIAKQR", &refs(), 0, 0.8).unwrap_err().code(),
            "zero_k"
        );
        assert_eq!(
            screen("MKTAYIAKQR", &refs(), 3, f64::INFINITY)
                .unwrap_err()
                .code(),
            "non_finite_threshold"
        );
        assert_eq!(
            screen("MKTAZIAKQR", &refs(), 3, 0.8).unwrap_err().code(),
            "invalid_residue"
        );
    }

    #[test]
    fn report_serde_round_trips() {
        let report = screen("MKTAYIAKQR", &refs(), 3, 0.8).unwrap();
        let json = serde_json::to_string(&report).unwrap();
        let back: OffTargetReport = serde_json::from_str(&json).unwrap();
        assert_eq!(report, back);
    }
}
