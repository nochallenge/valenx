//! **Feature 1 — template search.**
//!
//! Given a query sequence and a set of candidate templates (each a
//! sequence with a known structure), score every candidate by global
//! sequence alignment to the query and rank them. The best-scoring
//! candidate above an identity threshold is the modelling template.
//!
//! Scoring reuses `valenx-align`'s Gotoh affine-gap global alignment
//! with BLOSUM62 — the same alignment a real comparative-modelling
//! pipeline uses to pick a template. A production pipeline would
//! additionally screen a profile-HMM / PSI-BLAST search against the
//! whole PDB; here the caller supplies the candidate set and this
//! module does the scoring and ranking.

use valenx_align::pairwise::global::gotoh;
use valenx_align::ScoringScheme;

use crate::error::{Result, StructPredictError};

/// A candidate modelling template: a label, its sequence, and the
/// PDB / mmCIF coordinate text of its structure.
#[derive(Clone, Debug, PartialEq)]
pub struct TemplateCandidate {
    /// A human-readable identifier (a PDB id, a file name, …).
    pub id: String,
    /// The template's one-letter amino-acid sequence.
    pub sequence: String,
    /// The template's structure as PDB or mmCIF text.
    pub structure_text: String,
}

impl TemplateCandidate {
    /// Builds a candidate.
    pub fn new(
        id: impl Into<String>,
        sequence: impl Into<String>,
        structure_text: impl Into<String>,
    ) -> Self {
        TemplateCandidate {
            id: id.into(),
            sequence: sequence.into(),
            structure_text: structure_text.into(),
        }
    }
}

/// A scored template hit — one entry of a ranked search result.
#[derive(Clone, Debug, PartialEq)]
pub struct TemplateHit {
    /// The candidate's index in the input slice.
    pub index: usize,
    /// The candidate's identifier.
    pub id: String,
    /// The alignment score (Gotoh, BLOSUM62).
    pub score: i32,
    /// Sequence identity of the query-template alignment, `[0, 1]`.
    pub identity: f64,
    /// Alignment coverage: aligned columns / query length, `[0, 1]`.
    pub coverage: f64,
}

impl TemplateHit {
    /// A composite confidence in `[0, 1]`: the geometric blend of
    /// identity and coverage. Higher is a more reliable template.
    pub fn confidence(&self) -> f64 {
        (self.identity * self.coverage).max(0.0).sqrt()
    }
}

/// Searches a set of candidate templates against a query sequence and
/// returns the hits ranked best-first.
///
/// Each candidate is globally aligned to the query; the hit carries
/// the score, sequence identity and coverage. Hits are sorted by
/// alignment score descending.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty query or an empty
/// candidate set; [`StructPredictError::NotFound`] if no candidate
/// has a non-empty sequence.
pub fn search_templates(query: &str, candidates: &[TemplateCandidate]) -> Result<Vec<TemplateHit>> {
    let query = query.trim();
    if query.is_empty() {
        return Err(StructPredictError::invalid("query", "empty sequence"));
    }
    if candidates.is_empty() {
        return Err(StructPredictError::invalid(
            "candidates",
            "no template candidates supplied",
        ));
    }
    let scheme = ScoringScheme::blosum62_default();
    let q = query.as_bytes();
    let mut hits = Vec::new();
    for (index, cand) in candidates.iter().enumerate() {
        let t = cand.sequence.trim();
        if t.is_empty() {
            continue;
        }
        let aln = gotoh(q, t.as_bytes(), &scheme)
            .map_err(|e| StructPredictError::invalid("alignment", e.to_string()))?;
        let mut identical = 0usize;
        let mut aligned = 0usize;
        for (&a, &b) in aln.row1.iter().zip(&aln.row2) {
            if a != b'-' && b != b'-' {
                aligned += 1;
                if a == b {
                    identical += 1;
                }
            }
        }
        let identity = if aligned > 0 {
            identical as f64 / aligned as f64
        } else {
            0.0
        };
        let coverage = aligned as f64 / q.len() as f64;
        hits.push(TemplateHit {
            index,
            id: cand.id.clone(),
            score: aln.score,
            identity,
            coverage,
        });
    }
    if hits.is_empty() {
        return Err(StructPredictError::not_found(
            "template",
            "every candidate had an empty sequence",
        ));
    }
    hits.sort_by(|a, b| b.score.cmp(&a.score));
    Ok(hits)
}

/// Picks the single best template above an identity floor.
///
/// Runs [`search_templates`] and returns the top hit, but only if its
/// sequence identity is at least `min_identity`. Below the floor the
/// model would be unreliable, so the function fails rather than
/// returning a misleading template.
///
/// # Errors
/// [`StructPredictError::NotFound`] if no template clears
/// `min_identity`; propagates [`search_templates`] errors.
pub fn best_template(
    query: &str,
    candidates: &[TemplateCandidate],
    min_identity: f64,
) -> Result<TemplateHit> {
    let hits = search_templates(query, candidates)?;
    let best = hits.into_iter().next().expect("non-empty after search");
    if best.identity < min_identity {
        return Err(StructPredictError::not_found(
            "template",
            format!(
                "best identity {:.1}% below floor {:.1}%",
                best.identity * 100.0,
                min_identity * 100.0
            ),
        ));
    }
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_template_scores_top() {
        let q = "ACDEFGHIKLMNPQRSTVWY";
        let cands = vec![
            TemplateCandidate::new("far", "WWWWWWWWWWWWWWWWWWWW", ""),
            TemplateCandidate::new("self", q, ""),
        ];
        let hits = search_templates(q, &cands).expect("search");
        assert_eq!(hits[0].id, "self");
        assert!((hits[0].identity - 1.0).abs() < 1e-9);
        assert!(hits[0].confidence() > 0.9);
    }

    #[test]
    fn empty_query_rejected() {
        assert!(search_templates("", &[]).is_err());
    }

    #[test]
    fn best_template_respects_floor() {
        let q = "ACDEFGHIKLMNPQRSTVWY";
        let cands = vec![TemplateCandidate::new("weak", "ACDEFWWWWWWWWWWWWWWW", "")];
        // ~25% identity — below a 60% floor.
        assert!(best_template(q, &cands, 0.6).is_err());
        // ... but clears a 10% floor.
        assert!(best_template(q, &cands, 0.1).is_ok());
    }
}
