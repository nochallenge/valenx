//! The PPI confidence score and the interactome screen.
//!
//! [`PpiScore`] folds the heterogeneous evidence for one candidate
//! protein-protein interaction into one comparable `[0, 1]` value, in
//! the exact style of [`valenx_binder_score`]: every component stays
//! visible, the fusion is a monotone weighted mean, and **the result is
//! never a validated "interacts" verdict** — [`PpiScore::requires_review`]
//! is always `true`.

use valenx_align::matrix::ScoringScheme;
use valenx_align::msa::{align as align_msa, Msa};
use valenx_biostruct::structure::Chain;

use crate::coevolution::{predict_contacts, CoevolutionResult, PairedMsa};
use crate::complementarity::{interface_complementarity, Complementarity};
use crate::error::PpiError;

/// A fused PPI-confidence score with its components kept visible.
///
/// Mirrors `valenx_binder_score::BinderScore`: a `[0, 1]` aggregate plus
/// the per-channel signals it was built from. The `complementarity`
/// channel is `None` when no structures were supplied (sequence-only
/// mode).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PpiScore {
    /// The fused confidence in `[0, 1]` (higher = stronger evidence of
    /// interaction). **Not** a probability and **not** a verdict.
    pub value: f64,
    /// The sequence coevolution signal in `[0, 1]`
    /// (see [`CoevolutionResult::signal`]).
    pub coevolution: f64,
    /// The geometric interface-complementarity signal in `[0, 1]`, or
    /// `None` when no structures were available for both partners.
    pub complementarity: Option<f64>,
}

impl PpiScore {
    /// Always `true`: a PPI score *ranks candidate interactions for a
    /// human to triage*; it never confirms that two proteins interact.
    /// Any downstream "these proteins bind" claim needs orthogonal
    /// evidence and wet-lab validation. Mirrors
    /// `BinderScore::requires_review`.
    pub fn requires_review(&self) -> bool {
        true
    }
}

/// Per-channel weights for [`PpiScore`] fusion (default `1.0`; `0` drops
/// the channel). Mirrors `valenx_binder_score::BinderWeights`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PpiWeights {
    /// Weight on the coevolution channel.
    pub coevolution: f64,
    /// Weight on the complementarity channel (ignored when absent).
    pub complementarity: f64,
}

impl Default for PpiWeights {
    fn default() -> Self {
        Self {
            coevolution: 1.0,
            complementarity: 1.0,
        }
    }
}

fn check_weight(what: &'static str, w: f64) -> Result<(), PpiError> {
    if !w.is_finite() || w < 0.0 {
        return Err(PpiError::BadWeight { what, value: w });
    }
    Ok(())
}

/// Optional structural evidence for a chain pair: the two chains whose
/// geometric interface complementarity reinforces the coevolution
/// signal. Supply `None` for sequence-only scoring.
#[derive(Clone, Copy, Debug)]
pub struct StructuralEvidence<'a> {
    /// Chain A coordinates.
    pub chain_a: &'a Chain,
    /// Chain B coordinates.
    pub chain_b: &'a Chain,
}

/// Score one candidate PPI from a paired MSA, optionally reinforced by
/// structural complementarity, using default weights.
///
/// This is the headline `score_pair(chain_a, chain_b)` entry point: pass
/// the paired alignment of the two chains' orthologues, and (optionally)
/// their structures. The result is a [`PpiScore`] whose
/// [`requires_review`](PpiScore::requires_review) is always `true`.
///
/// # Errors
/// Propagates coevolution / complementarity validation errors
/// fail-loud (empty or mismatched MSA, missing structure).
pub fn score_pair(
    paired: &PairedMsa,
    structures: Option<StructuralEvidence<'_>>,
) -> Result<PpiScore, PpiError> {
    score_pair_weighted(paired, structures, &PpiWeights::default())
}

/// Like [`score_pair`] but with caller-supplied [`PpiWeights`].
///
/// The fusion is the weighted mean of the present, non-zero-weight
/// channels — monotone in each (a stronger coevolution or
/// complementarity signal never lowers the score), exactly as
/// `valenx_binder_score::score` is monotone in its channels.
///
/// # Errors
/// - [`PpiError::BadWeight`] for a negative / non-finite weight.
/// - [`PpiError::NonFinite`] if a channel produced a non-finite signal.
/// - propagates coevolution / complementarity errors.
pub fn score_pair_weighted(
    paired: &PairedMsa,
    structures: Option<StructuralEvidence<'_>>,
    weights: &PpiWeights,
) -> Result<PpiScore, PpiError> {
    check_weight("coevolution", weights.coevolution)?;
    check_weight("complementarity", weights.complementarity)?;

    let coev: CoevolutionResult = predict_contacts(paired)?;
    let coevolution = coev.signal();
    if !coevolution.is_finite() {
        return Err(PpiError::NonFinite {
            what: "coevolution",
        });
    }

    let complementarity: Option<f64> = match structures {
        Some(ev) => {
            let c: Complementarity = interface_complementarity(ev.chain_a, ev.chain_b)?;
            if !c.value.is_finite() {
                return Err(PpiError::NonFinite {
                    what: "complementarity",
                });
            }
            Some(c.value)
        }
        None => None,
    };

    // Weighted mean of present channels (binder-score style).
    let mut acc = 0.0;
    let mut wsum = 0.0;
    if weights.coevolution > 0.0 {
        acc += coevolution * weights.coevolution;
        wsum += weights.coevolution;
    }
    if let Some(c) = complementarity {
        if weights.complementarity > 0.0 {
            acc += c * weights.complementarity;
            wsum += weights.complementarity;
        }
    }
    // Coevolution is always present and >= 0 weight is rejected above
    // only when negative; a zero coevolution weight with no structures
    // would leave nothing to score.
    if wsum <= 0.0 {
        return Err(PpiError::invalid(
            "weights",
            "all present channels had zero weight",
        ));
    }
    let value = acc / wsum;

    Ok(PpiScore {
        value,
        coevolution,
        complementarity,
    })
}

/// One row of an interactome screen: a host protein, a pathogen protein,
/// and their [`PpiScore`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interaction {
    /// Index of the host protein in the screen's `host` list.
    pub host: usize,
    /// Index of the pathogen protein in the screen's `pathogen` list.
    pub pathogen: usize,
    /// The fused score for this pair.
    pub score: PpiScore,
}

/// A ranked host × pathogen interactome screen.
#[derive(Clone, Debug, PartialEq)]
pub struct RankedInteractions {
    /// All scored host-pathogen pairs, sorted by descending
    /// [`PpiScore::value`].
    pub ranked: Vec<Interaction>,
    /// Number of host proteins screened.
    pub n_host: usize,
    /// Number of pathogen proteins screened.
    pub n_pathogen: usize,
}

impl RankedInteractions {
    /// The top-`k` ranked interactions.
    pub fn top(&self, k: usize) -> &[Interaction] {
        &self.ranked[..k.min(self.ranked.len())]
    }

    /// Always `true`: an interactome screen *prioritises pairs for
    /// follow-up*, it never declares any pair a confirmed interaction.
    pub fn requires_review(&self) -> bool {
        true
    }
}

/// One protein entry in an interactome screen: a label plus the aligned
/// orthologue rows used to build paired MSAs. All entries on both sides
/// **must share the same row order** (organism `k` is the same species
/// across every entry) for the pairing to be biologically meaningful.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScreenEntry {
    /// A human-readable identifier for this protein.
    pub label: String,
    /// Aligned orthologue rows (one per organism, in the shared order).
    pub msa: Msa,
}

impl ScreenEntry {
    /// Build a screen entry from a label and an MSA.
    pub fn new(label: impl Into<String>, msa: Msa) -> Self {
        ScreenEntry {
            label: label.into(),
            msa,
        }
    }
}

/// Run a sequence-only host × pathogen all-vs-all interactome screen.
///
/// For each `(host i, pathogen j)` it pairs their orthologue MSAs and
/// computes a coevolution-only [`PpiScore`], then ranks every pair by
/// descending score. Entries that cannot be paired (depth mismatch,
/// too-few sequences, etc.) propagate the error fail-loud rather than
/// being silently dropped — a partial screen that hides failures would
/// be a wrong number.
///
/// Structural complementarity is intentionally *not* part of the
/// all-vs-all screen (it needs a docked pose per pair); use
/// [`score_pair`] with [`StructuralEvidence`] to reinforce a shortlisted
/// pair afterwards.
///
/// # Errors
/// Propagates any [`PpiError`] from pairing or scoring a pair.
pub fn interactome_screen(
    host: &[ScreenEntry],
    pathogen: &[ScreenEntry],
) -> Result<RankedInteractions, PpiError> {
    if host.is_empty() || pathogen.is_empty() {
        return Err(PpiError::invalid(
            "screen",
            "host and pathogen lists must both be non-empty",
        ));
    }
    let mut ranked = Vec::with_capacity(host.len() * pathogen.len());
    for (hi, h) in host.iter().enumerate() {
        for (pj, p) in pathogen.iter().enumerate() {
            let paired = PairedMsa::new(h.msa.clone(), p.msa.clone())?;
            let score = score_pair(&paired, None)?;
            ranked.push(Interaction {
                host: hi,
                pathogen: pj,
                score,
            });
        }
    }
    ranked.sort_by(|a, b| {
        b.score
            .value
            .partial_cmp(&a.score.value)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.host.cmp(&b.host))
            .then(a.pathogen.cmp(&b.pathogen))
    });
    Ok(RankedInteractions {
        ranked,
        n_host: host.len(),
        n_pathogen: pathogen.len(),
    })
}

/// Convenience: build a paired MSA by aligning two sets of raw
/// (unaligned) orthologue sequences with the protein default scheme,
/// then validating the pairing.
///
/// `a_seqs` and `b_seqs` must have the same length (one orthologue pair
/// per organism, same order). Each side is aligned independently with
/// [`align_msa`] under [`ScoringScheme::blosum62_default`]; the two
/// resulting alignments are paired by row.
///
/// # Errors
/// - [`PpiError::Invalid`] if the two input lists differ in length.
/// - alignment failures are surfaced as [`PpiError::Invalid`] tagged
///   `"align"` (the underlying [`valenx_align`] error message is
///   carried through).
/// - [`PairedMsa::new`] validation errors.
pub fn build_paired_msa(a_seqs: &[&[u8]], b_seqs: &[&[u8]]) -> Result<PairedMsa, PpiError> {
    if a_seqs.len() != b_seqs.len() {
        return Err(PpiError::invalid(
            "paired_seqs",
            format!(
                "host has {} sequences, pathogen has {} — must pair 1:1",
                a_seqs.len(),
                b_seqs.len()
            ),
        ));
    }
    let scheme = ScoringScheme::blosum62_default();
    let a = align_msa(a_seqs, &scheme)
        .map_err(|e| PpiError::invalid("align", format!("chain A: {e}")))?;
    let b = align_msa(b_seqs, &scheme)
        .map_err(|e| PpiError::invalid("align", format!("chain B: {e}")))?;
    PairedMsa::new(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn msa(rows: &[&[u8]]) -> Msa {
        Msa::new(rows.iter().map(|r| r.to_vec()).collect()).unwrap()
    }

    #[test]
    fn score_is_always_review_flagged() {
        let a = msa(&[b"MA", b"MA", b"MT", b"MT"]);
        let b = msa(&[b"CK", b"CK", b"GK", b"GK"]);
        let paired = PairedMsa::new(a, b).unwrap();
        let s = score_pair(&paired, None).unwrap();
        assert!(s.requires_review());
        assert!(s.complementarity.is_none());
        assert!((0.0..=1.0).contains(&s.value));
    }

    #[test]
    fn coupled_pair_outscores_independent_pair() {
        // Coupled: B-col0 tracks A-col1.
        let coupled = PairedMsa::new(
            msa(&[b"MA", b"MA", b"MT", b"MT"]),
            msa(&[b"CK", b"CK", b"GK", b"GK"]),
        )
        .unwrap();
        // Independent: B varies on its own, uncorrelated with A.
        let independent = PairedMsa::new(
            msa(&[b"MA", b"MA", b"MT", b"MT"]),
            msa(&[b"CK", b"GK", b"CK", b"GK"]),
        )
        .unwrap();
        let sc = score_pair(&coupled, None).unwrap().value;
        let si = score_pair(&independent, None).unwrap().value;
        assert!(sc > si, "coupled {sc} should beat independent {si}");
    }

    #[test]
    fn rejects_bad_weight() {
        let paired =
            PairedMsa::new(msa(&[b"MA", b"MA", b"MT"]), msa(&[b"CK", b"CK", b"GK"])).unwrap();
        let w = PpiWeights {
            coevolution: -1.0,
            complementarity: 1.0,
        };
        let err = score_pair_weighted(&paired, None, &w).unwrap_err();
        assert_eq!(err.code(), "bad_weight");
    }

    #[test]
    fn screen_rejects_empty_side() {
        let host = vec![ScreenEntry::new("h", msa(&[b"AA", b"AA", b"AA"]))];
        let err = interactome_screen(&host, &[]).unwrap_err();
        assert_eq!(err.code(), "invalid");
    }
}
