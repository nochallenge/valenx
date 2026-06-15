//! # valenx-lrp5
//!
//! A **second-target cross-reactivity demonstration** for the biodesign pipeline,
//! mirroring the myostatin demo on a different protein family. The therapeutic
//! goal is an **LRP5-specific binder** (LRP5 is a Wnt co-receptor; gain-of-function
//! variants cause high bone mass, making it a bone-density drug target). The hard
//! part is *specificity*: LRP5's closest paralog **LRP6** is a near-twin Wnt
//! co-receptor, so a binder that also hits LRP6 would perturb a different pathway —
//! the textbook off-target hazard, analogous to GDF-8/GDF-11 in myostatin.
//!
//! This crate screens LRP5 against a panel — the paralogs LRP6 and LRP4 plus an
//! unrelated negative control (GDF-8 / myostatin) — using the two similarity
//! measures in [`valenx_offtarget`], computed from **verified UniProt sequences**
//! (see [`sequences`]).
//!
//! ## The metric lesson
//!
//! [`valenx_offtarget::best_ungapped_identity`] is the right tool for *short,
//! alignable* mature domains (it gave GDF-8/GDF-11 = 89.9% in the myostatin demo).
//! For *full-length, multidomain* paralogs like LRP5/LRP6 — which are ~71%
//! identical but riddled with insertions/deletions — an ungapped slide finds only
//! a few percent and badly **understates** the homology. The alignment-free
//! [`valenx_offtarget::kmer_jaccard`] is robust to those indels and correctly
//! ranks the paralog far above the unrelated control. This crate reports both, so
//! the contrast is explicit rather than hidden.
//!
//! ## Honesty
//!
//! Research/educational. The sequences and both computed numbers are real and
//! verifiable. k-mer Jaccard is a transparent similarity heuristic, **not** a
//! clinical cross-reactivity prediction — a real specificity call needs structural
//! modelling and wet-lab binding assays. Nothing here clears a candidate as safe.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod sequences;

use valenx_offtarget::{best_ungapped_identity, kmer_jaccard};

use sequences::{GDF8_O14793, LRP4_O75096, LRP5_O75197, LRP6_O75581};

/// k-mer length for the alignment-free Jaccard metric (penta-peptides: long
/// enough that unrelated proteins share few, short enough for sensitivity).
pub const KMER_K: usize = 5;

/// LRP5's computed similarity to one reference protein.
#[derive(Debug, Clone, PartialEq)]
pub struct CrossReactivity {
    /// The reference protein (name + UniProt accession).
    pub reference: String,
    /// Best ungapped fractional identity in `[0, 1]` — conservative; misses
    /// indels, so it understates full-length paralog homology.
    pub ungapped_identity: f64,
    /// Alignment-free k-mer Jaccard overlap in `[0, 1]` — robust to indels; the
    /// metric to trust for full-length multidomain proteins.
    pub kmer_jaccard: f64,
}

/// Screen LRP5 against its paralogs (LRP6, LRP4) and an unrelated negative
/// control (GDF-8), computing both similarity metrics from the verified
/// sequences. Returned in panel order (LRP6, LRP4, GDF-8).
pub fn screen_lrp5_off_targets() -> Vec<CrossReactivity> {
    let panel = [
        ("LRP6 (O75581, closest paralog)", LRP6_O75581),
        ("LRP4 (O75096, paralog)", LRP4_O75096),
        ("GDF8/myostatin (O14793, unrelated control)", GDF8_O14793),
    ];
    panel
        .iter()
        .map(|(name, seq)| CrossReactivity {
            reference: (*name).to_string(),
            ungapped_identity: best_ungapped_identity(LRP5_O75197, seq)
                .expect("verified sequences are valid"),
            kmer_jaccard: kmer_jaccard(LRP5_O75197, seq, KMER_K)
                .expect("verified sequences are valid"),
        })
        .collect()
}

/// The panel entry with the highest k-mer Jaccard — the worst predicted
/// off-target liability for an LRP5 binder.
pub fn worst_off_target() -> CrossReactivity {
    screen_lrp5_off_targets()
        .into_iter()
        .max_by(|a, b| {
            a.kmer_jaccard
                .partial_cmp(&b.kmer_jaccard)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("panel is non-empty")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequences_have_verified_lengths() {
        assert_eq!(LRP5_O75197.len(), 1615);
        assert_eq!(LRP6_O75581.len(), 1613);
        assert_eq!(LRP4_O75096.len(), 1905);
        assert_eq!(GDF8_O14793.len(), 375);
    }

    #[test]
    fn jaccard_ranks_paralog_above_unrelated_control() {
        let hits = screen_lrp5_off_targets();
        assert_eq!(hits.len(), 3);
        let lrp6 = &hits[0];
        let lrp4 = &hits[1];
        let gdf8 = &hits[2];

        // All metrics in range.
        for h in &hits {
            assert!((0.0..=1.0).contains(&h.kmer_jaccard));
            assert!((0.0..=1.0).contains(&h.ungapped_identity));
        }

        // The closest paralog (LRP6) is the worst off-target by k-mer Jaccard,
        // far above the unrelated GDF-8 control; LRP4 (a more distant paralog)
        // sits in between or above the control.
        assert!(
            lrp6.kmer_jaccard > gdf8.kmer_jaccard,
            "LRP6 ({:.4}) should exceed GDF8 control ({:.4})",
            lrp6.kmer_jaccard,
            gdf8.kmer_jaccard
        );
        assert!(
            lrp6.kmer_jaccard >= lrp4.kmer_jaccard,
            "LRP6 ({:.4}) should be the closest paralog, >= LRP4 ({:.4})",
            lrp6.kmer_jaccard,
            lrp4.kmer_jaccard
        );
        assert!(
            lrp4.kmer_jaccard > gdf8.kmer_jaccard,
            "even the more distant paralog LRP4 ({:.4}) should beat the unrelated control ({:.4})",
            lrp4.kmer_jaccard,
            gdf8.kmer_jaccard
        );
    }

    #[test]
    fn worst_off_target_is_the_paralog_lrp6() {
        let worst = worst_off_target();
        assert!(worst.reference.contains("LRP6"));
    }

    #[test]
    fn ungapped_identity_understates_fulllength_homology() {
        // The metric lesson, asserted: ungapped identity is far below the k-mer
        // Jaccard for the LRP5/LRP6 paralog pair (it cannot bridge the indels).
        let lrp6 = &screen_lrp5_off_targets()[0];
        assert!(
            lrp6.ungapped_identity < lrp6.kmer_jaccard,
            "ungapped identity ({:.4}) should understate the paralog homology that \
             k-mer Jaccard ({:.4}) captures",
            lrp6.ungapped_identity,
            lrp6.kmer_jaccard
        );
    }
}
