//! Feature 27 — the edit-strategy advisor.
//!
//! Given a desired genomic change, which editing modality should you
//! use — a CRISPR nuclease cut, a base editor, or a prime editor? Each
//! has a sweet spot:
//!
//! - **Base editing** — the cleanest route for a single
//!   `C•G↔T•A` / `A•T↔G•C` transition, *if* a PAM places the target
//!   base in an editor's window. No double-strand break, no indels.
//! - **Prime editing** — handles *any* of the 12 substitutions plus
//!   small insertions / deletions, no donor; lower efficiency, more
//!   design work.
//! - **Nuclease + HDR** — the route for a large knock-in (a tag, a
//!   cassette); needs a donor and a double-strand break.
//! - **Nuclease + NHEJ** — the route for a gene **knockout** (the
//!   desired change *is* loss of function).
//!
//! This module classifies a [`DesiredChange`], scores each approach
//! with a transparent rubric, and returns a ranked
//! [`StrategyAdvice`].
//!
//! ## v1 scope
//!
//! The advisor reasons from the *kind and size* of the change and the
//! transition class — a transparent decision rubric, not a trained
//! recommender. It does not itself scan a locus for PAMs; the workflow
//! [`driver`](crate::workflow::driver) does that and feeds the
//! advisor a `pam_available` hint.

use serde::{Deserialize, Serialize};

/// The kind of genomic change a user wants to make.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum DesiredChange {
    /// A single-nucleotide substitution `from → to`.
    PointMutation {
        /// The reference base (`A`/`C`/`G`/`T`).
        from: u8,
        /// The desired base.
        to: u8,
    },
    /// A small insertion (`size` bp).
    SmallInsertion {
        /// Insertion size in base pairs.
        size: usize,
    },
    /// A small deletion (`size` bp).
    SmallDeletion {
        /// Deletion size in base pairs.
        size: usize,
    },
    /// A large knock-in — a tag, a reporter, a cassette (`size` bp).
    LargeKnockIn {
        /// Knock-in size in base pairs.
        size: usize,
    },
    /// A gene knockout — the desired outcome is loss of function.
    GeneKnockout,
}

impl DesiredChange {
    /// `true` when the change is a base-editor-reachable transition
    /// (`C↔T` or `A↔G`, either strand).
    pub fn is_base_editable_transition(&self) -> bool {
        match self {
            DesiredChange::PointMutation { from, to } => {
                let f = from.to_ascii_uppercase();
                let t = to.to_ascii_uppercase();
                matches!(
                    (f, t),
                    (b'C', b'T')
                        | (b'T', b'C')
                        | (b'A', b'G')
                        | (b'G', b'A')
                )
            }
            _ => false,
        }
    }

    /// A short human-readable label.
    pub fn label(&self) -> String {
        match self {
            DesiredChange::PointMutation { from, to } => {
                format!("{}->{} point mutation", *from as char, *to as char)
            }
            DesiredChange::SmallInsertion { size } => format!("{size}-bp insertion"),
            DesiredChange::SmallDeletion { size } => format!("{size}-bp deletion"),
            DesiredChange::LargeKnockIn { size } => format!("{size}-bp knock-in"),
            DesiredChange::GeneKnockout => "gene knockout".to_string(),
        }
    }
}

/// An editing modality the advisor can recommend.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EditApproach {
    /// A base editor (CBE / ABE).
    BaseEditing,
    /// A prime editor (PE2 / PE3 / PE3b / PEmax).
    PrimeEditing,
    /// A CRISPR nuclease cut repaired by HDR with a donor template.
    NucleaseHdr,
    /// A CRISPR nuclease cut repaired by NHEJ (knockout).
    NucleaseNhej,
}

impl EditApproach {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            EditApproach::BaseEditing => "base editing",
            EditApproach::PrimeEditing => "prime editing",
            EditApproach::NucleaseHdr => "nuclease + HDR knock-in",
            EditApproach::NucleaseNhej => "nuclease + NHEJ knockout",
        }
    }

    /// `true` when the approach makes a double-strand break.
    pub fn makes_dsb(self) -> bool {
        matches!(self, EditApproach::NucleaseHdr | EditApproach::NucleaseNhej)
    }
}

/// One scored approach in the advisor's ranking.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ScoredApproach {
    /// The editing modality.
    pub approach: EditApproach,
    /// A suitability score in `[0, 1]` — higher is a better fit.
    pub suitability: f64,
    /// A one-line rationale for this score.
    pub rationale: String,
}

/// The advisor's ranked recommendation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StrategyAdvice {
    /// The change that was analysed.
    pub change: DesiredChange,
    /// Every approach, scored and sorted best-first.
    pub ranked: Vec<ScoredApproach>,
    /// A one-line summary naming the recommended approach.
    pub summary: String,
}

impl StrategyAdvice {
    /// The single recommended approach (the top of the ranking).
    pub fn recommended(&self) -> EditApproach {
        self.ranked
            .first()
            .map(|s| s.approach)
            .unwrap_or(EditApproach::NucleaseNhej)
    }
}

/// Advises on an editing strategy for a desired change (feature 27).
///
/// Scores each modality with a transparent rubric keyed on the kind /
/// size of the change and the transition class. `pam_available` is a
/// hint (from a prior locus scan) that a PAM-adjacent guide exists —
/// when `false`, the PAM-dependent approaches are demoted.
///
/// The result is always a *complete* ranking: every approach is
/// scored, so the caller can see the trade-offs, not just the winner.
pub fn advise_strategy(change: &DesiredChange, pam_available: bool) -> StrategyAdvice {
    let mut ranked: Vec<ScoredApproach> = vec![
        score_base_editing(change, pam_available),
        score_prime_editing(change, pam_available),
        score_nuclease_hdr(change, pam_available),
        score_nuclease_nhej(change, pam_available),
    ];
    ranked.sort_by(|a, b| {
        b.suitability
            .partial_cmp(&a.suitability)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let best = &ranked[0];
    let summary = format!(
        "For a {}, the recommended approach is {} (suitability {:.2}). {}",
        change.label(),
        best.approach.name(),
        best.suitability,
        best.rationale,
    );
    StrategyAdvice {
        change: change.clone(),
        ranked,
        summary,
    }
}

/// Scores base editing for a change.
fn score_base_editing(change: &DesiredChange, pam: bool) -> ScoredApproach {
    let (score, why) = match change {
        DesiredChange::PointMutation { .. } if change.is_base_editable_transition() => {
            if pam {
                (
                    0.95,
                    "a base-editor-reachable transition with a PAM in range — \
                     the cleanest, indel-free route",
                )
            } else {
                (
                    0.45,
                    "a base-editable transition, but no PAM places the base in \
                     an editor window here",
                )
            }
        }
        DesiredChange::PointMutation { .. } => (
            0.10,
            "this substitution is not a C↔T / A↔G transition — base editors \
             cannot make it",
        ),
        _ => (
            0.05,
            "base editors install point transitions only — not insertions, \
             deletions, knock-ins or knockouts",
        ),
    };
    ScoredApproach {
        approach: EditApproach::BaseEditing,
        suitability: score,
        rationale: why.to_string(),
    }
}

/// Scores prime editing for a change.
fn score_prime_editing(change: &DesiredChange, pam: bool) -> ScoredApproach {
    let (score, why) = match change {
        DesiredChange::PointMutation { .. } => {
            if pam {
                (
                    0.80,
                    "prime editing installs any substitution without a \
                     double-strand break — a strong, versatile choice",
                )
            } else {
                (0.40, "prime editing fits, but needs a PAM-adjacent spacer")
            }
        }
        DesiredChange::SmallInsertion { size } | DesiredChange::SmallDeletion { size } => {
            if *size <= 44 && pam {
                (
                    0.85,
                    "prime editing is the leading route for a small indel — \
                     no donor, no double-strand break",
                )
            } else if *size <= 44 {
                (0.45, "fits a small indel, but needs a PAM-adjacent spacer")
            } else {
                (
                    0.30,
                    "the indel is larger than a pegRNA RT template comfortably \
                     encodes — twin-pegRNA or HDR may be better",
                )
            }
        }
        DesiredChange::LargeKnockIn { .. } => (
            0.20,
            "a large knock-in exceeds a standard pegRNA's capacity — HDR is \
             the usual route (twin-pegRNA can extend prime editing)",
        ),
        DesiredChange::GeneKnockout => (
            0.35,
            "prime editing can install a premature stop, but an NHEJ \
             knockout is simpler when the goal is loss of function",
        ),
    };
    ScoredApproach {
        approach: EditApproach::PrimeEditing,
        suitability: score,
        rationale: why.to_string(),
    }
}

/// Scores nuclease + HDR for a change.
fn score_nuclease_hdr(change: &DesiredChange, pam: bool) -> ScoredApproach {
    let base = match change {
        DesiredChange::LargeKnockIn { .. } => 0.90,
        DesiredChange::PointMutation { .. } => 0.55,
        DesiredChange::SmallInsertion { .. } | DesiredChange::SmallDeletion { .. } => 0.55,
        DesiredChange::GeneKnockout => 0.20,
    };
    let score = if pam { base } else { base * 0.5 };
    let why = match change {
        DesiredChange::LargeKnockIn { .. } => {
            "HDR with a donor template is the established route for a large \
             knock-in"
        }
        DesiredChange::GeneKnockout => {
            "HDR can disrupt a gene, but NHEJ is simpler for a pure knockout"
        }
        _ => {
            "HDR can make this change precisely, but it needs a donor and a \
             double-strand break and is inefficient in non-dividing cells"
        }
    };
    ScoredApproach {
        approach: EditApproach::NucleaseHdr,
        suitability: score,
        rationale: why.to_string(),
    }
}

/// Scores nuclease + NHEJ for a change.
fn score_nuclease_nhej(change: &DesiredChange, pam: bool) -> ScoredApproach {
    let (base, why) = match change {
        DesiredChange::GeneKnockout => (
            0.92,
            "an NHEJ frameshift is the standard, efficient route to a gene \
             knockout",
        ),
        DesiredChange::SmallDeletion { .. } => (
            0.55,
            "paired cuts give a defined NHEJ deletion, though the junction \
             is imprecise",
        ),
        _ => (
            0.10,
            "NHEJ produces uncontrolled indels — unsuitable for a precise \
             change",
        ),
    };
    let score = if pam { base } else { base * 0.5 };
    ScoredApproach {
        approach: EditApproach::NucleaseNhej,
        suitability: score,
        rationale: why.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transition_classifier() {
        assert!(DesiredChange::PointMutation { from: b'C', to: b'T' }
            .is_base_editable_transition());
        assert!(DesiredChange::PointMutation { from: b'A', to: b'G' }
            .is_base_editable_transition());
        // A transversion is not base-editable.
        assert!(!DesiredChange::PointMutation { from: b'C', to: b'G' }
            .is_base_editable_transition());
        assert!(!DesiredChange::GeneKnockout.is_base_editable_transition());
    }

    #[test]
    fn base_editing_wins_for_a_transition_with_pam() {
        let change = DesiredChange::PointMutation { from: b'A', to: b'G' };
        let advice = advise_strategy(&change, true);
        assert_eq!(advice.recommended(), EditApproach::BaseEditing);
    }

    #[test]
    fn prime_editing_wins_for_a_transversion() {
        // C->G is not base-editable; prime editing should top the rank.
        let change = DesiredChange::PointMutation { from: b'C', to: b'G' };
        let advice = advise_strategy(&change, true);
        assert_eq!(advice.recommended(), EditApproach::PrimeEditing);
    }

    #[test]
    fn nhej_wins_for_a_knockout() {
        let advice = advise_strategy(&DesiredChange::GeneKnockout, true);
        assert_eq!(advice.recommended(), EditApproach::NucleaseNhej);
    }

    #[test]
    fn hdr_wins_for_a_large_knock_in() {
        let change = DesiredChange::LargeKnockIn { size: 800 };
        let advice = advise_strategy(&change, true);
        assert_eq!(advice.recommended(), EditApproach::NucleaseHdr);
    }

    #[test]
    fn prime_editing_wins_for_a_small_insertion() {
        let change = DesiredChange::SmallInsertion { size: 12 };
        let advice = advise_strategy(&change, true);
        assert_eq!(advice.recommended(), EditApproach::PrimeEditing);
    }

    #[test]
    fn no_pam_demotes_pam_dependent_approaches() {
        let change = DesiredChange::PointMutation { from: b'A', to: b'G' };
        let with_pam = advise_strategy(&change, true);
        let without = advise_strategy(&change, false);
        // Base editing's score must drop when no PAM is available.
        let be_with = with_pam
            .ranked
            .iter()
            .find(|s| s.approach == EditApproach::BaseEditing)
            .unwrap()
            .suitability;
        let be_without = without
            .ranked
            .iter()
            .find(|s| s.approach == EditApproach::BaseEditing)
            .unwrap()
            .suitability;
        assert!(be_without < be_with);
    }

    #[test]
    fn ranking_is_complete_and_sorted() {
        let advice = advise_strategy(&DesiredChange::GeneKnockout, true);
        assert_eq!(advice.ranked.len(), 4); // all four approaches
        for w in advice.ranked.windows(2) {
            assert!(w[0].suitability >= w[1].suitability);
        }
        for s in &advice.ranked {
            assert!((0.0..=1.0).contains(&s.suitability));
        }
    }

    #[test]
    fn summary_names_the_recommendation() {
        let advice = advise_strategy(&DesiredChange::GeneKnockout, true);
        assert!(advice.summary.contains("knockout"));
        assert!(advice.summary.contains(advice.recommended().name()));
    }

    #[test]
    fn change_labels() {
        assert!(DesiredChange::GeneKnockout.label().contains("knockout"));
        assert!(DesiredChange::LargeKnockIn { size: 500 }
            .label()
            .contains("500"));
    }
}
