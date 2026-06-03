//! Feature 22 — poly-A tail and 5′ cap-analog selection.
//!
//! Two end modifications set how long a therapeutic mRNA survives and
//! how well it is translated:
//!
//! - the **poly-A tail** — its length tunes half-life. Too short and
//!   the mRNA is degraded quickly; ~100–150 nt is the therapeutic
//!   sweet spot; beyond ~150 nt there is little extra benefit and
//!   plasmid-encoded long tails are unstable. A short internal linker
//!   ("segmented" poly-A) can stabilise the template.
//! - the **5′ cap** — the cap chemistry decides ribosome recruitment
//!   and innate-immune visibility. A 2′-O-methylated **cap 1** evades
//!   the RIG-I / IFIT sensors; a co-transcriptional **CleanCap**
//!   gives near-100 % capping in one IVT step.
//!
//! This module recommends a poly-A length for a use case
//! ([`recommend_poly_a`]) and a cap analog ([`recommend_cap`]), and
//! scores a chosen pair ([`score_ends`]).
//!
//! ## v1 scope
//!
//! The recommendations encode the well-established qualitative rules
//! (the length / half-life relationship, cap-1 immune evasion); they
//! are not a quantitative pharmacokinetic model.

use crate::error::{GeneditingError, Result};
use crate::mrna::construct::CapType;
use serde::{Deserialize, Serialize};

/// The intended use of an mRNA construct — drives the poly-A / cap
/// recommendation.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MrnaUseCase {
    /// A prophylactic or therapeutic **vaccine** — a transient,
    /// strongly-translated antigen; cap-1 chemistry is essential.
    Vaccine,
    /// **Protein-replacement therapy** — sustained expression of a
    /// missing protein; favours the longer end of the tail range.
    ProteinReplacement,
    /// **Transient reprogramming / gene editing** delivery (e.g. mRNA
    /// encoding a nuclease) — a short, sharp pulse of expression.
    TransientEditing,
}

impl MrnaUseCase {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            MrnaUseCase::Vaccine => "vaccine",
            MrnaUseCase::ProteinReplacement => "protein replacement",
            MrnaUseCase::TransientEditing => "transient editing delivery",
        }
    }
}

/// A poly-A tail recommendation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PolyARecommendation {
    /// Recommended poly-A tail length in nucleotides.
    pub length: usize,
    /// `true` when a segmented (linker-interrupted) tail is suggested
    /// to keep the encoding template stable.
    pub segmented: bool,
    /// A one-line rationale.
    pub rationale: String,
}

/// Recommends a poly-A tail length for a use case (feature 22).
///
/// Vaccines and editing-delivery mRNA want a moderate ~100–120 nt
/// tail; protein-replacement mRNA benefits from the longer ~150 nt end
/// of the useful range. A segmented tail is suggested for the longer
/// tails, where a plasmid-encoded continuous A-run is hard to maintain.
pub fn recommend_poly_a(use_case: MrnaUseCase) -> PolyARecommendation {
    match use_case {
        MrnaUseCase::Vaccine => PolyARecommendation {
            length: 100,
            segmented: false,
            rationale: "~100 nt — strong translation with a transient \
                        half-life, the established vaccine-mRNA tail."
                .to_string(),
        },
        MrnaUseCase::ProteinReplacement => PolyARecommendation {
            length: 150,
            segmented: true,
            rationale: "~150 nt — the long end of the useful range for \
                        sustained expression; a segmented tail keeps \
                        the encoding template stable."
                .to_string(),
        },
        MrnaUseCase::TransientEditing => PolyARecommendation {
            length: 100,
            segmented: false,
            rationale: "~100 nt — a sharp expression pulse is enough \
                        for a nuclease / editor mRNA; a longer tail \
                        only prolongs an editing window you want short."
                .to_string(),
        },
    }
}

/// A 5′ cap recommendation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CapRecommendation {
    /// Recommended cap chemistry.
    pub cap: CapType,
    /// A one-line rationale.
    pub rationale: String,
}

/// Recommends a 5′ cap analog for a use case (feature 22).
///
/// Every therapeutic use case here wants a 2′-O-methylated cap-1
/// chemistry for innate-immune evasion; [`MrnaUseCase::Vaccine`] and
/// [`MrnaUseCase::TransientEditing`] are matched to a co-transcriptional
/// CleanCap (one-step IVT, near-complete capping).
pub fn recommend_cap(use_case: MrnaUseCase) -> CapRecommendation {
    match use_case {
        MrnaUseCase::Vaccine | MrnaUseCase::TransientEditing => CapRecommendation {
            cap: CapType::CleanCap,
            rationale: "CleanCap — co-transcriptional cap 1; near-100 % \
                        capping in a single IVT step and full RIG-I / \
                        IFIT immune evasion."
                .to_string(),
        },
        MrnaUseCase::ProteinReplacement => CapRecommendation {
            cap: CapType::Cap1,
            rationale: "Cap 1 — 2'-O-methylated cap for innate-immune \
                        evasion; enzymatic capping gives reliable cap-1 \
                        for a sustained-expression construct."
                .to_string(),
        },
    }
}

/// A score for a chosen poly-A / cap pair.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EndScore {
    /// Poly-A length component of the score in `[0, 1]`.
    pub poly_a_score: f64,
    /// Cap component of the score in `[0, 1]`.
    pub cap_score: f64,
    /// Combined end-modification score in `[0, 1]`.
    pub combined: f64,
}

/// Scores a chosen poly-A length and cap chemistry (feature 22).
///
/// The poly-A score peaks across the ~100–150 nt band and falls off
/// outside it; the cap score rewards an innate-immune-silent cap-1
/// chemistry.
///
/// # Errors
/// [`GeneditingError::Invalid`] for an absurd poly-A length
/// (`> 500` nt).
pub fn score_ends(poly_a_len: usize, cap: CapType) -> Result<EndScore> {
    if poly_a_len > 500 {
        return Err(GeneditingError::invalid(
            "poly_a_len",
            "poly-A tail length is beyond any physiological value",
        ));
    }
    // Poly-A: a trapezoid — flat ~1.0 across 100..=150, falling outside.
    let l = poly_a_len as f64;
    let poly_a_score = if (100..=150).contains(&poly_a_len) {
        1.0
    } else if poly_a_len < 100 {
        (l / 100.0).clamp(0.0, 1.0)
    } else {
        // Slow decay above 150 (extra length adds nothing, long
        // templates are unstable).
        (1.0 - (l - 150.0) / 350.0).clamp(0.0, 1.0)
    };
    // Cap: cap-1 chemistries score full, cap-0 / ARCA lower.
    let cap_score = if cap.is_innate_immune_silent() {
        1.0
    } else {
        0.55
    };
    let combined = (0.5 * poly_a_score + 0.5 * cap_score).clamp(0.0, 1.0);
    Ok(EndScore {
        poly_a_score,
        cap_score,
        combined,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vaccine_gets_a_moderate_tail() {
        let r = recommend_poly_a(MrnaUseCase::Vaccine);
        assert_eq!(r.length, 100);
        assert!(!r.segmented);
    }

    #[test]
    fn protein_replacement_gets_a_longer_segmented_tail() {
        let r = recommend_poly_a(MrnaUseCase::ProteinReplacement);
        assert!(r.length >= 140);
        assert!(r.segmented);
    }

    #[test]
    fn all_use_cases_get_a_cap1_chemistry() {
        for uc in [
            MrnaUseCase::Vaccine,
            MrnaUseCase::ProteinReplacement,
            MrnaUseCase::TransientEditing,
        ] {
            let c = recommend_cap(uc);
            assert!(c.cap.is_innate_immune_silent(), "{} needs cap 1", uc.name());
        }
    }

    #[test]
    fn vaccine_cap_is_cleancap() {
        assert_eq!(recommend_cap(MrnaUseCase::Vaccine).cap, CapType::CleanCap);
    }

    #[test]
    fn end_score_peaks_in_band() {
        let in_band = score_ends(120, CapType::Cap1).unwrap();
        let too_short = score_ends(40, CapType::Cap1).unwrap();
        let too_long = score_ends(400, CapType::Cap1).unwrap();
        assert!((in_band.poly_a_score - 1.0).abs() < 1e-9);
        assert!(too_short.poly_a_score < in_band.poly_a_score);
        assert!(too_long.poly_a_score < in_band.poly_a_score);
    }

    #[test]
    fn cap1_outscores_cap0() {
        let c1 = score_ends(120, CapType::Cap1).unwrap();
        let c0 = score_ends(120, CapType::Cap0).unwrap();
        assert!(c1.cap_score > c0.cap_score);
        assert!(c1.combined > c0.combined);
    }

    #[test]
    fn score_in_unit_range() {
        for len in [0usize, 50, 120, 200, 500] {
            let s = score_ends(len, CapType::CleanCap).unwrap();
            assert!((0.0..=1.0).contains(&s.combined));
        }
    }

    #[test]
    fn rejects_absurd_poly_a() {
        assert!(score_ends(9999, CapType::Cap1).is_err());
    }
}
