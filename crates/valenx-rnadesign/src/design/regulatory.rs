//! Feature 9 — regulatory-element design (5′ / 3′ UTRs, Kozak context).
//!
//! The untranslated regions and the Kozak context around the start
//! codon decide how well a coding mRNA is translated and how long it
//! survives. This module is an orchestration layer over
//! [`valenx_genediting`]'s UTR tooling — the Kozak scorer, the
//! AU-rich-element scanner and the reference-UTR sequences are *that*
//! crate's code, never re-implemented here.
//!
//! [`design_regulatory`] selects 5′ / 3′ UTRs for a CDS (the
//! well-expressing reference UTRs, or caller-supplied ones), scores the
//! Kozak context and the 3′UTR stability, and bundles the result as a
//! [`RegulatoryDesign`].
//!
//! ## v1 scope
//!
//! UTR "design" here is **selection + scoring**: it chooses the
//! reference UTRs known to express well and reports how they score
//! against the supplied CDS. It does not synthesise a novel UTR from
//! scratch. Every score is the transparent heuristic of the
//! `valenx-genediting` UTR module (a position-weighted Kozak match, an
//! ARE count) — not a trained translation-initiation model.

use crate::error::{Result, RnaDesignError};
use serde::{Deserialize, Serialize};
use valenx_genediting::mrna::utr::{analyze_utr3, analyze_utr5, reference_utr3, reference_utr5};

/// Parameters for [`design_regulatory`].
#[derive(Clone, Debug, Default)]
pub struct RegulatoryParams {
    /// An optional caller-supplied 5′UTR (DNA or RNA). When `None`, the
    /// well-expressing reference 5′UTR is used.
    pub utr5: Option<Vec<u8>>,
    /// An optional caller-supplied 3′UTR. When `None`, the reference
    /// 3′UTR is used.
    pub utr3: Option<Vec<u8>>,
}

/// The result of a regulatory-element design (feature 9).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RegulatoryDesign {
    /// The selected 5′UTR (RNA, `A C G U`).
    pub utr5: Vec<u8>,
    /// The selected 3′UTR (RNA, `A C G U`).
    pub utr3: Vec<u8>,
    /// The Kozak-context score in `[0, 1]` for the chosen 5′UTR against
    /// the supplied CDS.
    pub kozak_score: f64,
    /// `true` when the strongest Kozak determinants are met (a purine
    /// at −3 and a `G` at +4).
    pub strong_kozak: bool,
    /// Number of upstream `AUG`s (uORF starts) in the 5′UTR.
    pub uorf_count: usize,
    /// The 3′UTR stability score in `[0, 1]` — higher = a more stable,
    /// ARE-free trailer.
    pub utr3_stability: f64,
    /// Number of AU-rich / destabilising elements found in the 3′UTR.
    pub utr3_destabilizers: usize,
    /// Human-readable notes describing the regulatory design.
    pub notes: Vec<String>,
}

impl RegulatoryDesign {
    /// `true` when the regulatory design is solid: a strong Kozak, no
    /// uORFs, and a clean (ARE-free) 3′UTR.
    pub fn is_high_quality(&self) -> bool {
        self.strong_kozak && self.uorf_count == 0 && self.utr3_destabilizers == 0
    }
}

/// Designs the regulatory elements for a coding mRNA (feature 9).
///
/// `cds` is the coding sequence whose start codon the Kozak context
/// wraps (DNA or RNA). The 5′ / 3′ UTRs come from `params` if supplied,
/// otherwise the well-expressing reference UTRs are used; both are
/// scored against the CDS.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if the CDS is shorter than 4 nt (too
///   short to read the Kozak `+4` position).
/// - [`RnaDesignError::Upstream`] if the UTR analysers reject the
///   input.
pub fn design_regulatory(cds: &[u8], params: &RegulatoryParams) -> Result<RegulatoryDesign> {
    if cds.len() < 4 {
        return Err(RnaDesignError::goal(
            "protein",
            "the CDS is too short to evaluate the Kozak context",
        ));
    }

    let utr5_supplied = params.utr5.is_some();
    let utr3_supplied = params.utr3.is_some();
    let utr5 = params.utr5.clone().unwrap_or_else(reference_utr5);
    let utr3 = params.utr3.clone().unwrap_or_else(reference_utr3);

    // Score both UTRs against the CDS via the genediting UTR module.
    let a5 = analyze_utr5(&utr5, cds)?;
    let a3 = analyze_utr3(&utr3)?;

    let mut notes = Vec::new();
    notes.push(format!(
        "5'UTR: {} ({} nt) — Kozak score {:.2} ({}), {} uORF(s).",
        if utr5_supplied {
            "caller-supplied"
        } else {
            "well-expressing reference"
        },
        a5.length,
        a5.kozak_score,
        if a5.strong_kozak { "strong" } else { "weak" },
        a5.uorf_count,
    ));
    notes.push(format!(
        "3'UTR: {} ({} nt) — stability {:.2}, {} destabilising element(s).",
        if utr3_supplied {
            "caller-supplied"
        } else {
            "ARE-free reference"
        },
        a3.length,
        a3.stability_score,
        a3.destabilizer_count(),
    ));
    if a5.uorf_count > 0 {
        notes.push(
            "The 5'UTR carries upstream AUG(s) — each can sequester scanning \
             ribosomes and depress main-ORF translation."
                .to_string(),
        );
    }
    if !a3.is_clean() {
        notes.push(
            "The 3'UTR carries an AU-rich element — it recruits the mRNA-decay \
             machinery and shortens half-life; prefer an ARE-free 3'UTR."
                .to_string(),
        );
    }
    notes.push(
        "UTR scoring uses transparent heuristics (a position-weighted Kozak match, \
         an ARE count) — not a trained translation model."
            .to_string(),
    );

    Ok(RegulatoryDesign {
        utr5,
        utr3,
        kozak_score: a5.kozak_score,
        strong_kozak: a5.strong_kozak,
        uorf_count: a5.uorf_count,
        utr3_stability: a3.stability_score,
        utr3_destabilizers: a3.destabilizer_count(),
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_reference_utrs() {
        let d = design_regulatory(b"ATGGCCGCCTAA", &RegulatoryParams::default()).unwrap();
        assert!(!d.utr5.is_empty());
        assert!(!d.utr3.is_empty());
        assert!((0.0..=1.0).contains(&d.kozak_score));
        assert!((0.0..=1.0).contains(&d.utr3_stability));
        assert!(!d.notes.is_empty());
    }

    #[test]
    fn reference_utrs_are_decent() {
        // The built-in reference UTRs should score reasonably.
        let d = design_regulatory(b"ATGGCCGCCGCCTAA", &RegulatoryParams::default()).unwrap();
        // The reference 3'UTR is ARE-free.
        assert_eq!(d.utr3_destabilizers, 0);
    }

    #[test]
    fn accepts_supplied_utrs() {
        let params = RegulatoryParams {
            utr5: Some(b"GGGAAAGCCACC".to_vec()),
            utr3: Some(b"GCGCGCGCGCGC".to_vec()),
        };
        let d = design_regulatory(b"ATGGCCTAA", &params).unwrap();
        assert!(d.notes.iter().any(|n| n.contains("caller-supplied")));
    }

    #[test]
    fn rejects_too_short_cds() {
        assert!(design_regulatory(b"AT", &RegulatoryParams::default()).is_err());
    }

    #[test]
    fn high_quality_flag() {
        let d = design_regulatory(b"ATGGCCGCCGCCTAA", &RegulatoryParams::default()).unwrap();
        assert_eq!(
            d.is_high_quality(),
            d.strong_kozak && d.uorf_count == 0 && d.utr3_destabilizers == 0
        );
    }

    #[test]
    fn flags_uorf_in_5utr() {
        // A 5'UTR with an internal AUG.
        let params = RegulatoryParams {
            utr5: Some(b"GGGAUGAAAGCCACC".to_vec()),
            utr3: None,
        };
        let d = design_regulatory(b"ATGGCCTAA", &params).unwrap();
        assert!(d.uorf_count >= 1);
    }
}
