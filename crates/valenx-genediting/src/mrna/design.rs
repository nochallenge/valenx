//! The top-level mRNA-design driver and [`MrnaDesignReport`].
//!
//! This module is the mRNA half of feature 30: it runs the
//! mRNA-therapeutic features as one workflow and bundles the result.
//! Given a CDS and a use case it:
//!
//! 1. codon-optimises the CDS for the host ([`crate::mrna::codon`]);
//! 2. optionally depletes uridine for an `m1Ψ` construct
//!    ([`crate::mrna::uridine`]);
//! 3. picks reference UTRs and checks them ([`crate::mrna::utr`]);
//! 4. minimises start-codon-region structure
//!    ([`crate::mrna::structure`]);
//! 5. recommends a poly-A length and a cap ([`crate::mrna::tailcap`]);
//! 6. assembles a validated [`MrnaConstruct`]
//!    ([`crate::mrna::construct`]);
//! 7. plans the modified-nucleoside substitution.
//!
//! and returns an [`MrnaDesignReport`] with the construct plus every
//! intermediate score.
//!
//! ## v1 scope
//!
//! The driver chains the per-feature passes in a fixed order; it does
//! not jointly optimise (a uridine pass can slightly perturb the CAI
//! the codon pass set, etc.) — each pass's effect is reported so the
//! caller sees the trade-offs. Every score is the transparent
//! heuristic of its source module.

use crate::error::Result;
use crate::mrna::codon::{optimize_cds, ExpressionHost};
use crate::mrna::construct::{CapType, MrnaConstruct, MrnaConstructBuilder};
use crate::mrna::structure::minimize_start_structure;
use crate::mrna::tailcap::{recommend_cap, recommend_poly_a, MrnaUseCase};
use crate::mrna::uridine::{
    minimize_uridine, plan_modification, ModificationPlan, ModifiedNucleoside,
};
use crate::mrna::utr::{analyze_utr3, analyze_utr5, reference_utr3, reference_utr5};
use serde::{Deserialize, Serialize};

/// A request for an end-to-end mRNA-therapeutic design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MrnaDesignRequest {
    /// The coding sequence to express (DNA or RNA; `AUG`…stop, length
    /// divisible by three).
    pub cds: Vec<u8>,
    /// The expression host.
    pub host: ExpressionHost,
    /// The intended use — drives the poly-A / cap recommendation.
    pub use_case: MrnaUseCase,
    /// `true` to run a uridine-depleting pass (recommended for an
    /// `m1Ψ` therapeutic construct).
    pub deplete_uridine: bool,
    /// The modified nucleoside the construct will be made with.
    pub nucleoside: ModifiedNucleoside,
    /// How many CDS codons the structure-minimisation pass may scan
    /// near the start codon.
    pub structure_scan_codons: usize,
}

impl MrnaDesignRequest {
    /// A request with therapeutic defaults: human host, `m1Ψ`, uridine
    /// depletion on, 6-codon structure scan.
    pub fn new(cds: impl Into<Vec<u8>>, use_case: MrnaUseCase) -> Self {
        MrnaDesignRequest {
            cds: cds.into(),
            host: ExpressionHost::Human,
            use_case,
            deplete_uridine: true,
            nucleoside: ModifiedNucleoside::N1MethylPseudouridine,
            structure_scan_codons: 6,
        }
    }
}

/// The bundled result of an end-to-end mRNA-therapeutic design
/// (feature 30, mRNA half).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MrnaDesignReport {
    /// The assembled, validated mRNA construct.
    pub construct: MrnaConstruct,
    /// CAI of the CDS *before* optimisation.
    pub cai_before: f64,
    /// CAI of the final (optimised, possibly uridine-depleted) CDS.
    pub cai_after: f64,
    /// Uridine fraction of the final CDS.
    pub uridine_fraction: f64,
    /// Start-region structural openness in `[0, 1]` of the final
    /// construct (`1.0` = a fully open start codon).
    pub start_openness: f64,
    /// The 5′UTR Kozak-context score in `[0, 1]`.
    pub kozak_score: f64,
    /// The 3′UTR stability score in `[0, 1]`.
    pub utr3_stability: f64,
    /// The recommended poly-A tail length.
    pub poly_a_len: usize,
    /// The recommended cap chemistry.
    pub cap: CapType,
    /// The modified-nucleoside substitution plan.
    pub modification: ModificationPlan,
    /// An overall design score in `[0, 1]` — a transparent blend of
    /// the per-feature scores.
    pub overall_score: f64,
    /// Human-readable design notes (one line per major decision).
    pub notes: Vec<String>,
}

impl MrnaDesignReport {
    /// `true` when the design is solid on every axis (good CAI, an
    /// open start codon, a stable 3′UTR, a strong Kozak).
    pub fn is_high_quality(&self) -> bool {
        self.overall_score >= 0.7
    }
}

/// Runs an end-to-end mRNA-therapeutic design (feature 30, mRNA half).
///
/// Chains codon optimisation, optional uridine depletion, UTR
/// selection, structure minimisation and end-modification selection,
/// then assembles a validated [`MrnaConstruct`] and returns a full
/// [`MrnaDesignReport`].
///
/// # Errors
/// Propagates [`crate::error::GeneditingError`] from any pass — most
/// often [`crate::error::GeneditingError::InvalidTarget`] for a CDS
/// that is not a valid coding sequence.
pub fn design_mrna(req: &MrnaDesignRequest) -> Result<MrnaDesignReport> {
    let mut notes: Vec<String> = Vec::new();

    // 1) Codon-optimise.
    let codon = optimize_cds(&req.cds, req.host)?;
    notes.push(format!(
        "Codon-optimised for {}: CAI {:.3} -> {:.3}.",
        req.host.name(),
        codon.cai_before,
        codon.cai_after
    ));
    let mut cds = codon.optimized_cds.clone();
    let cai_before = codon.cai_before;
    let mut cai_after = codon.cai_after;

    // 2) Optional uridine depletion.
    if req.deplete_uridine {
        let u = minimize_uridine(&cds, Some(req.host))?;
        notes.push(format!(
            "Uridine depletion: U fraction {:.3} -> {:.3} (for an {} construct).",
            u.uridine_before,
            u.uridine_after,
            req.nucleoside.name()
        ));
        cds = u.optimized_cds;
        // Recompute CAI after the uridine pass (it may shift slightly).
        cai_after = crate::mrna::codon::cds_cai(&cds, req.host)?;
    }

    // 3) Reference UTRs.
    let utr5 = reference_utr5();
    let utr3 = reference_utr3();
    let u5 = analyze_utr5(&utr5, &cds)?;
    let u3 = analyze_utr3(&utr3)?;
    notes.push(format!(
        "Reference UTRs: Kozak score {:.2} ({}), 3'UTR stability {:.2} ({} ARE/destabilisers).",
        u5.kozak_score,
        if u5.strong_kozak { "strong" } else { "weak" },
        u3.stability_score,
        u3.destabilizer_count(),
    ));

    // 4) Structure minimisation near the start codon.
    let smin = minimize_start_structure(&utr5, &cds, req.structure_scan_codons)?;
    notes.push(format!(
        "Start-region structure minimisation: openness {:.2} -> {:.2} ({} variants tried).",
        smin.openness_before, smin.openness_after, smin.variants_tried,
    ));
    cds = smin.optimized_cds;
    let start_openness = smin.openness_after;

    // 5) End-modification recommendations.
    let poly_a = recommend_poly_a(req.use_case);
    let cap = recommend_cap(req.use_case);
    notes.push(format!(
        "Poly-A {} nt{}; cap {}.",
        poly_a.length,
        if poly_a.segmented { " (segmented)" } else { "" },
        cap.cap.name(),
    ));

    // 6) Assemble + validate the construct.
    let construct = MrnaConstructBuilder::new()
        .cap(cap.cap)
        .utr5(&utr5)
        .cds(&cds)
        .utr3(&utr3)
        .poly_a(poly_a.length)
        .build()?;

    // 7) Modified-nucleoside plan over the whole transcript.
    let modification = plan_modification(&construct.transcript(), req.nucleoside)?;
    notes.push(modification.rationale.clone());

    // Overall score: a transparent blend of the per-feature scores.
    let overall = (0.25 * cai_after
        + 0.25 * start_openness
        + 0.20 * u5.kozak_score
        + 0.20 * u3.stability_score
        + 0.10 * if modification.reduces_immunogenicity { 1.0 } else { 0.4 })
    .clamp(0.0, 1.0);

    Ok(MrnaDesignReport {
        construct,
        cai_before,
        cai_after,
        uridine_fraction: modification.uridine_fraction,
        start_openness,
        kozak_score: u5.kozak_score,
        utr3_stability: u3.stability_score,
        poly_a_len: poly_a.length,
        cap: cap.cap,
        modification,
        overall_score: overall,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // A modest CDS: ATG + a few codons + stop.
    fn cds() -> Vec<u8> {
        b"ATGGCCCTGCTGGAAGAATAA".to_vec()
    }

    #[test]
    fn designs_an_mrna_end_to_end() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::Vaccine);
        let report = design_mrna(&req).unwrap();
        // The construct round-trips as a valid mRNA.
        assert!(report.construct.codon_count() >= 3);
        assert!(!report.construct.transcript().is_empty());
        // Every score in range.
        assert!((0.0..=1.0).contains(&report.overall_score));
        assert!((0.0..=1.0).contains(&report.start_openness));
        assert!(report.cai_after > 0.0);
    }

    #[test]
    fn notes_describe_each_pass() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::ProteinReplacement);
        let report = design_mrna(&req).unwrap();
        assert!(report.notes.iter().any(|n| n.contains("Codon-optimised")));
        assert!(report.notes.iter().any(|n| n.contains("Poly-A")));
        assert!(report.notes.iter().any(|n| n.contains("structure")));
    }

    #[test]
    fn vaccine_uses_cleancap() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::Vaccine);
        let report = design_mrna(&req).unwrap();
        assert_eq!(report.cap, CapType::CleanCap);
        assert_eq!(report.poly_a_len, 100);
    }

    #[test]
    fn protein_replacement_gets_a_longer_tail() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::ProteinReplacement);
        let report = design_mrna(&req).unwrap();
        assert!(report.poly_a_len >= 140);
    }

    #[test]
    fn uridine_depletion_can_be_disabled() {
        let mut req = MrnaDesignRequest::new(cds(), MrnaUseCase::Vaccine);
        req.deplete_uridine = false;
        let report = design_mrna(&req).unwrap();
        assert!(!report.notes.iter().any(|n| n.contains("Uridine depletion")));
    }

    #[test]
    fn uridine_depletion_is_noted_when_enabled() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::Vaccine);
        let report = design_mrna(&req).unwrap();
        assert!(report.notes.iter().any(|n| n.contains("Uridine depletion")));
    }

    #[test]
    fn rejects_an_invalid_cds() {
        // Not a multiple of 3.
        let req = MrnaDesignRequest::new(b"ATGCT".to_vec(), MrnaUseCase::Vaccine);
        assert!(design_mrna(&req).is_err());
    }

    #[test]
    fn m1psi_construct_reduces_immunogenicity() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::Vaccine);
        let report = design_mrna(&req).unwrap();
        assert!(report.modification.reduces_immunogenicity);
    }

    #[test]
    fn high_quality_flag_tracks_overall_score() {
        let req = MrnaDesignRequest::new(cds(), MrnaUseCase::Vaccine);
        let report = design_mrna(&req).unwrap();
        assert_eq!(report.is_high_quality(), report.overall_score >= 0.7);
    }
}
