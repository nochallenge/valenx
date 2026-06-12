//! Feature 6 — coding-mRNA sequence design.
//!
//! Given a protein to express, build a complete five-part therapeutic
//! mRNA: a 5′ cap, a 5′UTR, the codon-optimised CDS, a 3′UTR and a
//! poly-A tail. This module is an orchestration layer over
//! [`valenx_genediting::design_mrna`] — the construct model, the codon
//! optimiser, the structure-minimisation pass and the uridine planner
//! are *that* crate's code, never re-implemented here.
//!
//! What this module adds is the **protein → CDS** front end the
//! `valenx-genediting` driver does not provide.
//!
//! ## The CDS optimiser — joint by default
//!
//! As of the RNA-design depth pass the CDS is produced by the
//! **LinearDesign joint optimiser** ([`crate::lineardesign`]) — it
//! optimises codon usage *and* secondary-structure stability together,
//! the way commercial mRNA design (LinearDesign, Zhang *et al.*,
//! *Nature* 2023) does, rather than the legacy "pick the host-optimal
//! codon, then fold" two-step. The caller picks the
//! stability/efficiency trade-off with
//! [`CodingDesignParams::joint_lambda`]. Setting
//! [`CodingDesignParams::joint_optimize`] to `false` falls back to the
//! legacy `valenx-bioseq` `codon_optimize` starting CDS.
//!
//! Either way a start codon is prepended if the protein lacks a leading
//! `M` and a stop codon appended if there is no trailing stop, so the
//! result is the valid `AUG…stop` CDS the construct builder requires.
//!
//! ## v1 scope
//!
//! With joint optimisation on, the CDS is the LinearDesign lattice-DP
//! optimum and the construct is assembled directly from the
//! `valenx-genediting` building blocks (reference UTRs, cap / poly-A
//! recommendation, the validating [`valenx_genediting::MrnaConstructBuilder`])
//! — the genediting driver's own max-CAI codon pass is bypassed so it
//! cannot overwrite the joint optimum. With joint optimisation off the
//! full [`valenx_genediting::design_mrna`] workflow runs unchanged. The
//! encoded protein is invariant under either path.

use crate::design::{DesignKind, RnaDesign};
use crate::error::{Result, RnaDesignError};
use crate::lineardesign::{linear_design, LinearDesignRequest};
use valenx_bioseq::cloning::codon_opt::{codon_optimize, Host};
use valenx_bioseq::ops::translate::GeneticCode;
use valenx_bioseq::{Seq, SeqKind};
use valenx_genediting::mrna::codon::ExpressionHost;
use valenx_genediting::mrna::construct::MrnaConstructBuilder;
use valenx_genediting::mrna::tailcap::{recommend_cap, recommend_poly_a, MrnaUseCase};
use valenx_genediting::mrna::uridine::ModifiedNucleoside;
use valenx_genediting::mrna::utr::{reference_utr3, reference_utr5};
use valenx_genediting::{design_mrna, MrnaDesignRequest};

/// Parameters for [`design_coding`].
#[derive(Copy, Clone, Debug)]
pub struct CodingDesignParams {
    /// The expression host (drives codon usage).
    pub host: ExpressionHost,
    /// The intended therapeutic use (drives poly-A / cap selection).
    pub use_case: MrnaUseCase,
    /// `true` to run the uridine-depletion pass (recommended for an
    /// `m1Ψ` therapeutic construct). Only consulted on the legacy
    /// (`joint_optimize == false`) path.
    pub deplete_uridine: bool,
    /// The modified nucleoside the physical RNA will be made with.
    pub nucleoside: ModifiedNucleoside,
    /// `true` (the default) to produce the CDS with the **LinearDesign
    /// joint optimiser** — codon usage *and* mRNA structural stability
    /// optimised together. `false` falls back to the legacy
    /// codon-optimise-then-fold pipeline via
    /// [`valenx_genediting::design_mrna`].
    pub joint_optimize: bool,
    /// The LinearDesign trade-off weight `λ ≥ 0` (consulted only when
    /// `joint_optimize` is `true`). `0` → most stable mRNA; large →
    /// most CAI-optimal; intermediate → the Pareto trade-off. The
    /// default `0.5` favours a stable, well-translated balance.
    pub joint_lambda: f64,
}

impl CodingDesignParams {
    /// Therapeutic defaults for a use case: human host, `m1Ψ`,
    /// LinearDesign joint optimisation on at `λ = 0.5`.
    pub fn new(use_case: MrnaUseCase) -> Self {
        CodingDesignParams {
            host: ExpressionHost::Human,
            use_case,
            deplete_uridine: true,
            nucleoside: ModifiedNucleoside::N1MethylPseudouridine,
            joint_optimize: true,
            joint_lambda: 0.5,
        }
    }
}

/// Maps the `valenx-genediting` host enum to the `valenx-bioseq` one.
fn bioseq_host(host: ExpressionHost) -> Host {
    match host {
        ExpressionHost::Human => Host::Human,
        ExpressionHost::EColi => Host::EColi,
    }
}

/// Designs a coding mRNA that expresses `protein` (feature 6).
///
/// With [`CodingDesignParams::joint_optimize`] (the default) the CDS is
/// produced by the LinearDesign joint optimiser, then the five-part
/// construct is assembled from the `valenx-genediting` building blocks;
/// with it `false` the full [`valenx_genediting::design_mrna`] workflow
/// runs. Either way the result is an [`RnaDesign`] whose `construct` is
/// the assembled five-part mRNA and whose `cds_span` locates the CDS.
///
/// # Errors
/// - [`RnaDesignError::Goal`] for an empty protein.
/// - [`RnaDesignError::Upstream`] if the codon optimiser, the
///   LinearDesign DP or the construct builder rejects the input.
pub fn design_coding(protein: &[u8], params: CodingDesignParams) -> Result<RnaDesign> {
    if protein.is_empty() {
        return Err(RnaDesignError::goal("protein", "protein sequence is empty"));
    }
    if params.joint_optimize {
        design_coding_joint(protein, params)
    } else {
        design_coding_legacy(protein, params)
    }
}

/// The LinearDesign joint-optimiser path: optimise the CDS for codon
/// usage *and* structural stability together, then assemble the
/// construct directly from the `valenx-genediting` building blocks.
fn design_coding_joint(protein: &[u8], params: CodingDesignParams) -> Result<RnaDesign> {
    // The LinearDesign optimiser adds the leading M / trailing stop
    // itself; pass the bare protein.
    let prepended_met = protein.first().map(|b| b.to_ascii_uppercase()) != Some(b'M');
    let appended_stop = protein.last().map(|b| b.to_ascii_uppercase()) != Some(b'*');

    let req = LinearDesignRequest::new(protein.to_vec(), params.host, params.joint_lambda);
    let ld = linear_design(&req)?;
    let cds = ld.cds.clone();

    // Assemble the five-part construct with the genediting building
    // blocks: reference UTRs, a recommended cap and poly-A tail.
    let utr5 = reference_utr5();
    let utr3 = reference_utr3();
    let poly_a = recommend_poly_a(params.use_case);
    let cap = recommend_cap(params.use_case);
    let construct = MrnaConstructBuilder::new()
        .cap(cap.cap)
        .utr5(&utr5)
        .cds(&cds)
        .utr3(&utr3)
        .poly_a(poly_a.length)
        .build()?;

    let transcript = construct.transcript();
    let cds_start = construct.cds_start();
    let cds_end = cds_start + construct.cds.len();

    let mut notes = Vec::new();
    notes.push(format!(
        "Designed a coding mRNA expressing a {}-residue protein for {} ({} host).",
        protein.len(),
        params.use_case.name(),
        params.host.name(),
    ));
    notes.push(format!(
        "CDS optimised by the LinearDesign joint optimiser (codon usage + structural \
         stability together) at λ={:.2}: predicted MFE {:.1} kcal/mol, CAI {:.3}{}.",
        params.joint_lambda,
        ld.mfe,
        ld.cai,
        if ld.exact {
            ", exact lattice optimum"
        } else {
            ""
        },
    ));
    if prepended_met {
        notes.push(
            "The protein did not begin with methionine — a start codon (AUG/Met) was \
             prepended to make a translatable CDS."
                .to_string(),
        );
    }
    if appended_stop {
        notes.push("A stop codon was appended to terminate translation.".to_string());
    }
    notes.push(format!(
        "Construct: cap {}, reference 5'/3' UTRs, poly-A {} nt; transcript {} nt \
         (CDS at {}..{}).",
        cap.cap.name(),
        poly_a.length,
        transcript.len(),
        cds_start,
        cds_end,
    ));
    notes.push(
        "This is an in-silico mRNA design — the MFE and CAI are predictions from an energy \
         model and codon-usage tables; the construct must be validated in the wet lab."
            .to_string(),
    );

    Ok(RnaDesign {
        sequence: transcript,
        kind: DesignKind::Coding,
        cds_span: Some((cds_start, cds_end)),
        construct: Some(construct),
        notes,
    })
}

/// The legacy path: codon-optimise the protein into a starting CDS,
/// then run the full [`valenx_genediting::design_mrna`] workflow.
fn design_coding_legacy(protein: &[u8], params: CodingDesignParams) -> Result<RnaDesign> {
    // 1) Build the CDS protein: prepend M if absent, append a stop if
    //    absent. `codon_optimize` reads `*` as the stop residue.
    let mut cds_protein: Vec<u8> = Vec::with_capacity(protein.len() + 2);
    let prepended_met = protein.first().map(|b| b.to_ascii_uppercase()) != Some(b'M');
    if prepended_met {
        cds_protein.push(b'M');
    }
    cds_protein.extend(protein.iter().map(|b| b.to_ascii_uppercase()));
    let appended_stop = protein.last().map(|b| b.to_ascii_uppercase()) != Some(b'*');
    if appended_stop {
        cds_protein.push(b'*');
    }

    let protein_seq = Seq::new(SeqKind::Protein, &cds_protein)?;

    // 2) Codon-optimise the protein into a starting DNA CDS.
    let cds_dna = codon_optimize(&protein_seq, bioseq_host(params.host))?;

    // 3) Run the full mRNA-design workflow.
    let mut req = MrnaDesignRequest::new(cds_dna.as_bytes().to_vec(), params.use_case);
    req.host = params.host;
    req.deplete_uridine = params.deplete_uridine;
    req.nucleoside = params.nucleoside;
    let report = design_mrna(&req)?;

    let construct = report.construct.clone();
    let transcript = construct.transcript();
    let cds_start = construct.cds_start();
    let cds_end = cds_start + construct.cds.len();

    let mut notes = Vec::new();
    notes.push(format!(
        "Designed a coding mRNA expressing a {}-residue protein for {} ({} host).",
        protein.len(),
        params.use_case.name(),
        params.host.name(),
    ));
    if prepended_met {
        notes.push(
            "The protein did not begin with methionine — a start codon (AUG/Met) was \
             prepended to make a translatable CDS."
                .to_string(),
        );
    }
    if appended_stop {
        notes.push("A stop codon was appended to terminate translation.".to_string());
    }
    notes.push(format!(
        "CDS codon-adaptation index {:.3} -> {:.3}; start-region openness {:.2}; \
         transcript {} nt (CDS at {}..{}).",
        report.cai_before,
        report.cai_after,
        report.start_openness,
        transcript.len(),
        cds_start,
        cds_end,
    ));
    // Carry the genediting driver's own per-pass notes through.
    notes.extend(report.notes.iter().cloned());
    notes.push(
        "This is an in-silico mRNA design — a strong predicted construct from heuristic \
         codon / structure models; it must be validated in the wet lab."
            .to_string(),
    );

    Ok(RnaDesign {
        sequence: transcript,
        kind: DesignKind::Coding,
        cds_span: Some((cds_start, cds_end)),
        construct: Some(construct),
        notes,
    })
}

/// Translates a coding design's CDS back to protein, for a round-trip
/// check. Returns the one-letter amino-acid sequence.
///
/// # Errors
/// [`RnaDesignError::Invalid`] if the design is not a coding design.
pub fn translate_cds(design: &RnaDesign) -> Result<Vec<u8>> {
    let cds = design
        .cds()
        .ok_or_else(|| RnaDesignError::invalid("design", "not a coding design"))?;
    let code = GeneticCode::standard();
    // The CDS is RNA; translate codon by codon (reverse-transcribe each).
    let mut protein = Vec::with_capacity(cds.len() / 3);
    for codon in cds.chunks(3) {
        if codon.len() < 3 {
            break;
        }
        let dna: [u8; 3] = [
            rna_to_dna(codon[0]),
            rna_to_dna(codon[1]),
            rna_to_dna(codon[2]),
        ];
        protein.push(code.translate_codon(&dna));
    }
    Ok(protein)
}

/// Reverse-transcribes a single RNA base to DNA (`U → T`).
fn rna_to_dna(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_a_coding_mrna() {
        let d = design_coding(b"MKVLA", CodingDesignParams::new(MrnaUseCase::Vaccine)).unwrap();
        assert!(d.is_coding());
        assert!(d.construct.is_some());
        assert!(d.cds_span.is_some());
        assert!(!d.sequence.is_empty());
        assert!(!d.notes.is_empty());
    }

    #[test]
    fn prepends_start_codon_when_missing() {
        // A protein not starting with M gets an M prepended.
        let d = design_coding(b"KVLA", CodingDesignParams::new(MrnaUseCase::Vaccine)).unwrap();
        assert!(d.notes.iter().any(|n| n.contains("start codon")));
        let protein = translate_cds(&d).unwrap();
        // Translated CDS: M K V L A * — starts with M.
        assert_eq!(protein.first(), Some(&b'M'));
    }

    #[test]
    fn keeps_existing_methionine() {
        let d = design_coding(b"MKVLA", CodingDesignParams::new(MrnaUseCase::Vaccine)).unwrap();
        // Protein already starts with M — no extra M note.
        assert!(!d.notes.iter().any(|n| n.contains("did not begin")));
        let protein = translate_cds(&d).unwrap();
        // M K V L A * — one M, not two.
        assert_eq!(&protein[..5], b"MKVLA");
    }

    #[test]
    fn cds_translates_back_to_the_protein() {
        let d = design_coding(b"MKVLAGD", CodingDesignParams::new(MrnaUseCase::Vaccine)).unwrap();
        let protein = translate_cds(&d).unwrap();
        // CDS protein = M K V L A G D * — strip the trailing stop.
        assert_eq!(&protein[..protein.len() - 1], b"MKVLAGD");
        assert_eq!(protein.last(), Some(&b'*'));
    }

    #[test]
    fn cds_span_locates_the_cds() {
        let d = design_coding(b"MKVL", CodingDesignParams::new(MrnaUseCase::Vaccine)).unwrap();
        let (s, e) = d.cds_span.unwrap();
        // CDS is M K V L * = 5 codons = 15 nt.
        assert_eq!(e - s, 15);
        // The CDS subsequence starts with the start codon.
        assert_eq!(&d.cds().unwrap()[..3], b"AUG");
    }

    #[test]
    fn rejects_empty_protein() {
        assert!(design_coding(b"", CodingDesignParams::new(MrnaUseCase::Vaccine)).is_err());
    }

    #[test]
    fn uridine_depletion_can_be_disabled() {
        let mut params = CodingDesignParams::new(MrnaUseCase::Vaccine);
        params.deplete_uridine = false;
        let d = design_coding(b"MKVLA", params).unwrap();
        assert!(d.is_coding());
    }

    #[test]
    fn translate_cds_rejects_non_coding() {
        let non_coding = RnaDesign {
            sequence: b"GGGGCCCC".to_vec(),
            kind: DesignKind::Coding,
            cds_span: None,
            construct: None,
            notes: Vec::new(),
        };
        assert!(translate_cds(&non_coding).is_err());
    }
}
