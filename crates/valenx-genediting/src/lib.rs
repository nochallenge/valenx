//! # valenx-genediting — gene-editing & mRNA-therapeutic design
//!
//! An *extension* that
//! unifies the scattered editing and mRNA pieces of the
//! computational-biology crates into proper design workflows.
//!
//! This is **research and therapeutic gene-editing + mRNA design
//! tooling** — the same software category as CHOPCHOP, CRISPOR,
//! PrimeDesign, BE-Designer and Benchling's design tools. It designs
//! guide RNAs, donor templates, pegRNAs and mRNA constructs *in
//! silico* and ranks them by transparent rule-based scores. It does
//! not perform edits, does not call out to wet-lab hardware, and makes
//! no claims about clinical outcomes.
//!
//! It builds on three sibling crates and never duplicates them:
//!
//! - [`valenx_genomics`] — its CRISPR PAM scanner
//!   ([`valenx_genomics::crispr::guide`]), on-target score and
//!   off-target enumerator ([`valenx_genomics::crispr::offtarget`]) are
//!   reused directly as the substrate for guide design.
//! - [`valenx_bioseq`] — the [`Seq`](valenx_bioseq::Seq)
//!   type, the genetic-code translation tables and the
//!   codon-optimisation / CAI machinery behind the mRNA module.
//! - [`valenx_rnastruct`] — the Zuker minimum-free-energy
//!   folder behind the mRNA secondary-structure check.
//!
//! ## What it does
//!
//! - **CRISPR nuclease editing** ([`crispr`]) — a nuclease database
//!   (SpCas9, SpCas9-NG, SaCas9, Cas12a / Cpf1, Cas12f, Cas13, xCas9),
//!   a guide-RNA design workflow, NHEJ knockout strategy design, HDR
//!   knock-in donor-template design and multiplex / gRNA-array design.
//! - **Base editing** ([`base_edit`]) — a CBE / ABE base-editor
//!   database, a reachable-edit finder, SNV-correcting guide design,
//!   editing-window + bystander analysis and a product-purity
//!   heuristic.
//! - **Prime editing** ([`prime_edit`]) — a PE2 / PE3 / PE3b / PEmax
//!   database, pegRNA design (spacer + PBS + RT template), a
//!   PBS-/RT-length optimisation scan, PE3 / PE3b nicking-guide design
//!   and a prime-editing efficiency heuristic.
//! - **mRNA therapeutic design** ([`mrna`]) — an mRNA construct model
//!   (cap, 5'UTR, CDS, 3'UTR, poly-A), CDS codon optimisation, Kozak-
//!   aware 5'UTR design, AU-rich-element-aware 3'UTR design,
//!   secondary-structure minimisation, uridine / pseudouridine
//!   planning, poly-A + cap-analog selection and self-amplifying /
//!   circular-RNA construct design.
//! - **Gene therapy & safety** ([`therapy`]) — gene-therapy
//!   expression-cassette design with AAV / lentivirus payload checks,
//!   informational delivery-vector planning and a safety-screen
//!   aggregator.
//! - **Workflow** ([`workflow`]) — an edit-strategy advisor, a
//!   variant-correction planner, a batch design driver and a typed
//!   LLM / MCP-controllable request / response surface.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, GeneditingError>`](error::GeneditingError). The error
//! type carries stable [`code`](error::GeneditingError::code) and
//! [`category`](error::GeneditingError::category) accessors for
//! telemetry.
//!
//! ## v1 scope — honest caveats
//!
//! This is a real working v1, not production parity with the reference
//! design servers. Each module documents its own simplifications
//! inline; the load-bearing ones are:
//!
//! - **Every efficiency / outcome predictor is a transparent
//!   rule-based heuristic, not a trained model.** Guide on-target /
//!   off-target scoring comes from `valenx-genomics`' Doench- and
//!   CFD-*style* feature-weighted heuristics; base-editing product
//!   purity, prime-editing efficiency and mRNA structure penalties are
//!   all hand-written feature scores with the right qualitative
//!   ranking. The real-world tools (DeepCBE, BE-Hive, the PrimeDesign
//!   efficiency models, the CodonBERT-class mRNA models) use trained
//!   neural networks; reimplementing trained weights is excluded by
//!   the project's standing "no trained-weights" rule, so every score
//!   here is documented as a heuristic.
//! - The nuclease / base-editor / prime-editor databases carry
//!   *representative published* PAM, window and chemistry parameters,
//!   not an exhaustive catalogue of every engineered variant.
//! - Donor-template and pegRNA designers operate on a linear
//!   reference window the caller supplies; they do not fetch a genome,
//!   resolve isoforms or model chromatin accessibility.
//! - mRNA secondary-structure minimisation reuses the `valenx-rnastruct`
//!   Zuker folder, whose Turner-2004 parameters are a faithful subset
//!   (see that crate's note), and scans a bounded set of synonymous
//!   variants — it is not an exhaustive sequence optimiser.

#![forbid(unsafe_code)]

pub mod base_edit;
pub mod crispr;
pub mod error;
pub mod mrna;
pub mod prime_edit;
pub mod therapy;
pub mod workflow;

mod sequtil;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, GeneditingError, Result};

pub use crispr::donor::{design_hdr_donor, DonorTemplate, HdrDonorRequest};
pub use crispr::donor_opt::{
    design_optimized_hdr_donor, recommend_arm_length, scan_splice_sites, OptimizedDonor,
    OptimizedDonorRequest, SpliceWarning,
};
pub use crispr::guide_design::{design_guides, GuideCandidate, GuideDesignRequest};
pub use crispr::knockout::{design_knockout, KnockoutRequest, KnockoutStrategy};
pub use crispr::multiplex::{design_multiplex, MultiplexDesign, MultiplexRequest};
pub use crispr::nuclease::{Nuclease, NucleaseClass, NucleaseId};
pub use crispr::offtarget_fm::{
    find_off_targets_genome, off_target_report, GenomeIndex, IndexedContig, OffTargetReport,
};

pub use base_edit::design::{design_base_edit, BaseEditDesign, BaseEditRequest};
pub use base_edit::editor::{BaseEditor, BaseEditorClass, BaseEditorId};

pub use prime_edit::editor::{PrimeEditor, PrimeEditorId};
pub use prime_edit::pegrna::{design_pegrna, PegRna, PegRnaRequest};

pub use mrna::construct::{MrnaConstruct, MrnaConstructBuilder};
pub use mrna::design::{design_mrna, MrnaDesignReport, MrnaDesignRequest};

pub use therapy::cassette::{design_cassette, CassetteRequest, ExpressionCassette};
pub use therapy::safety::{
    aggregate_safety, safety_screen, EditSafetyReport, EditScreenRequest, SafetyReport,
};
pub use therapy::safety_db::{
    cancer_driver_genes, essential_genes, safe_harbor_loci, ReferenceGeneDatabase,
    ReferenceGeneList, ReferenceListKind,
};

pub use workflow::advisor::{advise_strategy, EditApproach, StrategyAdvice};
pub use workflow::driver::{
    handle_request, run_editing_design, EditingReport, EditingRequest, GeneditingRequest,
    GeneditingResponse,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_error_reexported() {
        let e = GeneditingError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }
}
