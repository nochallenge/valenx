//! Group B — sequence design (features 5–9).
//!
//! The design step turns a [`crate::goal::DesignGoal`] into a concrete
//! RNA *candidate*: an actual sequence over `A C G U`. Five designers,
//! one per goal shape:
//!
//! - [`structural`] (feature 5) — inverse-fold to a target dot-bracket
//!   via `valenx-rnastruct`, multi-start, keeping the best candidate(s).
//! - [`coding`] (feature 6) — assemble a five-part mRNA construct and
//!   codon-optimise the CDS via `valenx-genediting`.
//! - [`riboswitch`] (feature 7) — design a sequence that can adopt two
//!   target structures (ligand-bound vs free).
//! - [`motif`] (feature 8) — start from a known aptamer / ribozyme /
//!   hairpin scaffold template.
//! - [`regulatory`] (feature 9) — design 5′ / 3′ UTRs and the Kozak
//!   context, delegating to `valenx-genediting` where it fits.
//!
//! Every designer produces an [`RnaDesign`] — the shared candidate
//! type this module defines. An [`RnaDesign`] carries the designed RNA
//! sequence, the design provenance ([`DesignKind`]), the coding-region
//! span when there is one, and the assembled [`MrnaConstruct`] for a
//! coding design.
//!
//! ## v1 scope — honest framing
//!
//! An [`RnaDesign`] is a *predicted* candidate. Structural designs are
//! produced by an inverse-folding adaptive walk against a
//! nearest-neighbor energy model; the achieved structure is that
//! model's MFE, not a measured structure. Coding designs use the
//! `valenx-genediting` heuristic codon / structure passes. The
//! candidate is strong and self-consistent in silico — it is never a
//! guarantee of in-vivo behaviour. See the crate root.

pub mod coding;
pub mod motif;
pub mod regulatory;
pub mod riboswitch;
pub mod structural;

use serde::{Deserialize, Serialize};
use valenx_genediting::mrna::construct::MrnaConstruct;
use valenx_rnastruct::{RnaSeq, Structure};

pub use coding::{design_coding, CodingDesignParams};
pub use motif::{motif_scaffold, MotifScaffold, ScaffoldKind};
pub use regulatory::{design_regulatory, RegulatoryDesign, RegulatoryParams};
pub use riboswitch::{design_riboswitch, RiboswitchDesign, TwoStateParams};
pub use structural::{design_structural, StructuralDesignParams};

/// How an [`RnaDesign`] was produced — its provenance.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum DesignKind {
    /// An inverse-folded structural RNA. Carries the target structure
    /// the design was inverse-folded against.
    Structural {
        /// The target secondary structure.
        target: Structure,
    },
    /// A codon-optimised, assembled coding mRNA.
    Coding,
    /// A two-state riboswitch design. Carries both target structures.
    Riboswitch {
        /// The structure adopted in the ligand-free state.
        free_target: Structure,
        /// The structure adopted in the ligand-bound state.
        bound_target: Structure,
    },
    /// A design seeded from a known functional-motif scaffold.
    MotifScaffold {
        /// The scaffold the design started from.
        scaffold: String,
    },
}

impl DesignKind {
    /// A short human-readable name for the design provenance.
    pub fn name(&self) -> &'static str {
        match self {
            DesignKind::Structural { .. } => "structural (inverse-folded)",
            DesignKind::Coding => "coding mRNA",
            DesignKind::Riboswitch { .. } => "riboswitch (two-state)",
            DesignKind::MotifScaffold { .. } => "functional-motif scaffold",
        }
    }
}

/// A concrete RNA design candidate — the shared output of every
/// designer in this module.
///
/// The sequence is stored as RNA (`A C G U`). For a coding design the
/// `construct` field carries the assembled five-part mRNA and
/// `cds_span` locates the CDS within the transcript; for a structural
/// design those are `None`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RnaDesign {
    /// The designed RNA sequence (`A C G U`).
    pub sequence: Vec<u8>,
    /// How the design was produced.
    pub kind: DesignKind,
    /// For a coding design, the `[start, end)` span of the CDS within
    /// [`sequence`](Self::sequence); `None` for a non-coding design.
    pub cds_span: Option<(usize, usize)>,
    /// For a coding design, the assembled five-part mRNA construct.
    pub construct: Option<MrnaConstruct>,
    /// Human-readable notes describing the design (one line per major
    /// decision the designer made).
    pub notes: Vec<String>,
}

impl RnaDesign {
    /// The design length in nucleotides.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// `true` when the design carries no nucleotides (never produced by
    /// a successful designer).
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// The designed sequence as a `&str` — valid UTF-8 because the
    /// designers only ever emit `A C G U`.
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }

    /// `true` when this is a coding design (it has a CDS and a
    /// construct).
    pub fn is_coding(&self) -> bool {
        self.construct.is_some()
    }

    /// The CDS subsequence of a coding design, or `None` for a
    /// non-coding design (or a coding design with an out-of-range span).
    pub fn cds(&self) -> Option<&[u8]> {
        let (s, e) = self.cds_span?;
        self.sequence.get(s..e)
    }

    /// Builds an [`RnaSeq`] (the `valenx-rnastruct` sequence wrapper)
    /// from this design's sequence, for folding.
    ///
    /// # Errors
    /// [`crate::error::RnaDesignError::Upstream`] if the sequence is not
    /// foldable (it always should be — designers emit only `A C G U`).
    pub fn to_rna_seq(&self) -> crate::error::Result<RnaSeq> {
        Ok(RnaSeq::parse(&self.sequence)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn demo_design() -> RnaDesign {
        RnaDesign {
            sequence: b"GGGGAAAACCCC".to_vec(),
            kind: DesignKind::Structural {
                target: Structure::from_dot_bracket("((((....))))").unwrap(),
            },
            cds_span: None,
            construct: None,
            notes: vec!["inverse-folded".to_string()],
        }
    }

    #[test]
    fn design_basics() {
        let d = demo_design();
        assert_eq!(d.len(), 12);
        assert!(!d.is_empty());
        assert_eq!(d.sequence_str(), "GGGGAAAACCCC");
        assert!(!d.is_coding());
        assert!(d.cds().is_none());
    }

    #[test]
    fn design_folds() {
        let d = demo_design();
        let seq = d.to_rna_seq().unwrap();
        assert_eq!(seq.len(), 12);
    }

    #[test]
    fn design_kind_names() {
        assert_eq!(DesignKind::Coding.name(), "coding mRNA");
        let s = DesignKind::Structural {
            target: Structure::empty(4),
        };
        assert!(s.name().contains("structural"));
    }

    #[test]
    fn cds_span_extracts_subsequence() {
        let mut d = demo_design();
        d.cds_span = Some((4, 8));
        assert_eq!(d.cds(), Some(&b"AAAA"[..]));
        // Out-of-range span yields None.
        d.cds_span = Some((4, 99));
        assert_eq!(d.cds(), None);
    }
}
