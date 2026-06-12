//! Feature 3 — NHEJ knockout strategy design.
//!
//! A gene knockout by non-homologous end joining (NHEJ) works because
//! the indels that NHEJ leaves at a Cas9 cut are usually *not* a
//! multiple of three, so they shift the reading frame and abolish the
//! protein downstream. A good knockout guide therefore cuts:
//!
//! - **early in the coding sequence** — an early frameshift truncates
//!   the whole protein; a frameshift in the last exon may leave a
//!   functional product;
//! - **inside (or just before) a functional domain** — disrupting a
//!   catalytic / binding domain is a fallback when the very 5′ end is
//!   untargetable;
//! - at a site where a frameshift cannot be "rescued" by a nearby
//!   in-frame ATG or by skipping the targeted exon.
//!
//! This module takes a [`GeneModel`] (a coding sequence laid out as
//! exons, optionally annotated with functional domains), designs guides
//! against each exon with [`crate::crispr::guide_design`], and ranks
//! them by a transparent **knockout score** that rewards an early,
//! domain-disrupting, high-efficiency cut.
//!
//! ## v1 scope
//!
//! The knockout score is a transparent feature-weighted heuristic
//! (position-in-CDS, domain overlap, on-target efficiency). It does
//! not predict the *indel-length distribution* of a particular cut
//! site — that is a trained-model problem (inDelphi / FORECasT) the
//! "no trained-weights" rule excludes; this module reasons about
//! *where* to cut, not the exact repair spectrum. For an
//! indel-spectrum classifier over actual edited reads, see
//! `valenx_genomics::crispr::editing` (the CRISPResso-class analyser),
//! which the [`crate::workflow`] driver wraps.

use crate::crispr::guide_design::{design_guides, GuideCandidate, GuideDesignRequest};
use crate::crispr::nuclease::NucleaseId;
use crate::error::{GeneditingError, Result};
use crate::sequtil::is_acgt;
use serde::{Deserialize, Serialize};

/// One exon of a gene model — a contiguous coding segment.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Exon {
    /// 0-based start of the exon within the gene-model sequence.
    pub start: usize,
    /// 0-based end (exclusive) of the exon within the sequence.
    pub end: usize,
    /// 1-based exon ordinal (exon 1 is the 5′-most coding exon).
    pub ordinal: usize,
}

impl Exon {
    /// Length of the exon in base pairs.
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// `true` when the exon spans no bases.
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// A protein functional domain projected onto the coding sequence.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FunctionalDomain {
    /// Domain name (e.g. `"kinase"`, `"DNA-binding"`).
    pub name: String,
    /// 0-based start of the domain within the gene-model sequence.
    pub start: usize,
    /// 0-based end (exclusive) of the domain.
    pub end: usize,
}

/// A gene model: a coding sequence laid out as exons, with optional
/// functional-domain annotations.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GeneModel {
    /// The full gene-model DNA sequence (forward strand, 5′→3′). For a
    /// knockout workflow this is the spliced CDS; coordinates of
    /// [`Exon`]s and [`FunctionalDomain`]s index into it.
    pub sequence: Vec<u8>,
    /// The coding exons, 5′→3′.
    pub exons: Vec<Exon>,
    /// Functional domains projected onto `sequence` (may be empty).
    pub domains: Vec<FunctionalDomain>,
}

impl GeneModel {
    /// Builds a gene model from a CDS plus exon boundaries.
    ///
    /// `boundaries` are the 0-based exon start offsets (the first must
    /// be `0`); each exon runs to the next boundary, the last to the
    /// end of the CDS.
    ///
    /// # Errors
    /// [`GeneditingError::InvalidTarget`] for a non-ACGT CDS, an empty
    /// boundary list, a first boundary other than `0`, or boundaries
    /// that are not strictly increasing / in range.
    pub fn from_cds(cds: impl Into<Vec<u8>>, boundaries: &[usize]) -> Result<Self> {
        let sequence: Vec<u8> = cds.into();
        if !is_acgt(&sequence) {
            return Err(GeneditingError::invalid_target(
                "cds",
                "coding sequence must be non-empty ACGT",
            ));
        }
        if boundaries.is_empty() || boundaries[0] != 0 {
            return Err(GeneditingError::invalid_target(
                "cds",
                "exon boundaries must start with 0",
            ));
        }
        let mut exons = Vec::new();
        for (i, &b) in boundaries.iter().enumerate() {
            let end = boundaries.get(i + 1).copied().unwrap_or(sequence.len());
            if b >= end || end > sequence.len() {
                return Err(GeneditingError::invalid_target(
                    "cds",
                    "exon boundaries must be strictly increasing and in range",
                ));
            }
            exons.push(Exon {
                start: b,
                end,
                ordinal: i + 1,
            });
        }
        Ok(GeneModel {
            sequence,
            exons,
            domains: Vec::new(),
        })
    }

    /// Adds a functional-domain annotation (builder style).
    pub fn with_domain(mut self, name: impl Into<String>, start: usize, end: usize) -> Self {
        self.domains.push(FunctionalDomain {
            name: name.into(),
            start,
            end,
        });
        self
    }

    /// Total coding length.
    pub fn cds_len(&self) -> usize {
        self.sequence.len()
    }
}

/// A request for an NHEJ knockout design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct KnockoutRequest {
    /// The gene to knock out.
    pub gene: GeneModel,
    /// Which nuclease to use.
    pub nuclease: NucleaseId,
    /// Skip the first `lead_in` bp of the CDS — a cut inside the first
    /// few codons can be rescued by a downstream in-frame ATG, so
    /// guides too close to the start codon are demoted (not forbidden).
    pub lead_in: usize,
    /// Avoid the last `tail_fraction` of the CDS — a frameshift in the
    /// final exon often leaves a partly functional product.
    pub tail_fraction: f64,
    /// Keep at most this many ranked guides (`0` = all).
    pub max_results: usize,
}

impl KnockoutRequest {
    /// A request with defaults: SpCas9, 30 bp lead-in, last 20 % of the
    /// CDS avoided.
    pub fn new(gene: GeneModel, nuclease: NucleaseId) -> Self {
        KnockoutRequest {
            gene,
            nuclease,
            lead_in: 30,
            tail_fraction: 0.20,
            max_results: 0,
        }
    }
}

/// One ranked knockout guide.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct KnockoutGuide {
    /// The underlying guide candidate.
    pub guide: GuideCandidate,
    /// 1-based exon the cut falls in.
    pub exon_ordinal: usize,
    /// 0-based cut site in CDS coordinates.
    pub cds_cut_site: usize,
    /// Fractional position of the cut within the CDS in `[0, 1]`.
    pub cds_fraction: f64,
    /// Name of the functional domain the cut disrupts, if any.
    pub domain_hit: Option<String>,
    /// The transparent knockout score in `[0, 1]` — higher is a better
    /// frameshift-knockout guide.
    pub knockout_score: f64,
}

/// The result of a knockout design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct KnockoutStrategy {
    /// Ranked knockout guides, best first.
    pub guides: Vec<KnockoutGuide>,
    /// A one-line rationale for the chosen strategy.
    pub rationale: String,
}

impl KnockoutStrategy {
    /// The single best knockout guide, if any.
    pub fn best(&self) -> Option<&KnockoutGuide> {
        self.guides.first()
    }
}

/// Transparent knockout score for a cut at CDS fraction `frac`.
///
/// - **Position term** — peaks for an early cut just past the lead-in,
///   falls toward the 3′ end, and is heavily demoted inside the
///   tail-fraction (a last-exon frameshift may not knock out).
/// - **Domain term** — a bonus when the cut lands inside an annotated
///   functional domain.
/// - **Efficiency term** — the guide's own on-target / specificity
///   design score carries through.
fn knockout_score(
    frac: f64,
    cds_cut: usize,
    lead_in: usize,
    tail_fraction: f64,
    domain_hit: bool,
    design_score: f64,
) -> f64 {
    // Position term: 1.0 at the lead-in, decaying linearly to ~0.25 at
    // the 3' end; severely penalised once inside the tail fraction.
    let pos = if frac >= 1.0 - tail_fraction {
        0.10
    } else if cds_cut < lead_in {
        0.45 // too close to the start codon — rescuable
    } else {
        1.0 - 0.75 * frac
    };
    let domain = if domain_hit { 0.15 } else { 0.0 };
    // Blend: position dominates, efficiency modulates, domain bonus.
    let raw = 0.55 * pos + 0.30 * design_score + domain;
    raw.clamp(0.0, 1.0)
}

/// Designs an NHEJ knockout strategy for a gene.
///
/// Scans every exon for guides, maps each cut into CDS coordinates,
/// scores it for frameshift-knockout quality and returns the guides
/// ranked best-first.
///
/// # Errors
/// - [`GeneditingError::Invalid`] for a `tail_fraction` outside
///   `[0, 1)`.
/// - [`GeneditingError::NoValidDesign`] if no exon yields a guide.
pub fn design_knockout(req: &KnockoutRequest) -> Result<KnockoutStrategy> {
    if !(0.0..1.0).contains(&req.tail_fraction) {
        return Err(GeneditingError::invalid(
            "tail_fraction",
            "must be in [0, 1)",
        ));
    }
    let cds_len = req.gene.cds_len();
    if cds_len == 0 {
        return Err(GeneditingError::invalid_target("cds", "empty gene model"));
    }
    let mut out: Vec<KnockoutGuide> = Vec::new();

    for exon in &req.gene.exons {
        if exon.is_empty() {
            continue;
        }
        let exon_seq = &req.gene.sequence[exon.start..exon.end];
        let g_req = GuideDesignRequest::new(exon_seq.to_vec(), req.nuclease);
        // A short exon may simply have no PAM site — that is not a
        // failure of the whole design, so swallow a no-design here.
        let guides = match design_guides(&g_req) {
            Ok(g) => g,
            Err(GeneditingError::NoValidDesign { .. }) => continue,
            Err(e) => return Err(e),
        };
        for g in guides {
            // Cut site in CDS coordinates: exon offset + within-exon cut.
            let cds_cut = exon.start + g.cut_site;
            let frac = cds_cut as f64 / cds_len as f64;
            let domain_hit = req
                .gene
                .domains
                .iter()
                .find(|d| cds_cut >= d.start && cds_cut < d.end)
                .map(|d| d.name.clone());
            let score = knockout_score(
                frac,
                cds_cut,
                req.lead_in,
                req.tail_fraction,
                domain_hit.is_some(),
                g.design_score,
            );
            out.push(KnockoutGuide {
                guide: g,
                exon_ordinal: exon.ordinal,
                cds_cut_site: cds_cut,
                cds_fraction: frac,
                domain_hit,
                knockout_score: score,
            });
        }
    }
    if out.is_empty() {
        return Err(GeneditingError::no_valid_design(
            "knockout",
            "no exon contains a PAM-adjacent guide",
        ));
    }
    out.sort_by(|a, b| {
        b.knockout_score
            .partial_cmp(&a.knockout_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if req.max_results > 0 && out.len() > req.max_results {
        out.truncate(req.max_results);
    }
    let best = &out[0];
    let rationale = format!(
        "Best knockout guide cuts in exon {} at CDS position {} ({:.0}% into the \
         coding sequence){}; an early NHEJ frameshift here truncates the protein.",
        best.exon_ordinal,
        best.cds_cut_site,
        best.cds_fraction * 100.0,
        match &best.domain_hit {
            Some(d) => format!(", disrupting the {d} domain"),
            None => String::new(),
        }
    );
    Ok(KnockoutStrategy {
        guides: out,
        rationale,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // A CDS with several SpCas9 PAM sites distributed along it.
    fn cds() -> Vec<u8> {
        let mut s = Vec::new();
        for _ in 0..6 {
            s.extend_from_slice(b"ACGTACGTACGTACGTACGTAGG");
        }
        s
    }

    #[test]
    fn gene_model_from_cds_lays_out_exons() {
        let g = GeneModel::from_cds(cds(), &[0, 46, 92]).unwrap();
        assert_eq!(g.exons.len(), 3);
        assert_eq!(g.exons[0].ordinal, 1);
        assert_eq!(g.exons[0].start, 0);
        assert_eq!(g.exons[1].start, 46);
        assert_eq!(g.exons[2].end, g.cds_len());
    }

    #[test]
    fn rejects_bad_boundaries() {
        assert!(GeneModel::from_cds(cds(), &[]).is_err());
        assert!(GeneModel::from_cds(cds(), &[5, 10]).is_err()); // must start at 0
        assert!(GeneModel::from_cds(cds(), &[0, 10, 5]).is_err()); // not increasing
    }

    #[test]
    fn designs_a_knockout() {
        let g = GeneModel::from_cds(cds(), &[0, 69]).unwrap();
        let req = KnockoutRequest::new(g, NucleaseId::SpCas9);
        let strat = design_knockout(&req).unwrap();
        assert!(!strat.guides.is_empty());
        assert!(strat.best().is_some());
        assert!(strat.rationale.contains("exon"));
    }

    #[test]
    fn early_cut_outscores_late_cut() {
        // Same guide quality, early vs late position.
        let early = knockout_score(0.10, 50, 30, 0.20, false, 0.6);
        let late = knockout_score(0.70, 350, 30, 0.20, false, 0.6);
        assert!(early > late);
    }

    #[test]
    fn tail_cut_is_heavily_demoted() {
        let mid = knockout_score(0.50, 250, 30, 0.20, false, 0.8);
        let tail = knockout_score(0.95, 475, 30, 0.20, false, 0.8);
        assert!(tail < mid);
        assert!(tail < 0.5);
    }

    #[test]
    fn domain_hit_adds_score() {
        let without = knockout_score(0.40, 200, 30, 0.20, false, 0.6);
        let with = knockout_score(0.40, 200, 30, 0.20, true, 0.6);
        assert!(with > without);
    }

    #[test]
    fn domain_annotation_is_detected() {
        // Domain covering CDS 40..120; a cut there should be flagged.
        let g = GeneModel::from_cds(cds(), &[0, 69])
            .unwrap()
            .with_domain("kinase", 40, 120);
        let req = KnockoutRequest::new(g, NucleaseId::SpCas9);
        let strat = design_knockout(&req).unwrap();
        // At least one guide cuts inside the annotated domain window.
        assert!(strat
            .guides
            .iter()
            .any(|g| g.domain_hit.as_deref() == Some("kinase")));
    }

    #[test]
    fn empty_gene_model_is_rejected_via_boundaries() {
        // An all-N "CDS" is rejected at construction.
        assert!(GeneModel::from_cds(b"NNNN".to_vec(), &[0]).is_err());
    }
}
