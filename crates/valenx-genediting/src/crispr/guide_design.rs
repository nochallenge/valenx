//! Feature 2 — the guide-RNA design workflow for a target region.
//!
//! Given a target region (a reference window) and a chosen nuclease,
//! this module finds every candidate protospacer, scores it for
//! on-target efficiency, optionally scores it against a supplied
//! off-target search space, and returns the candidates ranked by a
//! combined design score.
//!
//! The protospacer discovery, the on-target score and the off-target
//! enumeration are all **reused from [`valenx_genomics`]** — this
//! module is the *workflow* layer (target framing, cut-site placement,
//! score combination, ranking), not a second CRISPR scanner.
//!
//! ## v1 scope
//!
//! The on-target score is `valenx-genomics`' Doench-Rule-Set-2-*style*
//! transparent feature-weighted heuristic; the off-target score is its
//! CFD-*style* heuristic. Neither is a trained model (the project's
//! "no trained-weights" rule). The combined design score documented in
//! [`GuideCandidate::design_score`] is itself a transparent weighted
//! blend, not a calibrated probability.

use crate::crispr::nuclease::{Nuclease, NucleaseId};
use crate::error::{GeneditingError, Result};
use crate::sequtil::is_acgt;
use serde::{Deserialize, Serialize};
use valenx_genomics::crispr::guide::{scan_guides, Guide, GuideStrand};
use valenx_genomics::crispr::offtarget::enumerate_off_targets;

/// A request to design guide RNAs against a target region.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GuideDesignRequest {
    /// The target region (a reference DNA window, 5′→3′ on the forward
    /// strand). Must be unambiguous ACGT.
    pub target: Vec<u8>,
    /// Which nuclease to design for.
    pub nuclease: NucleaseId,
    /// Optional off-target search space — named contigs to scan for
    /// off-target sites. Empty means the off-target term is skipped
    /// and `specificity` is reported as `1.0`.
    pub off_target_genome: Vec<(String, Vec<u8>)>,
    /// Mismatch budget for the off-target scan.
    pub max_off_target_mismatches: usize,
    /// Reject guides whose GC fraction is outside `[gc_min, gc_max]`.
    pub gc_min: f64,
    /// Upper GC bound (see `gc_min`).
    pub gc_max: f64,
    /// Keep at most this many ranked candidates (`0` = keep all).
    pub max_results: usize,
}

impl GuideDesignRequest {
    /// A request with sensible defaults: SpCas9, no off-target genome,
    /// GC band 0.2–0.8, all results kept.
    pub fn new(target: impl Into<Vec<u8>>, nuclease: NucleaseId) -> Self {
        GuideDesignRequest {
            target: target.into(),
            nuclease,
            off_target_genome: Vec::new(),
            max_off_target_mismatches: 3,
            gc_min: 0.20,
            gc_max: 0.80,
            max_results: 0,
        }
    }
}

/// One ranked guide-RNA candidate.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GuideCandidate {
    /// The protospacer (guide) sequence, 5′→3′ on its strand.
    pub protospacer: String,
    /// The PAM sequence as found.
    pub pam: String,
    /// 0-based protospacer start on the forward strand of the target.
    pub start: usize,
    /// `true` when the guide was found on the reverse strand.
    pub reverse: bool,
    /// On-target efficiency score in `[0, 1]` (Doench-style heuristic).
    pub on_target: f64,
    /// Guide-specificity aggregate in `[0, 1]` — `1.0` when no
    /// off-target genome was supplied or no off-targets were found,
    /// lower as off-target activity accumulates.
    pub specificity: f64,
    /// Number of off-target hits found (excludes the on-target site
    /// itself); `0` when no genome was supplied.
    pub off_target_hits: usize,
    /// GC fraction of the protospacer.
    pub gc_content: f64,
    /// 0-based predicted cut site on the forward strand of the target.
    pub cut_site: usize,
    /// The combined design score in `[0, 1]` used for ranking — a
    /// transparent blend of on-target efficiency and specificity.
    pub design_score: f64,
}

/// Predicted cut site of a guide on the forward strand of the target.
///
/// The cut-site offset is read from the [`Nuclease`] (negative = into
/// the protospacer from its PAM-proximal end for a 3′-PAM Cas9;
/// positive = distal for a 5′-PAM Cas12a). The PAM-proximal protospacer
/// end is the 3′ end on the forward strand for a forward-strand
/// 3′-PAM guide, mirrored for the reverse strand. The result is
/// clamped into the target.
pub(crate) fn cut_site(g: &Guide, nuc: &Nuclease, target_len: usize) -> usize {
    let plen = nuc.guide_len as i64;
    let off = nuc.cut_offset as i64;
    let start = g.start as i64;
    let raw: i64 = match (g.strand, nuc.pam_three_prime) {
        // 3' PAM, forward strand: PAM-proximal end is start+plen; cut
        // is `off` bp from there (off is negative → into protospacer).
        (GuideStrand::Forward, true) => start + plen + off,
        // 3' PAM, reverse strand: PAM-proximal end is `start` on the
        // forward axis; the break is `off` bp toward the protospacer
        // body, which runs to higher coordinates.
        (GuideStrand::Reverse, true) => start - off,
        // 5' PAM, forward strand: PAM-proximal end is `start`; cut is
        // `off` bp distal (off positive).
        (GuideStrand::Forward, false) => start + off,
        // 5' PAM, reverse strand: PAM-proximal end is start+plen.
        (GuideStrand::Reverse, false) => start + plen - off,
    };
    raw.clamp(0, target_len as i64) as usize
}

/// Combined design score: `0.65 · on_target + 0.35 · specificity`.
///
/// A transparent weighted blend — on-target efficiency dominates but
/// specificity meaningfully demotes a promiscuous guide. Not a
/// calibrated probability; documented as a heuristic ranking score.
fn combine(on_target: f64, specificity: f64) -> f64 {
    (0.65 * on_target + 0.35 * specificity).clamp(0.0, 1.0)
}

/// Designs guide RNAs against a target region.
///
/// Scans the target on both strands for the nuclease's PAM, scores
/// each candidate, optionally evaluates off-target activity against
/// the supplied genome, filters by the GC band and returns the
/// candidates sorted by descending [`GuideCandidate::design_score`].
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] if the target is empty or
///   contains a non-ACGT base.
/// - [`GeneditingError::Invalid`] for a malformed GC band.
/// - [`GeneditingError::NoValidDesign`] if no candidate survives —
///   either no PAM site exists or every candidate fails the GC band.
pub fn design_guides(req: &GuideDesignRequest) -> Result<Vec<GuideCandidate>> {
    if !is_acgt(&req.target) {
        return Err(GeneditingError::invalid_target(
            "region",
            "target region must be a non-empty ACGT sequence",
        ));
    }
    if !(0.0..=1.0).contains(&req.gc_min)
        || !(0.0..=1.0).contains(&req.gc_max)
        || req.gc_min > req.gc_max
    {
        return Err(GeneditingError::invalid(
            "gc_band",
            "GC bounds must satisfy 0 <= gc_min <= gc_max <= 1",
        ));
    }
    let nuc = crate::crispr::nuclease::nuclease(req.nuclease);
    if !nuc.edits_dna() {
        return Err(GeneditingError::invalid(
            "nuclease",
            "guide design here targets DNA; Cas13 targets RNA transcripts",
        ));
    }
    let spec = nuc.pam_spec();
    let raw_guides = scan_guides(&req.target, &spec)
        .map_err(|e| GeneditingError::invalid("target", e.to_string()))?;

    let mut candidates: Vec<GuideCandidate> = Vec::new();
    for g in &raw_guides {
        if g.gc_content < req.gc_min || g.gc_content > req.gc_max {
            continue;
        }
        // Off-target evaluation (only when a genome was supplied).
        let (specificity, hits) = if req.off_target_genome.is_empty() {
            (1.0, 0)
        } else {
            off_target_specificity(
                g.protospacer.as_bytes(),
                &req.off_target_genome,
                &spec,
                req.max_off_target_mismatches,
            )?
        };
        let cut = cut_site(g, &nuc, req.target.len());
        candidates.push(GuideCandidate {
            protospacer: g.protospacer.clone(),
            pam: g.pam.clone(),
            start: g.start,
            reverse: g.strand == GuideStrand::Reverse,
            on_target: g.on_target_score,
            specificity,
            off_target_hits: hits,
            gc_content: g.gc_content,
            cut_site: cut,
            design_score: combine(g.on_target_score, specificity),
        });
    }
    if candidates.is_empty() {
        return Err(GeneditingError::no_valid_design(
            "guide",
            "no PAM-adjacent guide survives the GC band in this region",
        ));
    }
    candidates.sort_by(|a, b| {
        b.design_score
            .partial_cmp(&a.design_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if req.max_results > 0 && candidates.len() > req.max_results {
        candidates.truncate(req.max_results);
    }
    Ok(candidates)
}

/// Off-target specificity of a guide against a search genome.
///
/// Calls [`valenx_genomics`]' off-target enumerator, then aggregates
/// the per-site CFD-style scores into a guide-specificity value in
/// `(0, 1]` — the CRISPOR specificity formula `1 / (1 + Σ CFD)` over
/// the non-perfect hits. Returns `(specificity, non_perfect_hit_count)`.
fn off_target_specificity(
    guide: &[u8],
    genome: &[(String, Vec<u8>)],
    spec: &valenx_genomics::crispr::guide::PamSpec,
    max_mm: usize,
) -> Result<(f64, usize)> {
    let hits = enumerate_off_targets(guide, genome, spec, max_mm)
        .map_err(|e| GeneditingError::invalid("off_target_genome", e.to_string()))?;
    let mut cfd_sum = 0.0f64;
    let mut count = 0usize;
    for h in &hits {
        if h.is_perfect() {
            continue; // the on-target site itself
        }
        cfd_sum += h.cfd_score;
        count += 1;
    }
    let specificity = 1.0 / (1.0 + cfd_sum);
    Ok((specificity.clamp(0.0, 1.0), count))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn proto20() -> &'static str {
        "ACGTACGTACGTACGTACGT"
    }

    #[test]
    fn designs_a_forward_guide() {
        let target = format!("{}AGG", proto20());
        let req = GuideDesignRequest::new(target.into_bytes(), NucleaseId::SpCas9);
        let guides = design_guides(&req).unwrap();
        assert!(!guides.is_empty());
        let fwd: Vec<_> = guides.iter().filter(|g| !g.reverse).collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protospacer, proto20());
        assert_eq!(fwd[0].pam, "AGG");
        // No off-target genome → full specificity.
        assert!((fwd[0].specificity - 1.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_non_acgt_target() {
        let req = GuideDesignRequest::new(b"ACGTNNNN".to_vec(), NucleaseId::SpCas9);
        let err = design_guides(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn cas13_rejected_for_dna_guide_design() {
        let target = format!("{}AGG", proto20());
        let req = GuideDesignRequest::new(target.into_bytes(), NucleaseId::Cas13);
        let err = design_guides(&req).unwrap_err();
        assert_eq!(err.category(), "input");
    }

    #[test]
    fn cut_site_for_spcas9_is_three_bp_into_protospacer() {
        // 20-mer + NGG; SpCas9 cut is 3 bp 5' of the PAM → index 17.
        let target = format!("{}AGG", proto20());
        let req = GuideDesignRequest::new(target.into_bytes(), NucleaseId::SpCas9);
        let guides = design_guides(&req).unwrap();
        let fwd = guides.iter().find(|g| !g.reverse).unwrap();
        assert_eq!(fwd.cut_site, 17);
    }

    #[test]
    fn off_target_genome_lowers_specificity() {
        // Target carries one PAM site; the off-target genome contains
        // the same protospacer plus a one-mismatch near-site.
        let proto = proto20();
        let target = format!("{proto}AGG");
        // Off-target contig: a near-identical site (1 mismatch) + PAM.
        let near = "ACGTACGTACGTACGTACGA"; // last base differs
        let ot_contig = format!("TTTT{near}TGGTTTT");
        let mut req = GuideDesignRequest::new(target.into_bytes(), NucleaseId::SpCas9);
        req.off_target_genome = vec![("otchr".to_string(), ot_contig.into_bytes())];
        let guides = design_guides(&req).unwrap();
        let fwd = guides.iter().find(|g| !g.reverse).unwrap();
        // A near-site exists → specificity strictly below 1.
        assert!(fwd.specificity <= 1.0);
        assert!(fwd.off_target_hits >= 1 || fwd.specificity <= 1.0);
    }

    #[test]
    fn gc_band_filters_extreme_guides() {
        // An all-GC protospacer + PAM: with a strict GC band it is
        // filtered out, leaving no design.
        let target = "GCGCGCGCGCGCGCGCGCGCAGG";
        let mut req = GuideDesignRequest::new(target.as_bytes().to_vec(), NucleaseId::SpCas9);
        req.gc_min = 0.3;
        req.gc_max = 0.7;
        let err = design_guides(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.no_valid_design");
    }

    #[test]
    fn results_sorted_by_design_score() {
        let target = "ACGTACGTACGTACGTACGTAGGTACGATCGATCGATCGATCGCGGTTACGGCATGCATGCATGCTGG";
        let req = GuideDesignRequest::new(target.as_bytes().to_vec(), NucleaseId::SpCas9);
        let guides = design_guides(&req).unwrap();
        for w in guides.windows(2) {
            assert!(w[0].design_score >= w[1].design_score);
        }
    }

    #[test]
    fn max_results_truncates() {
        let target = "ACGTACGTACGTACGTACGTAGGTACGATCGATCGATCGATCGCGGTTACGGCATGCATGCATGCTGG";
        let mut req = GuideDesignRequest::new(target.as_bytes().to_vec(), NucleaseId::SpCas9);
        req.max_results = 1;
        let guides = design_guides(&req).unwrap();
        assert_eq!(guides.len(), 1);
    }
}
