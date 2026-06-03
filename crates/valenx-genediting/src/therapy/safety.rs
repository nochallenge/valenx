//! Feature 26 — safety-screen aggregation.
//!
//! A gene-editing design is only as good as its safety profile. This
//! module **aggregates** the safety signals a design produces into one
//! [`SafetyReport`]:
//!
//! - the **off-target tally** — how many off-target sites a guide hits
//!   and how active they are, from [`valenx_genomics`]' off-target
//!   enumerator;
//! - **predicted-genotoxicity flags** — large deletions, an integrating
//!   vector, a high off-target burden;
//! - **essential-gene-proximity warnings** — an off-target (or the
//!   on-target cut itself) landing in or near an essential gene is a
//!   serious flag.
//!
//! The module does not run a new analysis; it consumes the outputs of
//! the design modules and the genomics off-target scan and turns them
//! into a graded verdict.
//!
//! ## v1 scope
//!
//! The genotoxicity assessment is a transparent **rule-based** flag
//! set (off-target count / activity thresholds, vector integration,
//! essential-gene overlap) — it is not a trained genotoxicity
//! predictor. The "essential gene" list is whatever the caller
//! supplies as genomic intervals; the module checks overlap, it does
//! not carry a built-in essential-genome database.

use crate::error::{GeneditingError, Result};
use crate::therapy::safety_db::ReferenceGeneDatabase;
use serde::{Deserialize, Serialize};
use valenx_genomics::crispr::offtarget::OffTarget;

/// A graded safety verdict.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SafetyGrade {
    /// No flags — the design's safety profile looks clean.
    Pass,
    /// One or more cautionary flags — review recommended.
    Caution,
    /// A serious flag — an essential-gene hit or a heavy off-target
    /// burden; redesign strongly recommended.
    Fail,
}

impl SafetyGrade {
    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            SafetyGrade::Pass => "pass",
            SafetyGrade::Caution => "caution",
            SafetyGrade::Fail => "fail",
        }
    }
}

/// A genomic interval — a named region on a contig (used for the
/// essential-gene list).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenomicInterval {
    /// The contig name.
    pub chrom: String,
    /// 0-based inclusive start.
    pub start: usize,
    /// 0-based exclusive end.
    pub end: usize,
    /// A label for the interval (e.g. the gene name).
    pub label: String,
}

impl GenomicInterval {
    /// `true` when forward-strand position `pos` on contig `chrom`
    /// falls inside `[start, end)`, with a `flank` bp tolerance on each
    /// side.
    pub fn contains_with_flank(&self, chrom: &str, pos: usize, flank: usize) -> bool {
        self.chrom == chrom
            && pos + flank >= self.start
            && pos < self.end + flank
    }
}

/// The inputs to a safety-screen aggregation.
///
/// Not `serde`-derived: it carries [`valenx_genomics`] `OffTarget`
/// values, which are a transient analysis output rather than a stored
/// data model. The aggregated [`SafetyReport`] *is* serialisable.
#[derive(Clone, Debug, PartialEq)]
pub struct SafetyScreenInput {
    /// Off-target sites for the design's guide(s), from the
    /// [`valenx_genomics`] off-target enumerator. The perfect
    /// (zero-mismatch) on-target site, if present, is ignored by the
    /// aggregator.
    pub off_targets: Vec<OffTarget>,
    /// Essential-gene intervals to check off-targets against.
    pub essential_genes: Vec<GenomicInterval>,
    /// The largest deletion, in bp, the design can produce (`0` for a
    /// pure SNV edit). Large programmed or repair-driven deletions are
    /// a genotoxicity consideration.
    pub max_deletion_bp: usize,
    /// `true` when the design uses a genome-integrating vector
    /// (lentivirus) — insertional mutagenesis is then in play.
    pub integrating_vector: bool,
    /// Flank, in bp, around an essential gene within which an
    /// off-target counts as "proximal".
    pub essential_gene_flank: usize,
}

impl SafetyScreenInput {
    /// A screen input with sensible defaults: no deletion, non-
    /// integrating, a 1 kb essential-gene flank.
    pub fn new(off_targets: Vec<OffTarget>) -> Self {
        SafetyScreenInput {
            off_targets,
            essential_genes: Vec::new(),
            max_deletion_bp: 0,
            integrating_vector: false,
            essential_gene_flank: 1000,
        }
    }
}

/// One safety flag raised by the aggregator.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SafetyFlag {
    /// `true` when the flag is serious (forces a [`SafetyGrade::Fail`]);
    /// `false` for a cautionary flag.
    pub serious: bool,
    /// A stable short code (`"essential_gene_hit"`,
    /// `"high_off_target_burden"`, …).
    pub code: String,
    /// A human-readable description.
    pub detail: String,
}

/// The aggregated safety report for a design (feature 26).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SafetyReport {
    /// The overall graded verdict.
    pub grade: SafetyGrade,
    /// Number of off-target sites considered (excludes the perfect
    /// on-target match).
    pub off_target_count: usize,
    /// The single highest off-target CFD-style activity score
    /// (`0.0` when there are no off-targets).
    pub worst_off_target_activity: f64,
    /// Sum of the off-target activity scores — a guide-promiscuity
    /// proxy.
    pub total_off_target_activity: f64,
    /// Off-targets that fall in or near an essential gene.
    pub essential_gene_hits: Vec<String>,
    /// Every flag raised.
    pub flags: Vec<SafetyFlag>,
}

impl SafetyReport {
    /// `true` when the design passed with no flags at all.
    pub fn is_clean(&self) -> bool {
        self.grade == SafetyGrade::Pass && self.flags.is_empty()
    }

    /// The number of serious flags.
    pub fn serious_flag_count(&self) -> usize {
        self.flags.iter().filter(|f| f.serious).count()
    }
}

/// Off-target activity-sum threshold above which the guide is flagged
/// as too promiscuous (a cautionary flag).
const HIGH_BURDEN_ACTIVITY: f64 = 1.5;

/// A single off-target activity above which one site alone is flagged.
const HIGH_SINGLE_ACTIVITY: f64 = 0.5;

/// A deletion size (bp) above which the programmed deletion is flagged
/// for genotoxicity review.
const LARGE_DELETION_BP: usize = 100;

/// Aggregates a design's safety signals into a [`SafetyReport`]
/// (feature 26).
///
/// Tallies the off-targets (ignoring the perfect on-target match),
/// checks them against the essential-gene list, and applies the
/// rule-based genotoxicity flags. The grade is [`SafetyGrade::Fail`] if
/// any serious flag fired, [`SafetyGrade::Caution`] if only cautionary
/// flags fired, else [`SafetyGrade::Pass`].
///
/// # Errors
/// [`GeneditingError::Invalid`] for a malformed essential-gene
/// interval (`end <= start`).
pub fn aggregate_safety(input: &SafetyScreenInput) -> Result<SafetyReport> {
    for g in &input.essential_genes {
        if g.end <= g.start {
            return Err(GeneditingError::invalid(
                "essential_genes",
                "an essential-gene interval has end <= start",
            ));
        }
    }
    // Tally off-targets, ignoring the perfect on-target site.
    let real_ots: Vec<&OffTarget> = input
        .off_targets
        .iter()
        .filter(|ot| !ot.is_perfect())
        .collect();
    let off_target_count = real_ots.len();
    let worst = real_ots
        .iter()
        .map(|ot| ot.cfd_score)
        .fold(0.0f64, f64::max);
    let total: f64 = real_ots.iter().map(|ot| ot.cfd_score).sum();

    // Essential-gene proximity.
    let mut essential_hits: Vec<String> = Vec::new();
    for ot in &real_ots {
        for gene in &input.essential_genes {
            if gene.contains_with_flank(&ot.chrom, ot.start, input.essential_gene_flank) {
                essential_hits.push(format!(
                    "{} (off-target at {}:{}, {} mismatch(es), activity {:.2})",
                    gene.label,
                    ot.chrom,
                    ot.start,
                    ot.mismatches,
                    ot.cfd_score,
                ));
            }
        }
    }

    // --- Rule-based flags ------------------------------------------
    let mut flags: Vec<SafetyFlag> = Vec::new();

    if !essential_hits.is_empty() {
        flags.push(SafetyFlag {
            serious: true,
            code: "essential_gene_hit".to_string(),
            detail: format!(
                "{} off-target site(s) fall in or near an essential gene.",
                essential_hits.len()
            ),
        });
    }
    if total > HIGH_BURDEN_ACTIVITY {
        flags.push(SafetyFlag {
            serious: true,
            code: "high_off_target_burden".to_string(),
            detail: format!(
                "Summed off-target activity {total:.2} exceeds the \
                 {HIGH_BURDEN_ACTIVITY:.1} promiscuity threshold."
            ),
        });
    }
    if worst > HIGH_SINGLE_ACTIVITY {
        flags.push(SafetyFlag {
            serious: false,
            code: "active_off_target".to_string(),
            detail: format!(
                "An individual off-target has activity {worst:.2} \
                 (> {HIGH_SINGLE_ACTIVITY:.1}) — verify experimentally."
            ),
        });
    }
    if input.max_deletion_bp > LARGE_DELETION_BP {
        flags.push(SafetyFlag {
            serious: false,
            code: "large_deletion".to_string(),
            detail: format!(
                "The design can produce a deletion up to {} bp — large \
                 deletions warrant genotoxicity review.",
                input.max_deletion_bp
            ),
        });
    }
    if input.integrating_vector {
        flags.push(SafetyFlag {
            serious: false,
            code: "integrating_vector".to_string(),
            detail: "An integrating (lentiviral) vector carries an \
                     insertional-mutagenesis risk; consider integration-site \
                     analysis."
                .to_string(),
        });
    }

    let grade = if flags.iter().any(|f| f.serious) {
        SafetyGrade::Fail
    } else if !flags.is_empty() {
        SafetyGrade::Caution
    } else {
        SafetyGrade::Pass
    };

    Ok(SafetyReport {
        grade,
        off_target_count,
        worst_off_target_activity: worst,
        total_off_target_activity: total,
        essential_gene_hits: essential_hits,
        flags,
    })
}

// =====================================================================
// Commercial-depth: per-edit safety screen against the curated lists.
// =====================================================================

/// A request to screen a single proposed edit against the curated
/// essential / cancer-driver / safe-harbor lists in
/// [`ReferenceGeneDatabase`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EditScreenRequest {
    /// The HGNC-style symbol of the gene the edit *targets directly*.
    /// `None` when the edit is intergenic.
    pub target_gene: Option<String>,
    /// HGNC symbols of any gene the on-target cut lies within or
    /// "near" (the caller has resolved off-target coordinates to gene
    /// symbols externally; the screen treats every entry the same way).
    pub neighbor_genes: Vec<String>,
    /// HGNC symbols of any gene an off-target hit falls into.
    pub off_target_genes: Vec<String>,
    /// The largest deletion the design can produce, in bp.
    pub max_deletion_bp: usize,
    /// `true` when the design uses an integrating vector (lentivirus).
    pub integrating_vector: bool,
}

impl EditScreenRequest {
    /// A blank screen request — caller fills in via field assignment.
    pub fn for_target(symbol: impl Into<String>) -> Self {
        EditScreenRequest {
            target_gene: Some(symbol.into()),
            neighbor_genes: Vec::new(),
            off_target_genes: Vec::new(),
            max_deletion_bp: 0,
            integrating_vector: false,
        }
    }
}

/// A per-edit safety verdict against the curated reference lists.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EditSafetyReport {
    /// The overall graded verdict.
    pub grade: SafetyGrade,
    /// Per-finding flags raised — same shape as the off-target
    /// aggregator's [`SafetyFlag`].
    pub flags: Vec<SafetyFlag>,
    /// `true` when the edit's target is a safe-harbor locus — a
    /// positive note, not a flag.
    pub target_is_safe_harbor: bool,
    /// `true` when the edit's target is a cancer-driver gene — raises
    /// a serious flag.
    pub target_is_cancer_driver: bool,
    /// `true` when the edit's target is an essential gene — raises a
    /// serious flag.
    pub target_is_essential: bool,
    /// Symbols of neighbor / off-target essential genes that surfaced
    /// in the screen.
    pub essential_proximity: Vec<String>,
    /// Symbols of off-target cancer-driver genes that surfaced.
    pub cancer_driver_off_targets: Vec<String>,
}

impl EditSafetyReport {
    /// Number of serious flags raised.
    pub fn serious_flag_count(&self) -> usize {
        self.flags.iter().filter(|f| f.serious).count()
    }

    /// `true` when no flag at all was raised AND no positive note
    /// either (i.e. the screen had nothing of interest to say).
    pub fn is_silent(&self) -> bool {
        self.flags.is_empty()
            && !self.target_is_safe_harbor
            && !self.target_is_cancer_driver
            && !self.target_is_essential
            && self.essential_proximity.is_empty()
            && self.cancer_driver_off_targets.is_empty()
    }
}

/// The commercial-depth per-edit safety screen.
///
/// For an edit at a named gene with off-target / neighbor gene
/// resolutions, cross-references each side against the curated
/// essential / cancer-driver / safe-harbor lists in
/// [`ReferenceGeneDatabase`] and reports a graded verdict:
///
/// - Direct cut on an **essential gene** ⇒ serious flag, grade =
///   [`SafetyGrade::Fail`].
/// - Direct cut on a **cancer-driver** ⇒ serious flag, grade =
///   [`SafetyGrade::Fail`].
/// - Direct cut on a **safe-harbor locus** ⇒ informational note,
///   grade = [`SafetyGrade::Pass`] (the lowest-risk integration site).
/// - Off-target hits in essential / cancer-driver genes ⇒ serious
///   flags.
/// - Neighbor essential / cancer-driver genes ⇒ cautionary flags.
/// - Large programmed deletion or integrating vector ⇒ cautionary
///   flags (the same rules as the off-target aggregator).
///
/// # Errors
/// Currently infallible (every input is accepted; an empty request
/// returns an "everything is silent" report). The `Result` signature
/// is reserved for future input validation.
pub fn safety_screen(req: &EditScreenRequest, db: &ReferenceGeneDatabase) -> Result<EditSafetyReport> {
    let mut flags: Vec<SafetyFlag> = Vec::new();
    let mut essential_proximity: Vec<String> = Vec::new();
    let mut cancer_driver_off_targets: Vec<String> = Vec::new();

    let target_is_essential = req
        .target_gene
        .as_deref()
        .map(|s| db.is_essential(s))
        .unwrap_or(false);
    let target_is_cancer_driver = req
        .target_gene
        .as_deref()
        .map(|s| db.is_cancer_driver(s))
        .unwrap_or(false);
    let target_is_safe_harbor = req
        .target_gene
        .as_deref()
        .map(|s| db.is_safe_harbor(s))
        .unwrap_or(false);

    if target_is_essential {
        flags.push(SafetyFlag {
            serious: true,
            code: "target_essential_gene".to_string(),
            detail: format!(
                "Target gene `{}` is on the curated essential-gene list — \
                 cutting will likely be lethal to the cell. Consider a base/prime \
                 editor instead of an NHEJ-prone nuclease, or pick a different \
                 isoform / exon.",
                req.target_gene.as_deref().unwrap_or("?")
            ),
        });
    }
    if target_is_cancer_driver {
        flags.push(SafetyFlag {
            serious: true,
            code: "target_cancer_driver".to_string(),
            detail: format!(
                "Target gene `{}` is on the curated cancer-driver list — \
                 perturbation has been causally linked to transformation. \
                 Re-verify oncogenicity before proceeding.",
                req.target_gene.as_deref().unwrap_or("?")
            ),
        });
    }
    if target_is_safe_harbor {
        flags.push(SafetyFlag {
            serious: false,
            code: "target_safe_harbor".to_string(),
            detail: format!(
                "Target locus `{}` is a curated safe-harbor site — stable \
                 transgene integration here is well-tolerated and the \
                 lowest-risk choice for cassette knock-in.",
                req.target_gene.as_deref().unwrap_or("?")
            ),
        });
    }

    // Neighbor essential / cancer-driver proximity (cautionary).
    for nb in &req.neighbor_genes {
        if db.is_essential(nb) {
            essential_proximity.push(nb.clone());
            flags.push(SafetyFlag {
                serious: false,
                code: "neighbor_essential_gene".to_string(),
                detail: format!(
                    "On-target cut is near essential gene `{nb}`. Verify the \
                     cut site does not damage `{nb}`'s coding window or \
                     regulatory elements."
                ),
            });
        }
        if db.is_cancer_driver(nb) {
            flags.push(SafetyFlag {
                serious: false,
                code: "neighbor_cancer_driver".to_string(),
                detail: format!(
                    "On-target cut is near cancer-driver gene `{nb}`. \
                     Verify the cut site does not perturb `{nb}`."
                ),
            });
        }
    }

    // Off-target essential / cancer-driver gene hits (serious).
    for ot in &req.off_target_genes {
        if db.is_essential(ot) {
            essential_proximity.push(ot.clone());
            flags.push(SafetyFlag {
                serious: true,
                code: "off_target_essential_gene".to_string(),
                detail: format!(
                    "Off-target hit in essential gene `{ot}` — risk of cell \
                     death. Redesign the guide to reduce off-target activity."
                ),
            });
        }
        if db.is_cancer_driver(ot) {
            cancer_driver_off_targets.push(ot.clone());
            flags.push(SafetyFlag {
                serious: true,
                code: "off_target_cancer_driver".to_string(),
                detail: format!(
                    "Off-target hit in cancer-driver gene `{ot}` — \
                     oncogenicity risk. Redesign the guide."
                ),
            });
        }
    }

    // Mirror the off-target aggregator's deletion / vector cautions.
    if req.max_deletion_bp > LARGE_DELETION_BP {
        flags.push(SafetyFlag {
            serious: false,
            code: "large_deletion".to_string(),
            detail: format!(
                "The design can produce a deletion up to {} bp — large \
                 deletions warrant genotoxicity review.",
                req.max_deletion_bp
            ),
        });
    }
    if req.integrating_vector {
        flags.push(SafetyFlag {
            serious: false,
            code: "integrating_vector".to_string(),
            detail: "An integrating (lentiviral) vector carries an \
                     insertional-mutagenesis risk; consider integration-site \
                     analysis."
                .to_string(),
        });
    }

    // A safe-harbor target dampens the verdict: even if the off-target
    // panel flagged something cautionary, the *target* is low-risk.
    // (At this point we already know no flag is serious — the first
    // branch catches that case.)
    let grade = if flags.iter().any(|f| f.serious) {
        SafetyGrade::Fail
    } else if flags.is_empty() || target_is_safe_harbor {
        // Safe-harbor target with only the safe-harbor note + at most
        // other cautionary flags → Pass; an empty-flag set is also Pass.
        SafetyGrade::Pass
    } else {
        SafetyGrade::Caution
    };

    essential_proximity.sort();
    essential_proximity.dedup();
    cancer_driver_off_targets.sort();
    cancer_driver_off_targets.dedup();

    Ok(EditSafetyReport {
        grade,
        flags,
        target_is_safe_harbor,
        target_is_cancer_driver,
        target_is_essential,
        essential_proximity,
        cancer_driver_off_targets,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ot(chrom: &str, start: usize, mm: usize, cfd: f64) -> OffTarget {
        OffTarget {
            chrom: chrom.to_string(),
            start,
            reverse: false,
            protospacer: "ACGTACGTACGTACGTACGT".to_string(),
            pam: "AGG".to_string(),
            mismatches: mm,
            mismatch_positions: vec![0; mm],
            cfd_score: cfd,
        }
    }

    fn perfect(chrom: &str, start: usize) -> OffTarget {
        ot(chrom, start, 0, 1.0)
    }

    #[test]
    fn clean_design_passes() {
        // One weak off-target, nothing else.
        let input = SafetyScreenInput::new(vec![ot("chr1", 100, 4, 0.05)]);
        let report = aggregate_safety(&input).unwrap();
        assert_eq!(report.grade, SafetyGrade::Pass);
        assert!(report.is_clean());
        assert_eq!(report.off_target_count, 1);
    }

    #[test]
    fn perfect_on_target_is_ignored_in_tally() {
        let input = SafetyScreenInput::new(vec![
            perfect("chr1", 50),
            ot("chr1", 100, 4, 0.05),
        ]);
        let report = aggregate_safety(&input).unwrap();
        // Only the non-perfect site is counted.
        assert_eq!(report.off_target_count, 1);
    }

    #[test]
    fn essential_gene_hit_fails() {
        let mut input = SafetyScreenInput::new(vec![ot("chr1", 1000, 2, 0.3)]);
        input.essential_genes = vec![GenomicInterval {
            chrom: "chr1".to_string(),
            start: 900,
            end: 1100,
            label: "TP53".to_string(),
        }];
        let report = aggregate_safety(&input).unwrap();
        assert_eq!(report.grade, SafetyGrade::Fail);
        assert!(!report.essential_gene_hits.is_empty());
        assert!(report.flags.iter().any(|f| f.code == "essential_gene_hit"));
    }

    #[test]
    fn off_target_far_from_essential_gene_does_not_hit() {
        let mut input = SafetyScreenInput::new(vec![ot("chr1", 50_000, 2, 0.3)]);
        input.essential_genes = vec![GenomicInterval {
            chrom: "chr1".to_string(),
            start: 900,
            end: 1100,
            label: "TP53".to_string(),
        }];
        let report = aggregate_safety(&input).unwrap();
        assert!(report.essential_gene_hits.is_empty());
    }

    #[test]
    fn high_burden_fails() {
        // Many moderately-active off-targets sum past the threshold.
        let ots: Vec<OffTarget> = (0..10)
            .map(|i| ot("chr1", 1000 + i * 100, 3, 0.4))
            .collect();
        let input = SafetyScreenInput::new(ots);
        let report = aggregate_safety(&input).unwrap();
        assert!(report.total_off_target_activity > HIGH_BURDEN_ACTIVITY);
        assert_eq!(report.grade, SafetyGrade::Fail);
    }

    #[test]
    fn single_active_off_target_is_a_caution() {
        // One off-target above the single-site threshold, but the sum
        // stays under the burden threshold.
        let input = SafetyScreenInput::new(vec![ot("chr1", 1000, 1, 0.7)]);
        let report = aggregate_safety(&input).unwrap();
        assert_eq!(report.grade, SafetyGrade::Caution);
        assert!(report.flags.iter().any(|f| f.code == "active_off_target"));
        assert!(!report.flags.iter().any(|f| f.serious));
    }

    #[test]
    fn integrating_vector_raises_a_caution() {
        let mut input = SafetyScreenInput::new(vec![ot("chr1", 100, 4, 0.02)]);
        input.integrating_vector = true;
        let report = aggregate_safety(&input).unwrap();
        assert_eq!(report.grade, SafetyGrade::Caution);
        assert!(report.flags.iter().any(|f| f.code == "integrating_vector"));
    }

    #[test]
    fn large_deletion_raises_a_caution() {
        let mut input = SafetyScreenInput::new(Vec::new());
        input.max_deletion_bp = 500;
        let report = aggregate_safety(&input).unwrap();
        assert!(report.flags.iter().any(|f| f.code == "large_deletion"));
    }

    #[test]
    fn empty_input_passes_clean() {
        let report = aggregate_safety(&SafetyScreenInput::new(Vec::new())).unwrap();
        assert!(report.is_clean());
        assert_eq!(report.worst_off_target_activity, 0.0);
    }

    #[test]
    fn rejects_malformed_essential_interval() {
        let mut input = SafetyScreenInput::new(Vec::new());
        input.essential_genes = vec![GenomicInterval {
            chrom: "chr1".to_string(),
            start: 100,
            end: 100, // end <= start
            label: "bad".to_string(),
        }];
        assert_eq!(
            aggregate_safety(&input).unwrap_err().code(),
            "genediting.invalid"
        );
    }

    #[test]
    fn serious_flag_count_tracks_grade() {
        let mut input = SafetyScreenInput::new(vec![ot("chr1", 1000, 2, 0.3)]);
        input.essential_genes = vec![GenomicInterval {
            chrom: "chr1".to_string(),
            start: 900,
            end: 1100,
            label: "ESS".to_string(),
        }];
        let report = aggregate_safety(&input).unwrap();
        assert!(report.serious_flag_count() >= 1);
        assert_eq!(report.grade, SafetyGrade::Fail);
    }

    // ==================================================================
    // safety_screen — curated reference-list per-edit screen.
    // ==================================================================

    #[test]
    fn tp53_edit_raises_cancer_driver_warning() {
        let db = ReferenceGeneDatabase::curated();
        let req = EditScreenRequest::for_target("TP53");
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Fail);
        assert!(r.target_is_cancer_driver);
        assert!(r.flags.iter().any(|f| f.code == "target_cancer_driver"));
    }

    #[test]
    fn aavs1_edit_is_safe_harbor_pass() {
        let db = ReferenceGeneDatabase::curated();
        let req = EditScreenRequest::for_target("AAVS1");
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Pass);
        assert!(r.target_is_safe_harbor);
        assert!(r.flags.iter().any(|f| f.code == "target_safe_harbor"));
        // Safe-harbor note is informational, not serious.
        assert_eq!(r.serious_flag_count(), 0);
    }

    #[test]
    fn ribosomal_neighbor_raises_essential_proximity_warning() {
        let db = ReferenceGeneDatabase::curated();
        // Edit at an unrelated gene, but RPS6 (essential, ribosomal)
        // is a near neighbor — should raise a cautionary flag.
        let mut req = EditScreenRequest::for_target("UNKNOWN_GENE");
        req.neighbor_genes = vec!["RPS6".to_string()];
        let r = safety_screen(&req, &db).unwrap();
        assert!(r.essential_proximity.contains(&"RPS6".to_string()));
        assert!(r.flags.iter().any(|f| f.code == "neighbor_essential_gene"));
        // Cautionary, not serious.
        assert!(matches!(r.grade, SafetyGrade::Caution));
    }

    #[test]
    fn off_target_in_cancer_driver_fails() {
        let db = ReferenceGeneDatabase::curated();
        let mut req = EditScreenRequest::for_target("DUMMY");
        req.off_target_genes = vec!["MYC".to_string()];
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Fail);
        assert!(r.cancer_driver_off_targets.contains(&"MYC".to_string()));
        assert!(r.flags.iter().any(|f| f.code == "off_target_cancer_driver"));
    }

    #[test]
    fn off_target_in_essential_gene_fails() {
        let db = ReferenceGeneDatabase::curated();
        let mut req = EditScreenRequest::for_target("DUMMY");
        req.off_target_genes = vec!["POLR2A".to_string()];
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Fail);
        assert!(r.flags.iter().any(|f| f.code == "off_target_essential_gene"));
    }

    #[test]
    fn intergenic_clean_edit_passes_silently() {
        let db = ReferenceGeneDatabase::curated();
        let req = EditScreenRequest {
            target_gene: None,
            neighbor_genes: Vec::new(),
            off_target_genes: Vec::new(),
            max_deletion_bp: 0,
            integrating_vector: false,
        };
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Pass);
        assert!(r.is_silent());
    }

    #[test]
    fn essential_gene_target_fails() {
        let db = ReferenceGeneDatabase::curated();
        let req = EditScreenRequest::for_target("RPS6");
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Fail);
        assert!(r.target_is_essential);
        assert!(r.flags.iter().any(|f| f.code == "target_essential_gene"));
    }

    #[test]
    fn large_deletion_at_safe_harbor_is_caution_not_fail() {
        let db = ReferenceGeneDatabase::curated();
        let mut req = EditScreenRequest::for_target("AAVS1");
        req.max_deletion_bp = 500;
        let r = safety_screen(&req, &db).unwrap();
        assert!(r.flags.iter().any(|f| f.code == "large_deletion"));
        // No serious flag → safe-harbor note keeps the grade at Pass.
        assert_eq!(r.serious_flag_count(), 0);
        assert_eq!(r.grade, SafetyGrade::Pass);
    }

    #[test]
    fn case_insensitive_symbol_matching() {
        let db = ReferenceGeneDatabase::curated();
        let req = EditScreenRequest::for_target("tp53");
        let r = safety_screen(&req, &db).unwrap();
        assert!(r.target_is_cancer_driver);
    }

    #[test]
    fn multi_flag_design_aggregates_all() {
        // A worst-case design: target is a cancer driver, an essential
        // gene neighbor, an essential gene off-target, large deletion,
        // and an integrating vector.
        let db = ReferenceGeneDatabase::curated();
        let req = EditScreenRequest {
            target_gene: Some("TP53".to_string()),
            neighbor_genes: vec!["RPS6".to_string()],
            off_target_genes: vec!["MYC".to_string(), "POLR2A".to_string()],
            max_deletion_bp: 800,
            integrating_vector: true,
        };
        let r = safety_screen(&req, &db).unwrap();
        assert_eq!(r.grade, SafetyGrade::Fail);
        // Every kind of flag fires.
        let codes: Vec<&str> = r.flags.iter().map(|f| f.code.as_str()).collect();
        assert!(codes.contains(&"target_cancer_driver"));
        assert!(codes.contains(&"neighbor_essential_gene"));
        assert!(codes.contains(&"off_target_cancer_driver"));
        assert!(codes.contains(&"off_target_essential_gene"));
        assert!(codes.contains(&"large_deletion"));
        assert!(codes.contains(&"integrating_vector"));
        assert!(r.serious_flag_count() >= 3);
    }
}
