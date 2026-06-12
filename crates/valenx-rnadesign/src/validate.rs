//! Group D — in-silico validation (features 15–18).
//!
//! Before a design is handed to synthesis it is **validated in
//! silico**: every constraint is checked, the design is re-folded and
//! compared to its target, the Boltzmann ensemble is characterised, and
//! a battery of robustness sanity checks is run. The output is a
//! [`ValidationReport`] — a pass / warn / fail verdict per constraint,
//! plus the structural metrics.
//!
//! - [`fold_back_validate`] (feature 15) — folds the design with the
//!   MFE folder and compares the result to the target (base-pair
//!   distance, % structure match).
//! - [`ensemble_validate`] (feature 16) — the partition function:
//!   base-pair probabilities, the probability of the target structure,
//!   the ensemble defect.
//! - [`robustness_checks`] (feature 17) — mutational robustness, a
//!   melting / temperature sanity check, a co-transcriptional-folding
//!   sanity check.
//! - [`validate_design`] (feature 18) — runs everything and assembles
//!   the [`ValidationReport`].
//!
//! ## v1 scope — honest framing, READ THIS
//!
//! **A passing [`ValidationReport`] is not a guarantee.** Every metric
//! in this module is computed from a nearest-neighbor secondary-
//! structure energy model (`valenx-rnastruct`'s Turner-2004
//! parameters — a faithful representative subset). A design this module
//! reports as folding to its target, with a high target-structure
//! probability and a low ensemble defect, is a **strong in-silico
//! prediction** — a well-supported hypothesis — and nothing more.
//! Physical RNA must still be synthesised and lab-validated: folded
//! experimentally, assayed for function, tested in cells. Nothing here
//! is named "verified" or "guaranteed correct"; the verdict is named
//! [`ValidationVerdict::Pass`] in the *in-silico-prediction* sense
//! only. The co-transcriptional check is a v1 sanity heuristic (a
//! 5′-window fold progression), not a kinetic folding simulation.

use crate::design::{DesignKind, RnaDesign};
use crate::error::{Result, RnaDesignError};
use crate::goal::DesignConstraints;
use crate::optimize::{
    ensemble_defect, forbidden_motif_scan, gc_content, repeat_scan, synthesizability_scan,
};
use serde::{Deserialize, Serialize};
use valenx_rnastruct::{
    base_pair_distance, melting_curve, mfe, partition_function, structure_energy, RnaSeq, Structure,
};

/// The pass / warn / fail verdict of one check.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ValidationVerdict {
    /// The check passed (in the in-silico-prediction sense).
    Pass,
    /// The check raised a non-fatal concern.
    Warn,
    /// The check failed — the design violates a hard requirement.
    Fail,
}

impl ValidationVerdict {
    /// A short human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            ValidationVerdict::Pass => "pass",
            ValidationVerdict::Warn => "warn",
            ValidationVerdict::Fail => "fail",
        }
    }

    /// Combines two verdicts, keeping the *worst* (`Fail` > `Warn` >
    /// `Pass`) — the rule a [`ValidationReport`] uses to aggregate.
    pub fn worst(self, other: ValidationVerdict) -> ValidationVerdict {
        use ValidationVerdict::*;
        match (self, other) {
            (Fail, _) | (_, Fail) => Fail,
            (Warn, _) | (_, Warn) => Warn,
            _ => Pass,
        }
    }
}

/// The result of checking one constraint (feature 18).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConstraintCheck {
    /// What was checked (`"GC content"`, `"length"`,
    /// `"forbidden motifs"`, …).
    pub name: String,
    /// The verdict.
    pub verdict: ValidationVerdict,
    /// A human-readable explanation of the result.
    pub detail: String,
}

/// The result of folding the design back and comparing it to the target
/// (feature 15).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FoldBackResult {
    /// The dot-bracket of the design's MFE fold.
    pub achieved_dot_bracket: String,
    /// The dot-bracket of the target structure.
    pub target_dot_bracket: String,
    /// Base-pair distance between the achieved fold and the target.
    pub base_pair_distance: usize,
    /// Percentage of the target's base pairs the design recovers,
    /// `[0, 100]`.
    pub structure_match_percent: f64,
    /// The design's minimum free energy (kcal/mol).
    pub mfe_energy: f64,
}

impl FoldBackResult {
    /// `true` when the design folds *exactly* to the target.
    pub fn is_exact(&self) -> bool {
        self.base_pair_distance == 0
    }
}

/// The result of characterising the design's Boltzmann ensemble
/// (feature 16).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EnsembleResult {
    /// The ensemble free energy `G = −RT·ln Q` (kcal/mol).
    pub ensemble_free_energy: f64,
    /// The Boltzmann probability of the *target* structure, `[0, 1]`
    /// (`exp(−(E_target − G)/RT)`).
    pub target_probability: f64,
    /// The ensemble defect — the expected number of incorrectly-paired
    /// bases relative to the target.
    pub ensemble_defect: f64,
    /// The normalised ensemble defect, `[0, 1]`.
    pub normalized_defect: f64,
}

/// The result of the robustness sanity checks (feature 17).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RobustnessResult {
    /// Mutational robustness — the fraction of single point mutations
    /// whose MFE fold stays close to the original, `[0, 1]`.
    pub mutational_robustness: f64,
    /// The estimated melting temperature (°C), or `None` if the design
    /// has no cooperative transition.
    pub melting_temperature_c: Option<f64>,
    /// `true` when the co-transcriptional-folding sanity check found no
    /// red flag (no strong premature 5′ structure that would trap a
    /// mis-fold during transcription).
    pub cotranscriptional_ok: bool,
    /// Human-readable notes from the robustness checks.
    pub notes: Vec<String>,
}

/// The full in-silico validation report (feature 18).
///
/// # Honest framing
///
/// A [`ValidationReport`] with [`verdict`](Self::verdict) `Pass` means
/// the design passed every in-silico check **under a nearest-neighbor
/// energy model** — it is a strong predicted candidate, not a verified
/// or guaranteed-correct molecule. Physical RNA must still be
/// synthesised and validated in the wet lab. See the module
/// documentation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ValidationReport {
    /// The aggregate verdict — the worst of every constraint check.
    pub verdict: ValidationVerdict,
    /// The per-constraint check results.
    pub constraint_checks: Vec<ConstraintCheck>,
    /// The fold-back result (feature 15) — `None` for a design with no
    /// structural target (a bare coding mRNA).
    pub fold_back: Option<FoldBackResult>,
    /// The ensemble characterisation (feature 16) — `None` for a design
    /// with no structural target.
    pub ensemble: Option<EnsembleResult>,
    /// The robustness sanity checks (feature 17).
    pub robustness: Option<RobustnessResult>,
    /// Human-readable notes, including the in-silico-prediction
    /// disclaimer.
    pub notes: Vec<String>,
}

impl ValidationReport {
    /// `true` when the aggregate verdict is [`ValidationVerdict::Pass`].
    pub fn passed(&self) -> bool {
        self.verdict == ValidationVerdict::Pass
    }

    /// The number of constraint checks with a given verdict.
    pub fn count_with(&self, verdict: ValidationVerdict) -> usize {
        self.constraint_checks
            .iter()
            .filter(|c| c.verdict == verdict)
            .count()
    }
}

/// Folds the design back and compares it to its target structure
/// (feature 15).
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if the design and the target differ in
///   length.
/// - [`RnaDesignError::Upstream`] if the folder fails.
pub fn fold_back_validate(seq: &[u8], target: &Structure) -> Result<FoldBackResult> {
    if seq.len() != target.len() {
        return Err(RnaDesignError::invalid(
            "target",
            "design and target structure differ in length",
        ));
    }
    let rna = RnaSeq::parse(seq)?;
    let folded = mfe(&rna)?;
    let dist = base_pair_distance(&folded.structure, target)?;

    let target_pairs = target.n_pairs();
    let match_pct = if target_pairs == 0 {
        100.0
    } else {
        let achieved_pairs = folded.structure.n_pairs();
        let shared = (achieved_pairs + target_pairs).saturating_sub(dist) / 2;
        100.0 * shared as f64 / target_pairs as f64
    };

    Ok(FoldBackResult {
        achieved_dot_bracket: folded.structure.to_dot_bracket(),
        target_dot_bracket: target.to_dot_bracket(),
        base_pair_distance: dist,
        structure_match_percent: match_pct,
        mfe_energy: folded.energy,
    })
}

/// Characterises the design's Boltzmann ensemble (feature 16).
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if the design and the target differ in
///   length.
/// - [`RnaDesignError::Upstream`] if the partition function fails.
pub fn ensemble_validate(seq: &[u8], target: &Structure) -> Result<EnsembleResult> {
    if seq.len() != target.len() {
        return Err(RnaDesignError::invalid(
            "target",
            "design and target structure differ in length",
        ));
    }
    let rna = RnaSeq::parse(seq)?;
    let pf = partition_function(&rna)?;
    let g = pf.ensemble_free_energy();

    // P(target) = exp(-(E_target - G) / RT). If the target has a
    // non-canonical pair on this sequence its energy is undefined —
    // probability 0.
    let target_probability = match structure_energy(&rna, target) {
        Ok(e_target) => {
            // RT at 37 C in kcal/mol — match the partition function's
            // reference temperature.
            let rt = 1.987_204e-3 * pf.temperature_k();
            ((-(e_target - g) / rt).exp()).clamp(0.0, 1.0)
        }
        Err(_) => 0.0,
    };

    let defect = ensemble_defect(seq, target)?;
    let normalized = if seq.is_empty() {
        0.0
    } else {
        defect / seq.len() as f64
    };

    Ok(EnsembleResult {
        ensemble_free_energy: g,
        target_probability,
        ensemble_defect: defect,
        normalized_defect: normalized,
    })
}

/// Runs the robustness sanity checks on a design (feature 17).
///
/// `intended` is the structural target when the design has one. The
/// checks are: **mutational robustness** (sample single point mutations,
/// measure how often the fold stays close), a **melting / temperature
/// sanity check** (a melting curve), and a **co-transcriptional sanity
/// check** (a 5′-window fold progression — does a strong structure
/// appear so early it could trap a mis-fold during transcription).
///
/// # Errors
/// [`RnaDesignError::Upstream`] if a folding call fails.
pub fn robustness_checks(seq: &[u8], intended: Option<&Structure>) -> Result<RobustnessResult> {
    let mut notes = Vec::new();

    // --- mutational robustness ---------------------------------------
    // Sample point mutations spread across the sequence; for each, fold
    // the mutant and measure base-pair distance to the wild-type fold.
    let wt_rna = RnaSeq::parse(seq)?;
    let wt_fold = mfe(&wt_rna)?.structure;
    let n = seq.len();
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];
    let mut robust = 0usize;
    let mut sampled = 0usize;
    // Sample at most ~20 positions, evenly spaced, one mutation each.
    let stride = (n / 20).max(1);
    let mut pos = 0;
    while pos < n {
        let original = seq[pos].to_ascii_uppercase();
        // Mutate to the next base in the cycle (deterministic).
        let new = BASES
            .iter()
            .copied()
            .find(|&b| b != original)
            .unwrap_or(b'A');
        let mut mutant = seq.to_vec();
        mutant[pos] = new;
        if let Ok(m_rna) = RnaSeq::parse(&mutant) {
            if let Ok(m_fold) = mfe(&m_rna) {
                sampled += 1;
                let d = base_pair_distance(&m_fold.structure, &wt_fold)?;
                // "Robust" = the mutant fold is within 2 pairs of WT.
                if d <= 2 {
                    robust += 1;
                }
            }
        }
        pos += stride;
    }
    let mutational_robustness = if sampled > 0 {
        robust as f64 / sampled as f64
    } else {
        1.0
    };
    notes.push(format!(
        "Mutational robustness: {:.0}% of {sampled} sampled point mutations keep the \
         MFE fold within 2 base pairs of wild-type.",
        mutational_robustness * 100.0,
    ));

    // --- melting / temperature sanity --------------------------------
    let melting_temperature_c = if n >= 4 {
        match melting_curve(&wt_rna, 10.0, 95.0, 17) {
            Ok(curve) => {
                let tm = curve.tm();
                match tm {
                    Some(t) => notes.push(format!(
                        "Melting sanity: a cooperative transition near {t:.0} C.",
                    )),
                    None => notes.push(
                        "Melting sanity: no cooperative transition in 10-95 C — the design \
                         has little stable structure (expected for an unstructured RNA)."
                            .to_string(),
                    ),
                }
                tm
            }
            Err(_) => None,
        }
    } else {
        None
    };

    // --- co-transcriptional sanity -----------------------------------
    // As the 5' end emerges first during transcription, a very strong
    // structure in the first ~25 nt can kinetically trap a mis-fold. We
    // fold the growing 5' prefix at a few lengths and flag a prefix
    // that folds far more strongly than the equivalent region of the
    // intended structure would.
    let cotranscriptional_ok = cotranscriptional_sanity(seq, intended, &mut notes)?;

    Ok(RobustnessResult {
        mutational_robustness,
        melting_temperature_c,
        cotranscriptional_ok,
        notes,
    })
}

/// The co-transcriptional-folding sanity heuristic. Returns `true` when
/// no red flag is found.
fn cotranscriptional_sanity(
    seq: &[u8],
    intended: Option<&Structure>,
    notes: &mut Vec<String>,
) -> Result<bool> {
    let n = seq.len();
    if n < 25 {
        notes.push(
            "Co-transcriptional sanity: sequence too short to assess premature 5' folding."
                .to_string(),
        );
        return Ok(true);
    }
    // Fold the first 25 nt — the earliest stretch the polymerase emits.
    let prefix = RnaSeq::parse(&seq[..25])?;
    let prefix_mfe = mfe(&prefix)?.energy;
    // A premature 5' structure stronger than -12 kcal/mol is a kinetic
    // trapping risk — *unless* the intended fold genuinely wants a 5'
    // structure there.
    let prefix_wants_structure = match intended {
        Some(target) => {
            let mut paired = 0;
            for i in 0..25.min(target.len()) {
                if let Some(j) = target.partner(i) {
                    if j < 25 {
                        paired += 1;
                    }
                }
            }
            paired >= 8
        }
        None => false,
    };
    if prefix_mfe < -12.0 && !prefix_wants_structure {
        notes.push(format!(
            "Co-transcriptional sanity: the first 25 nt fold strongly ({prefix_mfe:.1} \
             kcal/mol) with no such structure in the target — a kinetic mis-folding risk \
             during transcription.",
        ));
        Ok(false)
    } else {
        notes.push(
            "Co-transcriptional sanity: no strong premature 5' structure — the design is \
             unlikely to kinetically trap a mis-fold as it is transcribed."
                .to_string(),
        );
        Ok(true)
    }
}

/// Runs the full in-silico validation of a design (feature 18).
///
/// Checks every constraint, folds the design back to its target (when
/// it has one), characterises the ensemble, runs the robustness sanity
/// checks, and assembles the [`ValidationReport`].
///
/// # Errors
/// [`RnaDesignError::Upstream`] if a folding / partition-function call
/// fails.
pub fn validate_design(
    design: &RnaDesign,
    constraints: &DesignConstraints,
) -> Result<ValidationReport> {
    let seq = &design.sequence;
    let mut checks = Vec::new();

    // --- constraint: length ------------------------------------------
    {
        let len = seq.len();
        let verdict = if constraints.length_ok(len) {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Fail
        };
        checks.push(ConstraintCheck {
            name: "length".to_string(),
            verdict,
            detail: format!(
                "design length {len} nt; required window [{}, {}]",
                fmt_bound(constraints.length_min),
                fmt_bound(constraints.length_max),
            ),
        });
    }

    // --- constraint: GC content --------------------------------------
    {
        let gc = gc_content(seq);
        let verdict = if constraints.gc_ok(gc) {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Fail
        };
        checks.push(ConstraintCheck {
            name: "GC content".to_string(),
            verdict,
            detail: format!(
                "GC {:.0}%; required {:.0}-{:.0}%",
                gc * 100.0,
                constraints.gc_min * 100.0,
                constraints.gc_max * 100.0,
            ),
        });
    }

    // --- constraint: forbidden motifs / restriction sites ------------
    {
        let scan = forbidden_motif_scan(seq, constraints);
        let verdict = if scan.is_clean() {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Fail
        };
        checks.push(ConstraintCheck {
            name: "forbidden motifs".to_string(),
            verdict,
            detail: if scan.is_clean() {
                "no forbidden motif or restriction site found".to_string()
            } else {
                format!("{} forbidden hit(s) found", scan.count())
            },
        });
    }

    // --- constraint: required subsequences ---------------------------
    {
        let rna = transcribe(seq);
        let mut missing = Vec::new();
        for req in &constraints.required_subsequences {
            let needle = transcribe(req.as_bytes());
            if !needle.is_empty() && !contains(&rna, &needle) {
                missing.push(req.clone());
            }
        }
        let verdict = if missing.is_empty() {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Fail
        };
        checks.push(ConstraintCheck {
            name: "required subsequences".to_string(),
            verdict,
            detail: if missing.is_empty() {
                format!(
                    "all {} required subsequence(s) present",
                    constraints.required_subsequences.len()
                )
            } else {
                format!("missing required subsequence(s): {}", missing.join(", "))
            },
        });
    }

    // --- constraint: repeats / low complexity ------------------------
    {
        let scan = repeat_scan(seq, constraints.max_homopolymer);
        let verdict = if scan.is_clean() {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Warn
        };
        checks.push(ConstraintCheck {
            name: "repeats / low complexity".to_string(),
            verdict,
            detail: if scan.is_clean() {
                "no repeat or low-complexity stretch found".to_string()
            } else {
                format!("{} repeat / low-complexity stretch(es) found", scan.count())
            },
        });
    }

    // --- constraint: synthesizability --------------------------------
    {
        let scan = synthesizability_scan(seq, constraints.max_homopolymer)?;
        let verdict = if scan.is_synthesizable() {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Warn
        };
        checks.push(ConstraintCheck {
            name: "synthesizability".to_string(),
            verdict,
            detail: if scan.is_synthesizable() {
                "no synthesis liability (GC, homopolymer, 5' structure all ok)".to_string()
            } else {
                format!("{} synthesis warning(s)", scan.warnings.len())
            },
        });
    }

    // --- structural validation (designs with a target) ---------------
    let intended = intended_structure(design);
    let (fold_back, ensemble) = if let Some(target) = intended.as_ref() {
        if target.len() == seq.len() {
            let fb = fold_back_validate(seq, target)?;
            let en = ensemble_validate(seq, target)?;
            // A structure-match check folds the fold-back into the
            // verdict.
            let fb_verdict = if fb.is_exact() {
                ValidationVerdict::Pass
            } else if fb.structure_match_percent >= 80.0 {
                ValidationVerdict::Warn
            } else {
                ValidationVerdict::Fail
            };
            checks.push(ConstraintCheck {
                name: "target structure".to_string(),
                verdict: fb_verdict,
                detail: format!(
                    "design folds to {:.0}% of the target ({} base pair(s) differ); \
                     target probability {:.1}%",
                    fb.structure_match_percent,
                    fb.base_pair_distance,
                    en.target_probability * 100.0,
                ),
            });
            (Some(fb), Some(en))
        } else {
            (None, None)
        }
    } else {
        (None, None)
    };

    // --- robustness checks -------------------------------------------
    let robustness = Some(robustness_checks(seq, intended.as_ref())?);
    if let Some(r) = &robustness {
        // A flagged co-transcriptional fold or low mutational robustness
        // is a non-fatal concern — never a hard fail, since robustness
        // is an advisory sanity check, not a constraint.
        let verdict = if r.cotranscriptional_ok && r.mutational_robustness >= 0.5 {
            ValidationVerdict::Pass
        } else {
            ValidationVerdict::Warn
        };
        checks.push(ConstraintCheck {
            name: "robustness".to_string(),
            verdict,
            detail: format!(
                "mutational robustness {:.0}%; co-transcriptional folding {}",
                r.mutational_robustness * 100.0,
                if r.cotranscriptional_ok {
                    "ok"
                } else {
                    "flagged"
                },
            ),
        });
    }

    // --- aggregate verdict -------------------------------------------
    let verdict = checks
        .iter()
        .map(|c| c.verdict)
        .fold(ValidationVerdict::Pass, ValidationVerdict::worst);

    let mut notes = Vec::new();
    notes.push(format!(
        "In-silico validation: {} pass, {} warn, {} fail across {} check(s).",
        count_verdict(&checks, ValidationVerdict::Pass),
        count_verdict(&checks, ValidationVerdict::Warn),
        count_verdict(&checks, ValidationVerdict::Fail),
        checks.len(),
    ));
    notes.push(
        "IMPORTANT: this is an in-silico prediction from a nearest-neighbor energy model \
         — a strong validated candidate, NOT a guarantee of in-vivo behaviour. The physical \
         RNA must be synthesised and validated in the wet lab."
            .to_string(),
    );

    Ok(ValidationReport {
        verdict,
        constraint_checks: checks,
        fold_back,
        ensemble,
        robustness,
        notes,
    })
}

/// The structural target a design should fold to, if any.
fn intended_structure(design: &RnaDesign) -> Option<Structure> {
    match &design.kind {
        DesignKind::Structural { target } => Some(target.clone()),
        DesignKind::Riboswitch { free_target, .. } => Some(free_target.clone()),
        _ => None,
    }
}

/// Transcribes raw bytes to uppercase RNA.
fn transcribe(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'T' => b'U',
            other => other,
        })
        .collect()
}

/// `true` if `needle` occurs in `haystack` (case-insensitive).
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || needle.len() > haystack.len() {
        return false;
    }
    haystack
        .windows(needle.len())
        .any(|w| w.eq_ignore_ascii_case(needle))
}

/// Counts the constraint checks with a given verdict.
fn count_verdict(checks: &[ConstraintCheck], verdict: ValidationVerdict) -> usize {
    checks.iter().filter(|c| c.verdict == verdict).count()
}

/// Formats a `0`-means-unbounded length bound for display.
fn fmt_bound(bound: usize) -> String {
    if bound == 0 {
        "-".to_string()
    } else {
        bound.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::DesignKind;

    fn structural_design(seq: &[u8], db: &str) -> RnaDesign {
        RnaDesign {
            sequence: seq.to_vec(),
            kind: DesignKind::Structural {
                target: Structure::from_dot_bracket(db).unwrap(),
            },
            cds_span: None,
            construct: None,
            notes: Vec::new(),
        }
    }

    #[test]
    fn verdict_worst_aggregates() {
        use ValidationVerdict::*;
        assert_eq!(Pass.worst(Warn), Warn);
        assert_eq!(Warn.worst(Fail), Fail);
        assert_eq!(Pass.worst(Pass), Pass);
        assert_eq!(Fail.worst(Pass), Fail);
    }

    #[test]
    fn fold_back_on_a_hairpin() {
        // A strong GC hairpin should fold close to its target.
        let seq = b"GGGGGGAAAACCCCCC";
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let fb = fold_back_validate(seq, &target).unwrap();
        assert!(fb.mfe_energy < 0.0);
        assert!((0.0..=100.0).contains(&fb.structure_match_percent));
        assert_eq!(fb.target_dot_bracket, "((((((....))))))");
    }

    #[test]
    fn fold_back_rejects_length_mismatch() {
        let target = Structure::from_dot_bracket("((....))").unwrap();
        assert!(fold_back_validate(b"GGGG", &target).is_err());
    }

    #[test]
    fn ensemble_validate_gives_a_probability() {
        let seq = b"GGGGGGAAAACCCCCC";
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let en = ensemble_validate(seq, &target).unwrap();
        assert!((0.0..=1.0).contains(&en.target_probability));
        assert!(en.ensemble_defect >= 0.0);
        assert!((0.0..=1.0).contains(&en.normalized_defect));
    }

    #[test]
    fn robustness_checks_run() {
        let seq = b"GGGGGGAAAACCCCCC";
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let r = robustness_checks(seq, Some(&target)).unwrap();
        assert!((0.0..=1.0).contains(&r.mutational_robustness));
        assert!(!r.notes.is_empty());
    }

    #[test]
    fn validate_a_structural_design() {
        let d = structural_design(b"GGGGGGAAAACCCCCC", "((((((....))))))");
        let report = validate_design(&d, &DesignConstraints::default()).unwrap();
        assert!(!report.constraint_checks.is_empty());
        assert!(report.fold_back.is_some());
        assert!(report.ensemble.is_some());
        assert!(report.robustness.is_some());
        // The disclaimer is present.
        assert!(report
            .notes
            .iter()
            .any(|n| n.contains("in-silico prediction")));
    }

    #[test]
    fn validate_flags_a_forbidden_motif() {
        let d = structural_design(b"GGGAAUAAACCC", "(((......)))");
        let constraints = DesignConstraints::default().forbid_motif("AAUAAA");
        let report = validate_design(&d, &constraints).unwrap();
        // The forbidden-motifs check should fail.
        let mc = report
            .constraint_checks
            .iter()
            .find(|c| c.name == "forbidden motifs")
            .unwrap();
        assert_eq!(mc.verdict, ValidationVerdict::Fail);
        assert_eq!(report.verdict, ValidationVerdict::Fail);
        assert!(!report.passed());
    }

    #[test]
    fn validate_flags_length_violation() {
        let d = structural_design(b"GGGGAAAACCCC", "((((....))))");
        let constraints = DesignConstraints::default().with_length_range(50, 100);
        let report = validate_design(&d, &constraints).unwrap();
        let lc = report
            .constraint_checks
            .iter()
            .find(|c| c.name == "length")
            .unwrap();
        assert_eq!(lc.verdict, ValidationVerdict::Fail);
    }

    #[test]
    fn count_with_tallies_verdicts() {
        let d = structural_design(b"GGGGGGAAAACCCCCC", "((((((....))))))");
        let report = validate_design(&d, &DesignConstraints::default()).unwrap();
        let total = report.count_with(ValidationVerdict::Pass)
            + report.count_with(ValidationVerdict::Warn)
            + report.count_with(ValidationVerdict::Fail);
        assert_eq!(total, report.constraint_checks.len());
    }
}
