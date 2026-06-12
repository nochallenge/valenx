//! Feature 15 — prime-editing strategy modelling and a transparent
//! efficiency heuristic.
//!
//! Picking a prime-editing strategy means choosing an editor
//! configuration (PE2 / PE3 / PE3b / PEmax) and pegRNA parameters
//! (PBS length, RT-template length) and estimating how well they will
//! work. This module provides:
//!
//! - [`length_scan_score`] — the heuristic the
//!   [`crate::prime_edit::pegrna`] length scan ranks PBS / RT-template
//!   combinations with;
//! - [`prime_efficiency`] — a whole-strategy efficiency heuristic that
//!   blends the editor architecture, the pegRNA quality and the edit
//!   size;
//! - [`model_strategy`] — a one-call strategy summary
//!   ([`StrategyModel`]).
//!
//! ## v1 scope — honest caveat
//!
//! **Both scores are transparent feature-weighted heuristics, not
//! trained models.** The real-world prime-editing efficiency
//! predictors (the PRIDICT / DeepPrime neural networks) are trained on
//! large screens; the project's "no trained-weights" rule excludes
//! reimplementing them. The heuristics here encode the
//! well-established qualitative rules:
//!
//! - a PBS melting temperature near ~30 °C anneals well — too short
//!   under-primes, too long mis-folds;
//! - an RT-template homology of ~10–16 nt is the usual sweet spot;
//! - long edits, long RT templates and homopolymer runs reduce
//!   efficiency;
//! - PE3 > PE2 and PEmax > PE3 by architecture.
//!
//! Every value lands in `[0, 1]` with the right qualitative ranking,
//! and is documented as a heuristic wherever it surfaces.

use crate::prime_edit::editor::{prime_editor, PrimeEditor, PrimeEditorId};
use crate::prime_edit::pegrna::{PegRna, PrimeEdit};
use crate::sequtil::max_homopolymer;
use serde::{Deserialize, Serialize};

/// A rough nearest-neighbour-free melting temperature (°C) for a short
/// duplex, the Wallace ("2+4") rule: `2·(A+T) + 4·(G+C)`. Adequate for
/// the PBS-length heuristic; not a substitute for a full
/// nearest-neighbour Tm.
fn wallace_tm(seq: &[u8]) -> f64 {
    let mut at = 0usize;
    let mut gc = 0usize;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'A' | b'T' | b'U' => at += 1,
            b'G' | b'C' => gc += 1,
            _ => {}
        }
    }
    2.0 * at as f64 + 4.0 * gc as f64
}

/// A triangular "preference" curve: `1.0` at `optimum`, falling
/// linearly to `0.0` at `optimum ± width`, clamped to `[0, 1]`.
fn triangular(value: f64, optimum: f64, width: f64) -> f64 {
    if width <= 0.0 {
        return 0.0;
    }
    (1.0 - (value - optimum).abs() / width).clamp(0.0, 1.0)
}

/// Transparent PBS / RT-template length-scan score in `[0, 1]`
/// (used by [`crate::prime_edit::pegrna::scan_pbs_rt`]).
///
/// Blends four feature terms:
///
/// - **PBS Tm** — a Wallace-rule melting temperature with a triangular
///   optimum at ~30 °C (a well-annealing primer);
/// - **RT-template length** — a triangular optimum at ~13 nt of
///   downstream homology;
/// - **PBS homopolymer penalty** — a long single-base run in the PBS
///   mis-primes;
/// - **RT homopolymer penalty** — likewise for the RT template.
///
/// `pbs` and `rt_template` are the designed RNA sequences;
/// `rt_homology` is the downstream-homology length the caller scanned.
pub fn length_scan_score(
    pbs: &[u8],
    rt_template: &[u8],
    _pbs_len: usize,
    rt_homology: usize,
) -> f64 {
    let tm = wallace_tm(pbs);
    let tm_term = triangular(tm, 30.0, 22.0);
    let rt_term = triangular(rt_homology as f64, 13.0, 12.0);

    let pbs_hp = max_homopolymer(pbs);
    let pbs_pen = if pbs_hp >= 5 { 0.20 } else { 0.0 };
    let rt_hp = max_homopolymer(rt_template);
    let rt_pen = if rt_hp >= 6 { 0.15 } else { 0.0 };

    let score = 0.55 * tm_term + 0.45 * rt_term - pbs_pen - rt_pen;
    score.clamp(0.0, 1.0)
}

/// Transparent whole-strategy prime-editing efficiency heuristic in
/// `[0, 1]`.
///
/// Blends:
///
/// - the **pegRNA quality** — the [`length_scan_score`] of the chosen
///   PBS / RT template;
/// - the **editor architecture** — PE2 / PE3 / PE3b / PEmax, via the
///   editor's `architecture_factor`, normalised so PEmax maps near
///   `1.0`;
/// - an **edit-size penalty** — substitutions edit most easily, then
///   small insertions / deletions; large edits are penalised.
///
/// It is a heuristic ranking score, **not** a predicted editing
/// percentage and **not** a trained model — see the module note.
pub fn prime_efficiency(peg: &PegRna, editor: PrimeEditorId) -> f64 {
    let ed = prime_editor(editor);
    let peg_quality = length_scan_score(
        &peg.pbs,
        &peg.rt_template,
        peg.pbs.len(),
        // The RT-template homology length is the RT template minus the
        // edit footprint; recover an approximate downstream-homology
        // length from the RT template length and the edit.
        rt_homology_estimate(peg),
    );
    // Architecture factor ranges 1.0 (PE2) .. 1.6 (PEmax); normalise to
    // a [0.62, 1.0] multiplier.
    let arch = (ed.architecture_factor / 1.6).clamp(0.0, 1.0);
    let size_penalty = edit_size_penalty(&peg.edit);

    let score = (0.60 * peg_quality + 0.40 * arch) * (1.0 - size_penalty);
    score.clamp(0.0, 1.0)
}

/// Estimates the downstream-homology length of a pegRNA's RT template
/// (RT-template length minus the edit footprint).
fn rt_homology_estimate(peg: &PegRna) -> usize {
    let footprint = match &peg.edit {
        PrimeEdit::Substitution { to, .. } => to.len(),
        PrimeEdit::Insertion { seq } => seq.len(),
        PrimeEdit::Deletion { .. } => 0,
    };
    peg.rt_template.len().saturating_sub(footprint)
}

/// Edit-size penalty in `[0, 1)` — `0` for a 1-bp substitution,
/// growing for larger edits.
fn edit_size_penalty(edit: &PrimeEdit) -> f64 {
    let size = match edit {
        PrimeEdit::Substitution { len, to } => (*len).max(to.len()),
        PrimeEdit::Insertion { seq } => seq.len(),
        PrimeEdit::Deletion { len } => *len,
    };
    // ~0 at 1 bp, ~0.5 at ~40 bp, asymptotic.
    let s = size as f64;
    (s / (s + 40.0)).min(0.85)
}

/// A modelled prime-editing strategy.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StrategyModel {
    /// The prime-editor configuration.
    pub editor: PrimeEditorId,
    /// The editor's display name.
    pub editor_name: String,
    /// `true` when the strategy uses a second nicking guide.
    pub uses_nicking_guide: bool,
    /// The whole-strategy efficiency heuristic in `[0, 1]`.
    pub efficiency: f64,
    /// A one-line rationale.
    pub rationale: String,
}

/// Models a prime-editing strategy for a designed pegRNA (feature 15).
///
/// Wraps [`prime_efficiency`] with a human-readable rationale that
/// names the editor configuration, the indel trade-off and the edit
/// kind.
pub fn model_strategy(peg: &PegRna, editor: PrimeEditorId) -> StrategyModel {
    let ed: PrimeEditor = prime_editor(editor);
    let eff = prime_efficiency(peg, editor);
    let rationale = format!(
        "{} ({}); installs a {}. Heuristic efficiency score {:.2} — \
         a transparent feature score (pegRNA quality + architecture - edit-size \
         penalty), not a trained-model prediction.",
        ed.name,
        if ed.nick_after_edit {
            "PE3b second nick — minimises indels"
        } else if ed.uses_nicking_guide {
            "PE3-style second nick — higher editing, more indels than PE3b"
        } else {
            "no second nick"
        },
        peg.edit.label(),
        eff,
    );
    StrategyModel {
        editor,
        editor_name: ed.name.clone(),
        uses_nicking_guide: ed.uses_nicking_guide,
        efficiency: eff,
        rationale,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prime_edit::pegrna::{design_pegrna, PegRnaRequest};

    fn reference() -> Vec<u8> {
        b"ACGTACGTACGTACGTACGTAGGCATGCATGCATGCATGCATGCATGC".to_vec()
    }

    fn peg(editor: PrimeEditorId) -> PegRna {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), editor);
        design_pegrna(&req).unwrap()
    }

    #[test]
    fn wallace_tm_counts_gc_more() {
        assert!(wallace_tm(b"GGGG") > wallace_tm(b"AAAA"));
        assert!((wallace_tm(b"AT") - 4.0).abs() < 1e-9);
        assert!((wallace_tm(b"GC") - 8.0).abs() < 1e-9);
    }

    #[test]
    fn triangular_peaks_at_optimum() {
        assert!((triangular(13.0, 13.0, 12.0) - 1.0).abs() < 1e-9);
        assert!(triangular(13.0, 13.0, 12.0) > triangular(20.0, 13.0, 12.0));
        assert_eq!(triangular(40.0, 13.0, 12.0), 0.0);
    }

    #[test]
    fn length_scan_score_in_unit_range() {
        let s = length_scan_score(b"ACGUACGUACGU", b"ACGUACGUACGUACGU", 12, 13);
        assert!((0.0..=1.0).contains(&s));
    }

    #[test]
    fn length_scan_penalises_homopolymer_pbs() {
        let clean = length_scan_score(b"ACGUACGUACGU", b"ACGUACGUACGU", 12, 13);
        let runny = length_scan_score(b"AAAAAACGUACG", b"ACGUACGUACGU", 12, 13);
        assert!(clean > runny);
    }

    #[test]
    fn efficiency_in_unit_range() {
        let e = prime_efficiency(&peg(PrimeEditorId::Pe2), PrimeEditorId::Pe2);
        assert!((0.0..=1.0).contains(&e));
    }

    #[test]
    fn pemax_outscores_pe2_for_the_same_pegrna() {
        let p = peg(PrimeEditorId::Pe2);
        let e2 = prime_efficiency(&p, PrimeEditorId::Pe2);
        let emax = prime_efficiency(&p, PrimeEditorId::PeMax);
        assert!(emax > e2);
    }

    #[test]
    fn large_edit_penalty_grows() {
        let small = edit_size_penalty(&PrimeEdit::snv(b'A'));
        let big = edit_size_penalty(&PrimeEdit::Insertion {
            seq: vec![b'A'; 30],
        });
        assert!(big > small);
    }

    #[test]
    fn model_strategy_describes_the_configuration() {
        let m = model_strategy(&peg(PrimeEditorId::Pe3b), PrimeEditorId::Pe3b);
        assert_eq!(m.editor, PrimeEditorId::Pe3b);
        assert!(m.uses_nicking_guide);
        assert!(m.rationale.contains("PE3b"));
        assert!((0.0..=1.0).contains(&m.efficiency));
    }
}
