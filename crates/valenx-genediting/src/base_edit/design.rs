//! Features 7–10 — base-editing target finding, guide design, window /
//! bystander analysis and a product-purity heuristic.
//!
//! Designing a base edit means answering four questions:
//!
//! 7. **Which edits are reachable?** Scan a target window for every
//!    `C` (CBE) or `A` (ABE) that some guide can place inside its
//!    activity window — [`find_editable_sites`].
//! 8. **Design the guide.** For a *specific* desired SNV, find the
//!    guide that puts that base in the editing window —
//!    [`design_base_edit`].
//! 9. **Window / bystander analysis.** A base editor edits *every*
//!    target base in the window, not just the intended one. Count the
//!    bystanders and report whether a clean edit is possible —
//!    [`analyze_window`].
//! 10. **Product-purity heuristic.** A transparent score for how
//!     "clean" an edit is — penalising bystanders, an off-centre
//!     intended base, and (for CBE) a `GC` motif that drives indels —
//!     [`product_purity`].
//!
//! ## v1 scope — honest caveat
//!
//! The product-purity score is a **transparent feature-weighted
//! heuristic**, not a trained model. The real-world predictors
//! (BE-Hive, DeepCBE / DeepABE) are trained neural networks; the
//! project's "no trained-weights" rule excludes reimplementing them.
//! The heuristic captures the well-established qualitative rules —
//! bystanders hurt purity, an edit dead-centre in the window is
//! cleanest, a CBE `GpC` context raises indel risk — and lands in
//! `[0, 1]` with the right ranking. It is documented as a heuristic
//! everywhere it surfaces.

use crate::base_edit::editor::{base_editor, BaseEditor, BaseEditorClass, BaseEditorId};
use crate::error::{GeneditingError, Result};
use crate::sequtil::is_acgt;
use serde::{Deserialize, Serialize};
use valenx_genomics::crispr::guide::{scan_guides, Guide, GuideStrand, PamSide, PamSpec};

/// 1-based protospacer position (counted from the PAM-distal 5′ end) of
/// a base at forward-strand index `ref_idx`, for a guide found at
/// `guide.start` with length `plen`.
///
/// For a forward-strand 3′-PAM guide, protospacer position 1 is the
/// PAM-distal base at `guide.start`; for a reverse-strand guide the
/// PAM-distal base is at `guide.start + plen - 1` on the forward axis.
fn protospacer_position(guide: &Guide, plen: usize, ref_idx: usize) -> Option<usize> {
    match guide.strand {
        GuideStrand::Forward => {
            if ref_idx < guide.start || ref_idx >= guide.start + plen {
                return None;
            }
            Some(ref_idx - guide.start + 1)
        }
        GuideStrand::Reverse => {
            let last = guide.start + plen; // exclusive
            if ref_idx < guide.start || ref_idx >= last {
                return None;
            }
            // PAM-distal end is the high-coordinate end on the forward
            // axis; position 1 is there.
            Some(last - 1 - ref_idx + 1)
        }
    }
}

/// One editable site discovered by [`find_editable_sites`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EditableSite {
    /// 0-based forward-strand index of the target base.
    pub ref_pos: usize,
    /// The reference base on the forward strand at `ref_pos`.
    pub ref_base: u8,
    /// `true` when the guide that reaches this base is on the reverse
    /// strand (so the editor acts on the `C`/`A` of the *reverse*
    /// strand, i.e. a forward-strand `G`/`T`).
    pub guide_reverse: bool,
    /// The protospacer (5′→3′) of the guide that reaches the site.
    pub protospacer: String,
    /// The PAM as found.
    pub pam: String,
    /// 0-based forward-strand start of that protospacer.
    pub guide_start: usize,
    /// 1-based protospacer position of the target base (from the
    /// PAM-distal 5′ end).
    pub window_position: usize,
}

/// Scans a target window for every base a base editor can reach.
///
/// For a CBE the targetable base on the protospacer (edited) strand is
/// `C`; on a forward-strand guide that is a forward-strand `C`, on a
/// reverse-strand guide a forward-strand `G`. (ABE: `A` / forward `T`.)
/// A site is returned only when some guide places it inside the
/// editor's activity window.
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] for a non-ACGT target.
pub fn find_editable_sites(target: &[u8], editor_id: BaseEditorId) -> Result<Vec<EditableSite>> {
    if !is_acgt(target) {
        return Err(GeneditingError::invalid_target(
            "region",
            "target window must be non-empty ACGT",
        ));
    }
    let editor = base_editor(editor_id);
    let spec = PamSpec {
        motif: editor.pam.clone(),
        protospacer_len: editor.guide_len,
        side: PamSide::ThreePrime,
    };
    let guides = scan_guides(target, &spec)
        .map_err(|e| GeneditingError::invalid("target", e.to_string()))?;

    let edited_strand_base = editor.target_base();
    // The forward-strand base that corresponds to an edited-strand
    // target base on a reverse-strand guide.
    let fwd_base_for_reverse = match editor.class {
        BaseEditorClass::Cbe => b'G', // reverse-strand C == forward G
        BaseEditorClass::Abe => b'T', // reverse-strand A == forward T
    };

    let mut sites = Vec::new();
    for g in &guides {
        let plen = editor.guide_len;
        for (offset, &raw) in target[g.start..g.start + plen].iter().enumerate() {
            let ref_pos = g.start + offset;
            let fwd = raw.to_ascii_uppercase();
            let is_target = match g.strand {
                GuideStrand::Forward => fwd == edited_strand_base,
                GuideStrand::Reverse => fwd == fwd_base_for_reverse,
            };
            if !is_target {
                continue;
            }
            let pos = match protospacer_position(g, plen, ref_pos) {
                Some(p) => p,
                None => continue,
            };
            if !editor.position_in_window(pos) {
                continue;
            }
            sites.push(EditableSite {
                ref_pos,
                ref_base: fwd,
                guide_reverse: g.strand == GuideStrand::Reverse,
                protospacer: g.protospacer.clone(),
                pam: g.pam.clone(),
                guide_start: g.start,
                window_position: pos,
            });
        }
    }
    Ok(sites)
}

/// Window / bystander analysis for one chosen guide and intended base
/// (feature 9).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WindowAnalysis {
    /// 1-based protospacer position of the intended target base.
    pub intended_position: usize,
    /// 1-based positions of *bystander* editable bases inside the
    /// window (other `C`s for a CBE / `A`s for an ABE).
    pub bystander_positions: Vec<usize>,
    /// `true` when the intended base is the only editable base in the
    /// window — a clean single-base edit is possible.
    pub clean: bool,
}

impl WindowAnalysis {
    /// Number of bystander bases in the window.
    pub fn bystander_count(&self) -> usize {
        self.bystander_positions.len()
    }
}

/// Analyses the activity window of a guide for bystander edits.
///
/// `protospacer` is the guide 5′→3′ on its strand; `intended_position`
/// is the 1-based PAM-distal-counted position of the base the caller
/// wants edited. Every *other* editable base inside the editor's
/// window is reported as a bystander.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGT protospacer.
/// - [`GeneditingError::Invalid`] if `intended_position` is outside the
///   protospacer or outside the editor's activity window.
pub fn analyze_window(
    protospacer: &[u8],
    editor_id: BaseEditorId,
    intended_position: usize,
) -> Result<WindowAnalysis> {
    if !is_acgt(protospacer) {
        return Err(GeneditingError::invalid_target(
            "locus",
            "protospacer must be non-empty ACGT",
        ));
    }
    let editor = base_editor(editor_id);
    if intended_position == 0 || intended_position > protospacer.len() {
        return Err(GeneditingError::invalid(
            "intended_position",
            "intended position is outside the protospacer",
        ));
    }
    if !editor.position_in_window(intended_position) {
        return Err(GeneditingError::invalid(
            "intended_position",
            "intended position is outside the editor's activity window",
        ));
    }
    let target = editor.target_base();
    let mut bystanders = Vec::new();
    for (i, &b) in protospacer.iter().enumerate() {
        let pos = i + 1; // 1-based, PAM-distal counting (forward guide)
        if pos == intended_position {
            continue;
        }
        if editor.position_in_window(pos) && b.to_ascii_uppercase() == target {
            bystanders.push(pos);
        }
    }
    Ok(WindowAnalysis {
        intended_position,
        bystander_positions: bystanders.clone(),
        clean: bystanders.is_empty(),
    })
}

/// A transparent product-purity heuristic in `[0, 1]` (feature 10).
///
/// Higher = a cleaner, more predictable edit. The score blends:
///
/// - **bystander penalty** — each bystander editable base in the window
///   sharply reduces purity (mixed alleles);
/// - **centring bonus** — an intended base near the centre of the
///   window edits more reliably than one at a window edge;
/// - **CBE context penalty** — for a CBE, a `G` immediately 5′ of the
///   target `C` (a `GpC` motif) raises uracil-glycosylase-driven
///   indel formation; this is penalised.
///
/// It is a heuristic, not a calibrated probability and not a trained
/// model — see the module note.
pub fn product_purity(
    protospacer: &[u8],
    editor_id: BaseEditorId,
    intended_position: usize,
) -> Result<f64> {
    let analysis = analyze_window(protospacer, editor_id, intended_position)?;
    let editor = base_editor(editor_id);

    // Bystander penalty: each bystander costs a chunk of purity.
    let bystander_penalty = 0.22 * analysis.bystander_count() as f64;

    // Centring bonus: distance from the window centre, normalised.
    let centre = (editor.window_start + editor.window_end) as f64 / 2.0;
    let half = (editor.window_width() as f64 / 2.0).max(1.0);
    let off_centre = ((intended_position as f64 - centre).abs() / half).min(1.0);
    let centre_bonus = 0.15 * (1.0 - off_centre);

    // CBE GpC context penalty (the base 5' of the target C).
    let mut context_penalty = 0.0;
    if editor.class == BaseEditorClass::Cbe && intended_position >= 2 {
        let preceding = protospacer[intended_position - 2].to_ascii_uppercase();
        if preceding == b'G' {
            context_penalty = 0.12;
        }
    }

    let score = 0.85 + centre_bonus - bystander_penalty - context_penalty;
    Ok(score.clamp(0.0, 1.0))
}

/// A request to design a base-editing guide for a specific SNV
/// (feature 8).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BaseEditRequest {
    /// The reference window around the variant (forward strand, ACGT).
    pub reference: Vec<u8>,
    /// 0-based forward-strand index of the base to edit.
    pub edit_pos: usize,
    /// The reference base expected at `edit_pos` (a sanity check).
    pub from_base: u8,
    /// The base the caller wants installed at `edit_pos`.
    pub to_base: u8,
    /// Which base editor to use.
    pub editor: BaseEditorId,
}

/// A designed base-editing guide for a specific SNV.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BaseEditDesign {
    /// The chosen guide (protospacer, 5′→3′ on its strand).
    pub protospacer: String,
    /// The PAM as found.
    pub pam: String,
    /// 0-based forward-strand start of the protospacer.
    pub guide_start: usize,
    /// `true` when the guide is on the reverse strand.
    pub guide_reverse: bool,
    /// 1-based protospacer position of the edited base.
    pub window_position: usize,
    /// Bystander analysis for the chosen guide.
    pub window: WindowAnalysis,
    /// The transparent product-purity heuristic in `[0, 1]`.
    pub product_purity: f64,
}

/// Designs a base-editing guide that installs a specific SNV.
///
/// Confirms the requested transition is one a base editor can make
/// (`C→T` / `G→A` for a CBE, `A→G` / `T→C` for an ABE), finds the
/// guide that places the target base inside the editor's window, runs
/// bystander analysis and scores product purity. When several guides
/// reach the base, the one with the highest product purity is chosen.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGT reference, an
///   out-of-range `edit_pos`, or a `from_base` that disagrees with the
///   reference.
/// - [`GeneditingError::Invalid`] if the requested transition is not
///   one the chosen editor class can install.
/// - [`GeneditingError::NoValidDesign`] if no guide places the target
///   base inside the activity window.
pub fn design_base_edit(req: &BaseEditRequest) -> Result<BaseEditDesign> {
    if !is_acgt(&req.reference) {
        return Err(GeneditingError::invalid_target(
            "region",
            "reference window must be non-empty ACGT",
        ));
    }
    if req.edit_pos >= req.reference.len() {
        return Err(GeneditingError::invalid_target(
            "variant",
            "edit position is outside the reference window",
        ));
    }
    let ref_base = req.reference[req.edit_pos].to_ascii_uppercase();
    if ref_base != req.from_base.to_ascii_uppercase() {
        return Err(GeneditingError::invalid_target(
            "variant",
            format!(
                "reference base at the edit position is `{}`, not the requested \
                 from-base `{}`",
                ref_base as char,
                req.from_base.to_ascii_uppercase() as char
            ),
        ));
    }
    let editor = base_editor(req.editor);
    // Which strand must the guide be on to make this edit?
    //   CBE: forward-strand C->T  ⇒ forward guide;  G->A ⇒ reverse guide.
    //   ABE: forward-strand A->G  ⇒ forward guide;  T->C ⇒ reverse guide.
    let from = req.from_base.to_ascii_uppercase();
    let to = req.to_base.to_ascii_uppercase();
    let need_reverse = match (editor.class, from, to) {
        (BaseEditorClass::Cbe, b'C', b'T') => false,
        (BaseEditorClass::Cbe, b'G', b'A') => true,
        (BaseEditorClass::Abe, b'A', b'G') => false,
        (BaseEditorClass::Abe, b'T', b'C') => true,
        _ => {
            return Err(GeneditingError::invalid(
                "transition",
                format!(
                    "{} cannot install {}->{} ({})",
                    editor.name,
                    from as char,
                    to as char,
                    editor.class.transition()
                ),
            ));
        }
    };

    // Collect every guide that puts `edit_pos` inside the window.
    let sites = find_editable_sites(&req.reference, req.editor)?;
    let mut best: Option<BaseEditDesign> = None;
    for s in &sites {
        if s.ref_pos != req.edit_pos || s.guide_reverse != need_reverse {
            continue;
        }
        let window = analyze_window(s.protospacer.as_bytes(), req.editor, s.window_position)?;
        let purity = product_purity(s.protospacer.as_bytes(), req.editor, s.window_position)?;
        let candidate = BaseEditDesign {
            protospacer: s.protospacer.clone(),
            pam: s.pam.clone(),
            guide_start: s.guide_start,
            guide_reverse: s.guide_reverse,
            window_position: s.window_position,
            window,
            product_purity: purity,
        };
        best = match best {
            None => Some(candidate),
            Some(b) if candidate.product_purity > b.product_purity => Some(candidate),
            Some(b) => Some(b),
        };
    }
    best.ok_or_else(|| {
        GeneditingError::no_valid_design(
            "base_edit",
            "no guide places the target base inside the editor's activity window",
        )
    })
}

/// Convenience: the editor whose window best centres a base at
/// PAM-distal-counted `pos`, among the editors of the right `class`.
/// Returns `None` if no editor of that class can reach `pos`.
pub fn best_editor_for_position(class: BaseEditorClass, pos: usize) -> Option<BaseEditor> {
    BaseEditorId::all()
        .into_iter()
        .map(base_editor)
        .filter(|e| e.class == class && e.position_in_window(pos))
        .min_by(|a, b| {
            let ca = ((a.window_start + a.window_end) as f64 / 2.0 - pos as f64).abs();
            let cb = ((b.window_start + b.window_end) as f64 / 2.0 - pos as f64).abs();
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_editable_c_for_cbe() {
        // A protospacer with a C at position 5 (window for BE4max is
        // 4-8), 20-mer + NGG PAM.
        let target = b"ACGTCACGTACGTACGTACGTGG"; // C at index 4 -> pos 5
        let sites = find_editable_sites(target, BaseEditorId::Be4Max).unwrap();
        assert!(sites
            .iter()
            .any(|s| s.window_position >= 4 && s.window_position <= 8));
    }

    #[test]
    fn rejects_non_acgt_target() {
        assert!(find_editable_sites(b"ACGTNNNN", BaseEditorId::Be3).is_err());
    }

    #[test]
    fn bystander_analysis_counts_extra_cs() {
        // Protospacer with C at window positions 4 and 6 (BE4max 4-8).
        // index 3 -> pos 4, index 5 -> pos 6.
        let proto = b"ACGCACGTACGTACGTACGT";
        let a = analyze_window(proto, BaseEditorId::Be4Max, 4).unwrap();
        // Position 6 is a C inside the window ⇒ one bystander.
        assert!(a.bystander_count() >= 1);
        assert!(!a.clean);
    }

    #[test]
    fn clean_window_has_no_bystanders() {
        // Only one C inside the 4-8 window: index 3 -> pos 4.
        let proto = b"ACGCATGTATGTATGTATGT";
        let a = analyze_window(proto, BaseEditorId::Be4Max, 4).unwrap();
        assert!(a.clean);
        assert_eq!(a.bystander_count(), 0);
    }

    #[test]
    fn intended_position_outside_window_rejected() {
        let proto = b"CCGTACGTACGTACGTACGT";
        // Position 1 is outside the 4-8 window.
        let err = analyze_window(proto, BaseEditorId::Be4Max, 1).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn purity_drops_with_bystanders() {
        let clean = b"ACGCATGTATGTATGTATGT"; // one C at pos 4
        let dirty = b"ACGCACGCACGTACGTACGT"; // multiple C in window
        let pc = product_purity(clean, BaseEditorId::Be4Max, 4).unwrap();
        let pd = product_purity(dirty, BaseEditorId::Be4Max, 4).unwrap();
        assert!(pc > pd);
    }

    #[test]
    fn purity_in_unit_range() {
        let proto = b"ACGCACGCACGCACGCACGC";
        let p = product_purity(proto, BaseEditorId::Be4Max, 4).unwrap();
        assert!((0.0..=1.0).contains(&p));
    }

    #[test]
    fn cbe_gpc_context_lowers_purity() {
        // Target C at pos 5 (index 4). Preceding base at index 3.
        // "ACGGC..." -> index 3 = G -> GpC penalty.
        let gpc = b"ACGGCATGTATGTATGTATG"; // G before the target C
        let apc = b"ACGACATGTATGTATGTATG"; // A before the target C
        let pg = product_purity(gpc, BaseEditorId::Be4Max, 5).unwrap();
        let pa = product_purity(apc, BaseEditorId::Be4Max, 5).unwrap();
        assert!(pg < pa);
    }

    #[test]
    fn designs_a_c_to_t_edit() {
        // Reference with a C at index 4; BE4max window 4-8 → the
        // forward guide starting at 0 puts it at protospacer pos 5.
        let reference = b"ACGTCATGTATGTATGTATGTGG".to_vec();
        let req = BaseEditRequest {
            reference,
            edit_pos: 4,
            from_base: b'C',
            to_base: b'T',
            editor: BaseEditorId::Be4Max,
        };
        let d = design_base_edit(&req).unwrap();
        assert!(!d.guide_reverse);
        assert!(d.window_position >= 4 && d.window_position <= 8);
    }

    #[test]
    fn rejects_impossible_transition() {
        let reference = b"ACGTCATGTATGTATGTATGTGG".to_vec();
        let req = BaseEditRequest {
            reference,
            edit_pos: 4,
            from_base: b'C',
            to_base: b'G', // a CBE cannot do C->G
            editor: BaseEditorId::Be4Max,
        };
        let err = design_base_edit(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn rejects_wrong_from_base() {
        let reference = b"ACGTCATGTATGTATGTATGTGG".to_vec();
        let req = BaseEditRequest {
            reference,
            edit_pos: 4,
            from_base: b'A', // index 4 is actually C
            to_base: b'G',
            editor: BaseEditorId::Abe8e,
        };
        let err = design_base_edit(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn no_design_when_base_unreachable() {
        // A C with no PAM placing it in any window.
        let reference = b"CAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        let req = BaseEditRequest {
            reference,
            edit_pos: 0,
            from_base: b'C',
            to_base: b'T',
            editor: BaseEditorId::Be4Max,
        };
        let err = design_base_edit(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.no_valid_design");
    }

    #[test]
    fn best_editor_for_position_picks_a_class_member() {
        let e = best_editor_for_position(BaseEditorClass::Abe, 6).unwrap();
        assert_eq!(e.class, BaseEditorClass::Abe);
        assert!(e.position_in_window(6));
        // Position 1 is in no editor's window.
        assert!(best_editor_for_position(BaseEditorClass::Cbe, 1).is_none());
    }
}
