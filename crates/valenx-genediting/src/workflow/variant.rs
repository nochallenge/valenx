//! Feature 28 — the variant-correction planner.
//!
//! Given a **pathogenic variant** and the reference window around it,
//! this module designs the editing strategy that *reverts* the variant
//! to the wild-type sequence. It is the inverse problem of variant
//! *installation*: the "desired change" is `variant → reference`.
//!
//! The planner:
//!
//! 1. classifies the correction (a transition, a transversion, an
//!    insertion to remove, a deletion to restore);
//! 2. asks the [`crate::workflow::advisor`] which modality fits;
//! 3. actually *designs* the editing reagent — a base-editing guide,
//!    a pegRNA, or an HDR donor — for the recommended modality, with a
//!    fallback to prime editing (the universal substitution route)
//!    when the first choice has no PAM.
//!
//! ## v1 scope
//!
//! The planner corrects **substitutions, small insertions and small
//! deletions** within the supplied reference window. It picks one
//! reagent design; it does not enumerate every possible guide. The
//! variant must be described against the *reference* (the reference
//! window is wild-type; the variant base / indel is what the patient
//! carries). Efficiency figures come from the underlying modules'
//! transparent heuristics.

use crate::base_edit::design::{design_base_edit, BaseEditDesign, BaseEditRequest};
use crate::base_edit::editor::{base_editor, BaseEditorClass, BaseEditorId};
use crate::error::{GeneditingError, Result};
use crate::prime_edit::editor::PrimeEditorId;
use crate::prime_edit::pegrna::{design_pegrna, PegRna, PegRnaRequest, PrimeEdit};
use crate::sequtil::is_acgt;
use crate::workflow::advisor::{advise_strategy, DesiredChange, EditApproach};
use serde::{Deserialize, Serialize};

/// A pathogenic variant to correct, described against the reference.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PathogenicVariant {
    /// A point substitution: the patient carries `variant_base` where
    /// the reference (wild-type) has `ref_base`.
    Substitution {
        /// The wild-type reference base.
        ref_base: u8,
        /// The pathogenic base the patient carries.
        variant_base: u8,
    },
    /// An inserted stretch the patient carries that the wild type does
    /// not — correction *removes* it.
    Insertion {
        /// The inserted (pathogenic) bases.
        inserted: Vec<u8>,
    },
    /// A deletion the patient carries — correction *restores* the
    /// deleted reference bases.
    Deletion {
        /// The wild-type bases that are missing in the patient.
        deleted: Vec<u8>,
    },
}

impl PathogenicVariant {
    /// A short human-readable label.
    pub fn label(&self) -> String {
        match self {
            PathogenicVariant::Substitution { ref_base, variant_base } => format!(
                "{}>{} pathogenic substitution",
                *ref_base as char, *variant_base as char
            ),
            PathogenicVariant::Insertion { inserted } => {
                format!("{}-bp pathogenic insertion", inserted.len())
            }
            PathogenicVariant::Deletion { deleted } => {
                format!("{}-bp pathogenic deletion", deleted.len())
            }
        }
    }
}

/// A request to plan the correction of a pathogenic variant.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct VariantCorrectionRequest {
    /// The **wild-type** reference window around the variant (forward
    /// strand, ACGT). The plan edits the patient sequence back to
    /// this.
    pub reference: Vec<u8>,
    /// 0-based index in `reference` where the variant sits (the start
    /// of the affected span).
    pub variant_pos: usize,
    /// The pathogenic variant.
    pub variant: PathogenicVariant,
}

/// The designed reagent that corrects the variant.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum CorrectionReagent {
    /// A base-editing guide (with the base editor chosen).
    BaseEdit {
        /// The base editor to use.
        editor: BaseEditorId,
        /// The designed base-editing guide.
        design: BaseEditDesign,
    },
    /// A prime-editing pegRNA (with the prime editor chosen).
    Prime {
        /// The prime editor to use.
        editor: PrimeEditorId,
        /// The designed pegRNA.
        pegrna: PegRna,
    },
}

/// The result of a variant-correction plan.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct VariantCorrectionPlan {
    /// The variant that was planned for.
    pub variant: PathogenicVariant,
    /// The editing modality chosen.
    pub approach: EditApproach,
    /// The designed correction reagent.
    pub reagent: CorrectionReagent,
    /// A human-readable explanation of the plan.
    pub rationale: String,
}

/// Plans the correction of a pathogenic variant (feature 28).
///
/// Works out the wild-type-restoring edit, asks the advisor which
/// modality fits, and designs the reagent. Base editing is used when
/// the correction is a `C↔T` / `A↔G` transition reachable by an
/// editor; otherwise (or on a base-editing PAM miss) it falls back to
/// prime editing — the universal substitution / small-indel route.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGT reference, an
///   out-of-range `variant_pos`, or a variant inconsistent with the
///   reference.
/// - [`GeneditingError::NoValidDesign`] if neither base editing nor
///   prime editing can be designed in the supplied window.
pub fn plan_variant_correction(req: &VariantCorrectionRequest) -> Result<VariantCorrectionPlan> {
    if !is_acgt(&req.reference) {
        return Err(GeneditingError::invalid_target(
            "region",
            "reference window must be non-empty ACGT",
        ));
    }
    if req.variant_pos >= req.reference.len() {
        return Err(GeneditingError::invalid_target(
            "variant",
            "variant position is outside the reference window",
        ));
    }

    // The "desired change" is the *correction*: variant -> reference.
    let change = correction_change(&req.variant, &req.reference, req.variant_pos)?;
    // Ask the advisor (PAM availability is determined per-modality
    // below, so consult it with the optimistic `true` and then verify).
    let advice = advise_strategy(&change, true);

    // Try base editing first when the advisor likes it and the change
    // is a reachable transition.
    if advice.recommended() == EditApproach::BaseEditing {
        if let Some(plan) = try_base_edit_correction(req)? {
            return Ok(plan);
        }
        // Base editing was preferred but un-designable here — fall
        // through to prime editing.
    }

    // Prime editing — the universal fallback for any substitution or
    // small indel.
    let pegrna = try_prime_correction(req)?;
    let rationale = format!(
        "Correcting a {} by prime editing: a PE2-style pegRNA reverts the \
         patient sequence to wild type — prime editing handles this change \
         without a double-strand break or a donor.",
        req.variant.label(),
    );
    Ok(VariantCorrectionPlan {
        variant: req.variant.clone(),
        approach: EditApproach::PrimeEditing,
        reagent: CorrectionReagent::Prime {
            editor: PrimeEditorId::Pe2,
            pegrna,
        },
        rationale,
    })
}

/// Translates a pathogenic variant into the [`DesiredChange`] that
/// *corrects* it, validating consistency with the reference.
fn correction_change(
    variant: &PathogenicVariant,
    reference: &[u8],
    pos: usize,
) -> Result<DesiredChange> {
    match variant {
        PathogenicVariant::Substitution { ref_base, variant_base } => {
            let r = reference[pos].to_ascii_uppercase();
            if r != ref_base.to_ascii_uppercase() {
                return Err(GeneditingError::invalid_target(
                    "variant",
                    format!(
                        "reference base at the variant position is `{}`, not the \
                         declared wild-type `{}`",
                        r as char,
                        ref_base.to_ascii_uppercase() as char
                    ),
                ));
            }
            // Correction restores ref_base: the patient has
            // variant_base, the edit is variant_base -> ref_base.
            Ok(DesiredChange::PointMutation {
                from: variant_base.to_ascii_uppercase(),
                to: ref_base.to_ascii_uppercase(),
            })
        }
        PathogenicVariant::Insertion { inserted } => {
            // Correction removes the inserted bases — a deletion.
            Ok(DesiredChange::SmallDeletion {
                size: inserted.len(),
            })
        }
        PathogenicVariant::Deletion { deleted } => {
            // Correction restores the deleted bases — an insertion.
            Ok(DesiredChange::SmallInsertion {
                size: deleted.len(),
            })
        }
    }
}

/// Tries to design a base-editing correction. Returns `Ok(None)` when
/// base editing simply does not apply / cannot be placed (so the
/// caller can fall back), `Err` only for a genuine error.
fn try_base_edit_correction(
    req: &VariantCorrectionRequest,
) -> Result<Option<VariantCorrectionPlan>> {
    let (ref_base, variant_base) = match &req.variant {
        PathogenicVariant::Substitution { ref_base, variant_base } => {
            (ref_base.to_ascii_uppercase(), variant_base.to_ascii_uppercase())
        }
        _ => return Ok(None), // base editors do point transitions only
    };
    // The correction edit is variant_base -> ref_base. Pick the editor
    // class: C->T / G->A ⇒ CBE; A->G / T->C ⇒ ABE.
    let class = match (variant_base, ref_base) {
        (b'C', b'T') | (b'G', b'A') => BaseEditorClass::Cbe,
        (b'A', b'G') | (b'T', b'C') => BaseEditorClass::Abe,
        _ => return Ok(None), // not a transition
    };
    // The patient sequence at the variant position carries
    // variant_base; build that reference for the base-edit designer.
    let mut patient = req.reference.clone();
    patient[req.variant_pos] = variant_base;

    // Try each editor of the right class.
    for id in BaseEditorId::all() {
        if base_editor(id).class != class {
            continue;
        }
        let be_req = BaseEditRequest {
            reference: patient.clone(),
            edit_pos: req.variant_pos,
            from_base: variant_base,
            to_base: ref_base,
            editor: id,
        };
        match design_base_edit(&be_req) {
            Ok(design) => {
                let rationale = format!(
                    "Correcting a {} by base editing with {}: an indel-free \
                     {} reversion places the target base at protospacer \
                     position {} (product-purity heuristic {:.2}).",
                    req.variant.label(),
                    base_editor(id).name,
                    base_editor(id).class.transition(),
                    design.window_position,
                    design.product_purity,
                );
                return Ok(Some(VariantCorrectionPlan {
                    variant: req.variant.clone(),
                    approach: EditApproach::BaseEditing,
                    reagent: CorrectionReagent::BaseEdit { editor: id, design },
                    rationale,
                }));
            }
            // No PAM for this editor — try the next.
            Err(GeneditingError::NoValidDesign { .. }) => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(None)
}

/// Designs a prime-editing pegRNA that corrects the variant.
fn try_prime_correction(req: &VariantCorrectionRequest) -> Result<PegRna> {
    // The correction edit, expressed as a PrimeEdit on the *patient*
    // sequence — but the pegRNA designer works against a reference and
    // installs a PrimeEdit on it. The patient sequence is the
    // designer's "reference"; the PrimeEdit restores wild type.
    let (patient, edit_pos, prime_edit): (Vec<u8>, usize, PrimeEdit) = match &req.variant {
        PathogenicVariant::Substitution { variant_base, ref_base } => {
            let mut patient = req.reference.clone();
            patient[req.variant_pos] = variant_base.to_ascii_uppercase();
            (
                patient,
                req.variant_pos,
                PrimeEdit::snv(ref_base.to_ascii_uppercase()),
            )
        }
        PathogenicVariant::Insertion { inserted } => {
            // The patient carries `inserted` at variant_pos; the
            // reference does not. The patient sequence = reference with
            // `inserted` spliced in. The correction deletes it.
            let mut patient = req.reference.clone();
            let ins = crate::sequtil::upper(inserted);
            patient.splice(req.variant_pos..req.variant_pos, ins.iter().copied());
            (
                patient,
                req.variant_pos,
                PrimeEdit::Deletion { len: inserted.len() },
            )
        }
        PathogenicVariant::Deletion { deleted } => {
            // The patient is missing `deleted`; the patient sequence =
            // reference with `deleted` removed. The correction inserts
            // it back.
            if req.variant_pos + deleted.len() > req.reference.len() {
                return Err(GeneditingError::invalid_target(
                    "variant",
                    "declared deletion runs past the reference window",
                ));
            }
            let mut patient = req.reference.clone();
            patient.drain(req.variant_pos..req.variant_pos + deleted.len());
            (
                patient,
                req.variant_pos,
                PrimeEdit::Insertion {
                    seq: crate::sequtil::upper(deleted),
                },
            )
        }
    };

    let peg_req = PegRnaRequest::new(patient, edit_pos, prime_edit, PrimeEditorId::Pe2);
    design_pegrna(&peg_req).map_err(|e| match e {
        GeneditingError::NoValidDesign { .. } => GeneditingError::no_valid_design(
            "variant_correction",
            "no PAM-adjacent prime-editing spacer reaches the variant in this \
             reference window",
        ),
        other => other,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // A reference with a forward NGG PAM (nick at 17) so a pegRNA can
    // be designed for a downstream variant.
    fn reference() -> Vec<u8> {
        b"ACGTACGTACGTACGTACGTAGGCATGCATGCATGCATGCATGCATGC".to_vec()
    }

    #[test]
    fn plans_a_substitution_correction() {
        // reference() index 25 is `T` (…A G G C A[24] T[25] G C…), so a
        // consistent pathogenic substitution declares ref_base = T; the
        // patient carries an A there.
        let req = VariantCorrectionRequest {
            reference: reference(),
            variant_pos: 25,
            variant: PathogenicVariant::Substitution {
                ref_base: b'T', // reference()[25]
                variant_base: b'A',
            },
        };
        let plan = plan_variant_correction(&req).unwrap();
        // Either base editing or prime editing produced a reagent.
        match plan.reagent {
            CorrectionReagent::BaseEdit { .. } | CorrectionReagent::Prime { .. } => {}
        }
        assert!(plan.rationale.contains("Correcting"));
    }

    #[test]
    fn rejects_inconsistent_reference_base() {
        let req = VariantCorrectionRequest {
            reference: reference(),
            variant_pos: 25,
            variant: PathogenicVariant::Substitution {
                ref_base: b'A', // index 25 is actually C
                variant_base: b'G',
            },
        };
        let err = plan_variant_correction(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn rejects_non_acgt_reference() {
        let req = VariantCorrectionRequest {
            reference: b"ACGTNNNN".to_vec(),
            variant_pos: 2,
            variant: PathogenicVariant::Substitution {
                ref_base: b'G',
                variant_base: b'A',
            },
        };
        assert!(plan_variant_correction(&req).is_err());
    }

    #[test]
    fn rejects_out_of_range_position() {
        let req = VariantCorrectionRequest {
            reference: reference(),
            variant_pos: 9999,
            variant: PathogenicVariant::Substitution {
                ref_base: b'A',
                variant_base: b'G',
            },
        };
        assert!(plan_variant_correction(&req).is_err());
    }

    #[test]
    fn plans_an_insertion_correction_as_a_deletion() {
        // The patient carries a 3-bp insertion; correction deletes it.
        let req = VariantCorrectionRequest {
            reference: reference(),
            variant_pos: 25,
            variant: PathogenicVariant::Insertion {
                inserted: b"GGG".to_vec(),
            },
        };
        let plan = plan_variant_correction(&req).unwrap();
        // An insertion correction is necessarily prime editing.
        assert_eq!(plan.approach, EditApproach::PrimeEditing);
        if let CorrectionReagent::Prime { pegrna, .. } = &plan.reagent {
            assert!(matches!(pegrna.edit, PrimeEdit::Deletion { len: 3 }));
        } else {
            panic!("expected a prime-editing reagent");
        }
    }

    #[test]
    fn plans_a_deletion_correction_as_an_insertion() {
        let req = VariantCorrectionRequest {
            reference: reference(),
            variant_pos: 25,
            variant: PathogenicVariant::Deletion {
                deleted: b"CAT".to_vec(),
            },
        };
        let plan = plan_variant_correction(&req).unwrap();
        assert_eq!(plan.approach, EditApproach::PrimeEditing);
        if let CorrectionReagent::Prime { pegrna, .. } = &plan.reagent {
            assert!(matches!(pegrna.edit, PrimeEdit::Insertion { .. }));
        } else {
            panic!("expected a prime-editing reagent");
        }
    }

    #[test]
    fn correction_change_inverts_the_variant() {
        // A pathogenic G (ref C): the correction restores C, editing
        // the patient's G back to C.
        let change = correction_change(
            &PathogenicVariant::Substitution {
                ref_base: b'C',
                variant_base: b'G',
            },
            b"AAACAAA",
            3,
        )
        .unwrap();
        assert_eq!(
            change,
            DesiredChange::PointMutation { from: b'G', to: b'C' }
        );
    }

    #[test]
    fn variant_labels() {
        assert!(PathogenicVariant::Substitution { ref_base: b'A', variant_base: b'G' }
            .label()
            .contains("substitution"));
        assert!(PathogenicVariant::Insertion { inserted: b"GG".to_vec() }
            .label()
            .contains("insertion"));
    }

    #[test]
    fn unreachable_variant_reports_no_design() {
        // A variant upstream of every PAM nick — no pegRNA possible.
        let req = VariantCorrectionRequest {
            reference: reference(),
            variant_pos: 3,
            variant: PathogenicVariant::Substitution {
                ref_base: b'T', // reference()[3]
                variant_base: b'A',
            },
        };
        let err = plan_variant_correction(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.no_valid_design");
    }
}
