//! **Feature 3 — backbone transfer.**
//!
//! Given a target-template alignment and the template's structure,
//! copy the template's backbone coordinates onto the target for every
//! equivalenced residue. This is the step that turns an alignment
//! into 3-D coordinates: a comparative model's core assumption is
//! that homologous residues occupy homologous positions, so the
//! template backbone *is* the model backbone for the aligned regions.
//!
//! The result is a [`crate::model::ProteinModel`] in which the
//! equivalenced residues carry the template's `N`/`CA`/`C`/`O`/`CB`
//! atoms and the non-equivalenced residues are gaps for the loop
//! modeller to fill.

use crate::error::{Result, StructPredictError};
use crate::homology::align::TargetTemplateAlignment;
use crate::model::{ModelResidue, ProteinModel};

/// Transfers the template backbone onto the target.
///
/// `target_sequence` is the target's one-letter sequence (defines the
/// model's residue identities and length); `alignment` is the
/// target-template alignment; `template` is the template structure;
/// `template_chain` names the chain to copy.
///
/// Every equivalenced target residue receives the corresponding
/// template residue's backbone atoms; non-equivalenced residues are
/// left as gaps.
///
/// # Errors
/// [`StructPredictError::Invalid`] if `target_sequence` length
/// disagrees with the alignment's target row;
/// [`StructPredictError::NotFound`] if the template chain is absent.
pub fn transfer_backbone(
    target_sequence: &str,
    alignment: &TargetTemplateAlignment,
    template: &valenx_biostruct::Structure,
    template_chain: &str,
) -> Result<ProteinModel> {
    let target_sequence = target_sequence.trim();
    if target_sequence.is_empty() {
        return Err(StructPredictError::invalid("target_sequence", "empty"));
    }
    // The number of residues in the gapped target row.
    let row_residues = alignment
        .target_row
        .iter()
        .filter(|&&b| b != b'-')
        .count();
    if row_residues != target_sequence.len() {
        return Err(StructPredictError::invalid(
            "target_sequence",
            format!(
                "length {} disagrees with alignment target row ({} residues)",
                target_sequence.len(),
                row_residues
            ),
        ));
    }

    let mut model = ProteinModel::from_sequence(target_sequence)?;

    // The template's amino-acid residues, in chain order — the
    // template index in the alignment is an index into this list.
    let tmpl_model = template.first_model();
    let chain = tmpl_model.chain(template_chain).ok_or_else(|| {
        StructPredictError::not_found("chain", format!("template chain `{template_chain}`"))
    })?;
    let tmpl_residues: Vec<&valenx_biostruct::Residue> = chain
        .residues
        .iter()
        .filter(|r| r.is_amino_acid())
        .collect();

    for &(ti, pi) in &alignment.equivalences {
        let (Some(model_res), Some(&tmpl_res)) =
            (model.residues.get_mut(ti), tmpl_residues.get(pi))
        else {
            continue;
        };
        copy_backbone(tmpl_res, model_res);
    }

    Ok(model)
}

/// Copies the `N`/`CA`/`C`/`O`/`CB` atom coordinates of a template
/// residue into a model residue.
fn copy_backbone(src: &valenx_biostruct::Residue, dst: &mut ModelResidue) {
    let pos = |name: &str| src.atom(name).map(|a| a.coord);
    dst.n = pos("N");
    dst.ca = pos("CA");
    dst.c = pos("C");
    dst.o = pos("O");
    dst.cb = pos("CB");
}

/// Counts the equivalenced residues whose backbone was actually
/// transferred (the template residue carried the backbone atoms).
pub fn transferred_count(model: &ProteinModel) -> usize {
    model.residues.iter().filter(|r| r.has_backbone()).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::homology::align::target_template_alignment;

    const TEMPLATE_PDB: &str = "\
ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1      2.000   1.400   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      1.300   2.400   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1      2.000  -0.800   1.200  1.00  0.00           C
ATOM      6  N   GLY A   2      3.300   1.400   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2      4.000   2.700   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2      5.500   2.500   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2      6.000   1.400   0.000  1.00  0.00           O
ATOM     10  N   CYS A   3      6.200   3.600   0.000  1.00  0.00           N
ATOM     11  CA  CYS A   3      7.650   3.700   0.000  1.00  0.00           C
ATOM     12  C   CYS A   3      8.100   5.100   0.000  1.00  0.00           C
ATOM     13  O   CYS A   3      7.300   6.000   0.000  1.00  0.00           O
ATOM     14  CB  CYS A   3      8.200   2.900   1.200  1.00  0.00           C
END
";

    #[test]
    fn transfers_template_backbone() {
        let s = valenx_biostruct::read_structure(TEMPLATE_PDB, "tmpl").expect("parse");
        // Target == template sequence "AGC" → full equivalence.
        let aln = target_template_alignment("AGC", "AGC", &[]).expect("align");
        let model = transfer_backbone("AGC", &aln, &s, "A").expect("transfer");
        assert_eq!(model.len(), 3);
        assert_eq!(transferred_count(&model), 3);
        assert!(model.is_complete());
        // The first residue's Cα matches the template.
        let ca = model.residues[0].ca.expect("ca");
        assert!((ca.x - 1.458).abs() < 1e-6);
    }

    #[test]
    fn unaligned_residues_become_gaps() {
        let s = valenx_biostruct::read_structure(TEMPLATE_PDB, "tmpl").expect("parse");
        // Target has an extra residue the 3-residue template lacks.
        let aln = target_template_alignment("AGWC", "AGC", &[]).expect("align");
        let model = transfer_backbone("AGWC", &aln, &s, "A").expect("transfer");
        assert_eq!(model.len(), 4);
        // One residue (the insertion) has no backbone.
        assert_eq!(transferred_count(&model), 3);
        assert!(!model.gaps().is_empty());
    }

    #[test]
    fn missing_chain_rejected() {
        let s = valenx_biostruct::read_structure(TEMPLATE_PDB, "tmpl").expect("parse");
        let aln = target_template_alignment("AGC", "AGC", &[]).expect("align");
        assert!(transfer_backbone("AGC", &aln, &s, "Z").is_err());
    }
}
