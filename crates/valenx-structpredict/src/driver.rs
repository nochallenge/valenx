//! **Feature 30 — top-level drivers + report.**
//!
//! This module is the crate's typed, single-call surface — the entry
//! points an MCP / LLM tool layer (or any caller) drives without
//! threading the individual algorithm modules together by hand:
//!
//! - [`predict_homology`] — full comparative-modelling pipeline:
//!   template search → alignment → backbone transfer → loop closure →
//!   sidechain placement → restraint satisfaction → refinement.
//! - [`predict_abinitio`] — full ab-initio pipeline: fragment library
//!   → Monte-Carlo fragment assembly → centroid-to-all-atom
//!   refinement.
//! - [`design_sequence`] — fixed-backbone protein design on a given
//!   structure.
//! - [`reconstruct_cryoem`] — classical cryo-EM 3-D reconstruction
//!   from a set of oriented projections, with an FSC resolution
//!   estimate.
//!
//! Every driver returns a [`StructPredictReport`] — a single typed
//! struct carrying the result plus the honest provenance and
//! quality numbers a caller needs to judge it.

use serde::{Deserialize, Serialize};

use crate::abinitio::protocol::coarse_to_fine;
use crate::cryoem::fsc::gold_standard_resolution;
use crate::cryoem::mrc::Volume3d;
use crate::cryoem::reconstruct::{reconstruct_3d, Projection};
use crate::design::fixbb::design_fixed_backbone;
use crate::design::search::ResiduePalette;
use crate::error::{Result, StructPredictError};
use crate::homology::align::target_template_alignment;
use crate::homology::loops::model_loops;
use crate::homology::restraints::{derive_restraints, satisfy_restraints};
use crate::homology::search::{search_templates, TemplateCandidate};
use crate::homology::sidechains::place_sidechains;
use crate::homology::transfer::transfer_backbone;
use crate::model::ProteinModel;
use crate::refine::quality::{assess_quality, QualityReport};

/// Which kind of job a [`StructPredictReport`] describes.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum JobKind {
    /// Comparative / homology modelling.
    Homology,
    /// Ab-initio fragment-assembly prediction.
    AbInitio,
    /// Fixed-backbone protein design.
    Design,
    /// Cryo-EM 3-D reconstruction.
    CryoEm,
}

/// The unified report returned by every top-level driver.
///
/// Not every field is meaningful for every job — a cryo-EM job leaves
/// the protein-model fields `None`, a design job leaves the volume
/// `None`. `notes` always carries an honest plain-text summary.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructPredictReport {
    /// Which kind of job produced this report.
    pub kind: JobKind,
    /// The predicted / designed protein model (homology, ab-initio,
    /// design jobs). `None` for a cryo-EM job.
    pub model: Option<ProteinModel>,
    /// The model as PDB `ATOM` records (when [`Self::model`] is set).
    pub pdb: Option<String>,
    /// The designed sequence (design jobs) — also the model's
    /// sequence for the others.
    pub sequence: Option<String>,
    /// A reference-free model-quality report (protein jobs).
    pub quality: Option<QualityReport>,
    /// The reconstructed density map (cryo-EM jobs).
    pub volume: Option<Volume3d>,
    /// The gold-standard FSC resolution in ångström (cryo-EM jobs
    /// supplied with two half-maps' projections).
    pub resolution: Option<f64>,
    /// An honest plain-text provenance + caveat summary.
    pub notes: String,
}

impl StructPredictReport {
    /// `true` when the job produced a usable primary result.
    pub fn is_ok(&self) -> bool {
        self.model.is_some() || self.volume.is_some()
    }
}

/// Runs the full comparative / homology-modelling pipeline.
///
/// `query` is the target sequence; `templates` is a set of candidate
/// templates (each a sequence + its structure text). The driver picks
/// the best template, transfers its backbone, closes loops, places
/// sidechains, satisfies template-derived spatial restraints, and
/// scores the model.
///
/// `template_chain` names the chain to use from the chosen template's
/// structure.
///
/// # Errors
/// [`StructPredictError::NotFound`] if no template is usable;
/// propagates pipeline-stage errors.
pub fn predict_homology(
    query: &str,
    templates: &[TemplateCandidate],
    template_chain: &str,
) -> Result<StructPredictReport> {
    let query = query.trim();
    if query.is_empty() {
        return Err(StructPredictError::invalid("query", "empty sequence"));
    }
    // 1. Template search.
    let hits = search_templates(query, templates)?;
    let best = &hits[0];
    let template = &templates[best.index];

    // 2. Parse the template structure and align.
    let structure = valenx_biostruct::read_structure(&template.structure_text, &template.id)
        .map_err(|e| StructPredictError::parse("pdb", e.to_string()))?;
    let alignment = target_template_alignment(query, &template.sequence, &[])?;

    // 3. Transfer the backbone.
    let mut model = transfer_backbone(query, &alignment, &structure, template_chain)?;

    // 4. Close the loops the template did not cover.
    let closures = model_loops(&mut model, 300, 0.3)?;

    // 5. Place sidechains.
    place_sidechains(&mut model)?;

    // 6. Satisfy template-derived spatial restraints (Cα level).
    // Build a template-coordinate model to measure restraints from.
    let template_model = ProteinModel::from_structure_chain(&structure, template_chain)?;
    let restraints = derive_restraints(&template_model, &alignment, 10.0)?;
    if !restraints.is_empty() {
        satisfy_restraints(&mut model, &restraints, 150)?;
        // Restraint satisfaction moved the Cα atoms — re-place Cβ.
        place_sidechains(&mut model)?;
    }

    // 7. Quality.
    let quality = assess_quality(&model).ok();

    let notes = format!(
        "Comparative model from template `{}` ({:.0}% identity, {:.0}% coverage). \
         {} loop region(s) built by CCD; {} template-derived spatial restraint(s) applied. \
         Classical homology modelling — reliable with a close template, NOT AlphaFold-accuracy.",
        best.id,
        best.identity * 100.0,
        best.coverage * 100.0,
        closures.len(),
        restraints.len()
    );

    Ok(StructPredictReport {
        kind: JobKind::Homology,
        pdb: Some(model.to_pdb()),
        sequence: Some(model.sequence()),
        model: Some(model),
        quality,
        volume: None,
        resolution: None,
        notes,
    })
}

/// Runs the full ab-initio fragment-assembly prediction pipeline.
///
/// `sequence` is the target. The driver builds a fragment library,
/// runs the centroid-to-all-atom coarse-to-fine protocol, and reports
/// the model with its quality.
///
/// `centroid_moves` / `repack_moves` size the two Monte-Carlo
/// searches; `seed` fixes every RNG so the prediction is
/// reproducible.
///
/// # Errors
/// [`StructPredictError::Invalid`] for too short a sequence;
/// propagates protocol errors.
pub fn predict_abinitio(
    sequence: &str,
    centroid_moves: usize,
    repack_moves: usize,
    seed: u64,
) -> Result<StructPredictReport> {
    let protocol = coarse_to_fine(sequence, centroid_moves, repack_moves, seed)?;
    let notes = format!(
        "Ab-initio model by Monte-Carlo fragment assembly (centroid score {:.2}; \
         {} Ramachandran outlier(s) removed in the all-atom pass). \
         Classical fragment assembly samples plausible folds but has no learnt \
         co-evolutionary signal — treat as a low-resolution hypothesis, NOT \
         AlphaFold-accurate.",
        protocol.centroid_score, protocol.outliers_removed
    );
    Ok(StructPredictReport {
        kind: JobKind::AbInitio,
        pdb: Some(protocol.model.to_pdb()),
        sequence: Some(protocol.model.sequence()),
        model: Some(protocol.model),
        quality: Some(protocol.quality),
        volume: None,
        resolution: None,
        notes,
    })
}

/// Runs fixed-backbone protein design on a structure.
///
/// `structure_text` is a PDB / mmCIF coordinate file; `chain` names
/// the chain whose backbone to redesign. The driver designs a new
/// sequence for that backbone with the classical physics-based design
/// search.
///
/// `moves` sizes the Monte-Carlo design search; `seed` fixes the RNG.
///
/// # Errors
/// [`StructPredictError::Parse`] for a malformed structure;
/// propagates design errors.
pub fn design_sequence(
    structure_text: &str,
    chain: &str,
    moves: usize,
    seed: u64,
) -> Result<StructPredictReport> {
    let structure = valenx_biostruct::read_structure(structure_text, "design_input")
        .map_err(|e| StructPredictError::parse("pdb", e.to_string()))?;
    let model = ProteinModel::from_structure_chain(&structure, chain)?;
    let palette = ResiduePalette::unrestricted(model.residues.len());
    let fixbb = design_fixed_backbone(&model, &palette, moves, seed)?;

    // Build the redesigned model (same backbone, new identities).
    let mut designed = model.clone();
    for (res, aa) in designed
        .residues
        .iter_mut()
        .zip(fixbb.designed_sequence.chars())
    {
        res.aa = aa;
    }
    let quality = assess_quality(&designed).ok();

    let notes = format!(
        "Fixed-backbone design on chain `{}` ({} residues, {:.0}% mutated from the \
         input sequence; design energy {:.2}). Classical physics-based design with a \
         knowledge-based score — real, but lower success rate on a hard backbone than \
         a trained network (ProteinMPNN); use the ProteinMPNN adapter for that.",
        chain,
        designed.len(),
        fixbb.mutation_fraction * 100.0,
        fixbb.design_energy
    );

    Ok(StructPredictReport {
        kind: JobKind::Design,
        pdb: Some(designed.to_pdb()),
        sequence: Some(fixbb.designed_sequence),
        model: Some(designed),
        quality,
        volume: None,
        resolution: None,
        notes,
    })
}

/// Runs classical cryo-EM 3-D reconstruction from oriented
/// projections.
///
/// `projections` is a set of 2-D projection images with known
/// orientations; `box_size` is their (square) size. The driver
/// reconstructs a 3-D density map by weighted back-projection.
///
/// When `half_maps` is supplied — a pair of independent
/// half-projection sets — the driver also reconstructs the two
/// half-maps and reports the gold-standard FSC resolution.
///
/// # Errors
/// Propagates reconstruction and FSC errors.
pub fn reconstruct_cryoem(
    projections: &[Projection],
    box_size: usize,
    half_maps: Option<(&[Projection], &[Projection])>,
) -> Result<StructPredictReport> {
    let recon = reconstruct_3d(projections, box_size)?;

    // Optional gold-standard resolution from two half-maps.
    let resolution = match half_maps {
        Some((a, b)) => {
            let ra = reconstruct_3d(a, box_size)?;
            let rb = reconstruct_3d(b, box_size)?;
            gold_standard_resolution(&ra.volume, &rb.volume).ok()
        }
        None => None,
    };

    let notes = format!(
        "Cryo-EM 3-D map reconstructed by weighted back-projection from {} oriented \
         projection(s){}. Classical signal-processing reconstruction — real \
         back-projection geometry; NOT the regularised-likelihood RELION pipeline, and \
         the modern neural particle pickers stay adapter-only.",
        recon.projections_used,
        match resolution {
            Some(r) => format!("; gold-standard FSC resolution ~{r:.1} Å"),
            None => String::new(),
        }
    );

    Ok(StructPredictReport {
        kind: JobKind::CryoEm,
        model: None,
        pdb: None,
        sequence: None,
        quality: None,
        volume: Some(recon.volume),
        resolution,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryoem::mrc::Image2d;

    const MINI_TEMPLATE: &str = "\
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
ATOM     15  N   ASP A   4      9.400   5.300   0.000  1.00  0.00           N
ATOM     16  CA  ASP A   4     10.100   6.600   0.000  1.00  0.00           C
ATOM     17  C   ASP A   4     11.600   6.400   0.000  1.00  0.00           C
ATOM     18  O   ASP A   4     12.100   5.300   0.000  1.00  0.00           O
ATOM     19  CB  ASP A   4     10.300   7.400   1.300  1.00  0.00           C
END
";

    #[test]
    fn homology_driver_builds_a_model() {
        let templates = vec![TemplateCandidate::new("1xyz", "AGCD", MINI_TEMPLATE)];
        // Query identical to the template → a clean transfer.
        let report = predict_homology("AGCD", &templates, "A").expect("homology");
        assert_eq!(report.kind, JobKind::Homology);
        assert!(report.is_ok());
        let model = report.model.as_ref().expect("model");
        assert_eq!(model.len(), 4);
        assert!(report.pdb.as_ref().unwrap().contains("ATOM"));
        assert!(report.notes.contains("1xyz"));
    }

    #[test]
    fn abinitio_driver_predicts_a_fold() {
        let report = predict_abinitio("EEEEAAAALLLLEEEE", 300, 150, 5).expect("abinitio");
        assert_eq!(report.kind, JobKind::AbInitio);
        let model = report.model.as_ref().expect("model");
        assert!(model.is_complete());
        assert!(report.quality.is_some());
        assert!(report.notes.contains("NOT"));
    }

    #[test]
    fn design_driver_redesigns_a_backbone() {
        let report = design_sequence(MINI_TEMPLATE, "A", 300, 9).expect("design");
        assert_eq!(report.kind, JobKind::Design);
        let seq = report.sequence.as_ref().expect("sequence");
        assert_eq!(seq.len(), 4); // 4-residue template
        assert!(report.notes.contains("ProteinMPNN"));
    }

    #[test]
    fn cryoem_driver_reconstructs_a_map() {
        // A handful of trivial projections of an 8³ box.
        let projections: Vec<Projection> = (0..8)
            .map(|k| {
                let mut img = Image2d::zeros(8, 8);
                // Bright centre.
                img.data[4 * 8 + 4] = 1.0;
                let tilt = std::f64::consts::PI * k as f64 / 8.0;
                Projection::new(img, 0.0, tilt, 0.0)
            })
            .collect();
        let report = reconstruct_cryoem(&projections, 8, None).expect("cryoem");
        assert_eq!(report.kind, JobKind::CryoEm);
        assert!(report.volume.is_some());
        assert!(report.model.is_none());
    }

    #[test]
    fn empty_query_rejected() {
        assert!(predict_homology("", &[], "A").is_err());
    }
}
