//! **Feature 18 — fixed-backbone design (Rosetta `fixbb`-class).**
//!
//! This is the headline protein-design entry point: given a backbone
//! and nothing else, return the sequence that best fits it. It is the
//! `fixbb` ("fixed backbone") protocol — the inverse of structure
//! prediction.
//!
//! The pipeline:
//!
//! 1. Predict the backbone's **secondary structure from its
//!    geometry** (a φ/ψ-based assignment) so the design score can use
//!    the SS-propensity term.
//! 2. Run the **combinatorial design search**
//!    ([`crate::design::combinatorial_design`]) — simulated annealing
//!    over sequence + rotamer space.
//! 3. Report the designed sequence, its score, and how much it
//!    improved on the starting sequence.
//!
//! The result is a real classical fixed-backbone design — the
//! sequence with the lowest knowledge-based design energy on the
//! given backbone. See the [`crate::design`] module note on how this
//! compares to ProteinMPNN.

use serde::{Deserialize, Serialize};

use crate::abinitio::ss::SecondaryStructure;
use crate::design::score::DesignScoreWeights;
use crate::design::search::{combinatorial_design, DesignSearchResult, ResiduePalette};
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;
use crate::refine::ramachandran::{is_allowed, model_phi_psi};

/// The outcome of a fixed-backbone design run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FixbbResult {
    /// The designed amino-acid sequence.
    pub designed_sequence: String,
    /// The design score's weighted total (lower is better).
    pub design_energy: f64,
    /// The starting sequence's design energy, for comparison.
    pub starting_energy: f64,
    /// The backbone's geometry-derived secondary structure (`HEC`).
    pub secondary_structure: String,
    /// Fraction of positions whose designed residue differs from the
    /// model's starting residue.
    pub mutation_fraction: f64,
}

/// Assigns three-state secondary structure to a backbone from its
/// φ/ψ geometry.
///
/// A residue whose (φ, ψ) lies in the helical basin is `H`, in the
/// strand basin `E`, else `C`. This is a fast, classical geometric
/// SS assignment — coarser than DSSP (which uses hydrogen bonds) but
/// adequate to drive the design score's SS-propensity term.
pub fn assign_ss_from_geometry(model: &ProteinModel) -> Vec<SecondaryStructure> {
    if !model.is_complete() {
        return vec![SecondaryStructure::Coil; model.residues.len()];
    }
    let pp = model_phi_psi(model);
    pp.iter()
        .enumerate()
        .map(|(i, &(phi, psi))| {
            if i == 0 || i + 1 == pp.len() {
                return SecondaryStructure::Coil;
            }
            if !is_allowed(phi, psi) {
                return SecondaryStructure::Coil;
            }
            // Helical basin: φ ≈ -60, ψ ≈ -45.
            if (-100.0..-30.0).contains(&phi) && (-80.0..10.0).contains(&psi) {
                SecondaryStructure::Helix
            } else if phi < -90.0 && psi > 90.0 {
                SecondaryStructure::Strand
            } else {
                SecondaryStructure::Coil
            }
        })
        .collect()
}

/// Designs a sequence for a fixed backbone.
///
/// `model` supplies the backbone. `palette` restricts the
/// per-position amino-acid choices — pass
/// [`ResiduePalette::unrestricted`] for full design. `moves` is the
/// Monte-Carlo budget; `seed` fixes the RNG.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a backbone with fewer than 2
/// Cα atoms, a palette length mismatch, or `moves == 0`.
pub fn design_fixed_backbone(
    model: &ProteinModel,
    palette: &ResiduePalette,
    moves: usize,
    seed: u64,
) -> Result<FixbbResult> {
    if model.ca_trace().len() < 2 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 2 Cα atoms to design",
        ));
    }
    let starting_sequence = model.sequence();
    let ss = assign_ss_from_geometry(model);
    let ss_string: String = ss.iter().map(|s| s.code()).collect();

    let weights = DesignScoreWeights::default();
    let search: DesignSearchResult =
        combinatorial_design(model, palette, &ss, weights, moves, seed)?;

    // Starting-sequence energy (for the improvement report).
    let starting_energy = search.initial_total;

    // Mutation fraction vs the starting sequence.
    let mutated = starting_sequence
        .chars()
        .zip(search.sequence.chars())
        .filter(|(a, b)| a != b)
        .count();
    let mutation_fraction = if starting_sequence.is_empty() {
        0.0
    } else {
        mutated as f64 / starting_sequence.len() as f64
    };

    Ok(FixbbResult {
        designed_sequence: search.sequence,
        design_energy: search.score.total,
        starting_energy,
        secondary_structure: ss_string,
        mutation_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::build_backbone_from_torsions;

    fn helix_backbone(n: usize) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        build_backbone_from_torsions(&mut m, &vec![(-63.0, -42.0); n]).expect("build");
        m
    }

    #[test]
    fn design_improves_on_the_starting_sequence() {
        let m = helix_backbone(16);
        let palette = ResiduePalette::unrestricted(m.residues.len());
        let res = design_fixed_backbone(&m, &palette, 700, 11).expect("design");
        assert!(
            res.design_energy <= res.starting_energy,
            "energy {} -> {}",
            res.starting_energy,
            res.design_energy
        );
        assert_eq!(res.designed_sequence.len(), m.residues.len());
    }

    #[test]
    fn helical_backbone_is_assigned_helix() {
        let m = helix_backbone(16);
        let ss = assign_ss_from_geometry(&m);
        let helix = ss
            .iter()
            .filter(|&&s| s == SecondaryStructure::Helix)
            .count();
        // The interior of an α-helix backbone should read as helix.
        assert!(helix > ss.len() / 2, "helix residues {helix}/{}", ss.len());
    }

    #[test]
    fn fixed_palette_yields_zero_mutations() {
        let m = helix_backbone(8);
        let palette = ResiduePalette::fixed_to(&m.sequence());
        let res = design_fixed_backbone(&m, &palette, 100, 1).expect("design");
        assert_eq!(res.mutation_fraction, 0.0);
        assert_eq!(res.designed_sequence, m.sequence());
    }

    #[test]
    fn too_short_rejected() {
        let m = ProteinModel::from_sequence("A").expect("model");
        let palette = ResiduePalette::unrestricted(1);
        assert!(design_fixed_backbone(&m, &palette, 10, 0).is_err());
    }
}
