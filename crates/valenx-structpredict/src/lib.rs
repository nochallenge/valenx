//! # valenx-structpredict — classical protein structure prediction, design & cryo-EM
//!
//! This crate delivers
//! **protein structure prediction, physics-based protein design and
//! cryo-EM reconstruction** — but with an important, honestly-stated
//! constraint.
//!
//! ## Classical methods only — and what that means
//!
//! AlphaFold, ESMFold, RoseTTAFold and RFdiffusion are **trained
//! neural networks**. Valenx ships those as *adapter hooks only* (the
//! project's standing "no llms" rule — a reimplemented network is
//! useless without its trained weights anyway). This crate instead
//! implements the **classical, physics-and-statistics-based
//! algorithms that predate and parallel the deep-learning era** and
//! that need *no trained weights*:
//!
//! - **Comparative / homology modelling** — the Modeller-class
//!   approach: find a homologous template by sequence alignment,
//!   transfer its backbone, close the loops, place the sidechains,
//!   and satisfy spatial restraints by optimisation.
//! - **Ab-initio fragment assembly** — the Rosetta-class approach:
//!   build short backbone-fragment libraries, score conformations
//!   with knowledge-based statistical potentials, and assemble a fold
//!   by simulated-annealing Monte-Carlo fragment insertion.
//! - **Physics-based protein design** — the Rosetta `fixbb`-class
//!   approach: search amino-acid sequence + rotamer space to minimise
//!   a statistical/physics energy on a fixed backbone.
//! - **Cryo-EM reconstruction** — the classical signal-processing
//!   core of RELION / EMAN: CTF modelling, particle picking, 2D class
//!   averaging, weighted-back-projection 3D reconstruction,
//!   projection-matching refinement and FSC resolution estimation.
//!
//! **These are real, genuinely useful algorithms** — comparative
//! modelling and fragment assembly were the field's workhorses for
//! two decades and remain useful when a good template exists or for
//! quick low-resolution models. **They are NOT AlphaFold-accuracy.**
//! AlphaFold-class accuracy comes specifically from a trained
//! network's learnt co-evolutionary signal; that accuracy is *by
//! design* out of scope here. When you need it, use the
//! AlphaFold / RoseTTAFold subprocess adapters; for a classical
//! model with no GPU and no weights, use this crate.
//!
//! ## What it does
//!
//! - **Homology modelling** ([`homology`]) — template search,
//!   target-template alignment, backbone transfer, CCD loop closure,
//!   rotamer sidechain placement, spatial-restraint assembly.
//! - **Ab-initio** ([`abinitio`]) — fragment libraries, GOR /
//!   Chou-Fasman secondary-structure prediction, knowledge-based
//!   statistical potentials, Monte-Carlo fragment assembly, a
//!   centroid→all-atom protocol, decoy clustering and selection.
//! - **Refinement** ([`refine`]) — DEE + simulated-annealing rotamer
//!   repacking, gradient energy minimisation (via the [`valenx_md`]
//!   force-field minimisers), Ramachandran backbone refinement, a
//!   clash / Ramachandran / packing quality score, and Kabsch
//!   superposition with RMSD / GDT.
//! - **Design** ([`design`]) — fixed-backbone sequence design, the
//!   combinatorial rotamer design search, a design-scoring potential,
//!   loop / motif grafting, interface design.
//! - **Cryo-EM** ([`cryoem`]) — MRC I/O + a particle-stack model, a
//!   CTF model + estimation, particle picking, 2D class averaging, 3D
//!   weighted back-projection, projection-matching refinement, FSC.
//! - **Drivers** ([`driver`]) — the top-level
//!   [`predict_homology`], [`predict_abinitio`],
//!   [`design_sequence`] and [`reconstruct_cryoem`] entry points
//!   plus a [`StructPredictReport`].
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, StructPredictError>`](error::StructPredictError). The
//! error type carries stable [`code`](error::StructPredictError::code)
//! and [`category`](error::StructPredictError::category) accessors for
//! telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with Modeller /
//! Rosetta / RELION. Each module documents its own simplifications
//! inline; the post-depth-pass status:
//!
//! - **Statistical potential** — the production default is now the
//!   DOPE-class distance-dependent atomic statistical potential
//!   (`E_ab(d) = −kT·ln(g_obs/g_ref)`) shipped in
//!   [`abinitio::dope`]. It encodes the published functional form
//!   over the canonical Cα-Cα, Cβ-Cβ and hydrophobic-Cα-Cα atom-pair
//!   tables, with the standard 0.5 Å bin width / 15 Å cutoff. Full
//!   Modeller-DOPE 158-atom-type coverage stays adapter-only T3.
//! - **Fragment libraries** — the production default is the
//!   PDB-curated-style library in [`abinitio::fragments`] with 14
//!   canonical backbone classes (α-interior / N-cap / C-cap,
//!   3₁₀-helix, π-helix, β-strand interior / edge, β-turn types
//!   I/II/I'/II', γ-turn classic / inverse, PPII) — published
//!   per-residue (φ, ψ, ω) means with published one-σ Ramachandran
//!   spreads. Full Rosetta 10⁴-PDB-chain mining stays adapter-only.
//! - **Refinement** — the production path adds a DOPE-driven
//!   simulated-annealing fragment-insertion MC refinement loop
//!   ([`refine::mc_refine`]) with a temperature schedule + per-cycle
//!   Ramachandran cleanup; the classical decoy generation +
//!   clustering pick lives in [`abinitio::cluster_decoys`].
//! - **Other** — CCD loop closure and the Monte-Carlo annealer are
//!   the genuine textbook algorithms; rotamer libraries are a compact
//!   backbone-independent set (not the Dunbrack backbone-dependent
//!   library); the cryo-EM reconstruction is real weighted
//!   back-projection / direct Fourier inversion on a discrete grid
//!   (not the regularised-likelihood RELION pipeline); FSC + the
//!   gold-standard criterion are exact. Nothing in this crate is a
//!   neural network.

#![forbid(unsafe_code)]

pub mod aa;
pub mod abinitio;
pub mod cryoem;
pub mod design;
pub mod driver;
pub mod error;
pub mod homology;
pub mod model;
pub mod refine;
pub mod rotamer;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, Result, StructPredictError};
pub use model::{ModelResidue, ProteinModel};

pub use driver::{
    design_sequence, predict_abinitio, predict_homology, reconstruct_cryoem, StructPredictReport,
};

// --- Re-exports of the new commercial-depth surface --------------------

pub use abinitio::{
    coarse_to_fine_with_refine, dope_energy, dope_score, fragment_class_count,
    fragment_class_labels, AssemblyScorer, DopePair, DopeScore, DopeWeights,
};
pub use refine::{mc_refine, McRefineOptions, McRefineResult};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_error_reexported() {
        let e = StructPredictError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn protein_model_reexported() {
        let m = ProteinModel::from_sequence("ACDEFG").expect("model");
        assert_eq!(m.len(), 6);
    }
}
