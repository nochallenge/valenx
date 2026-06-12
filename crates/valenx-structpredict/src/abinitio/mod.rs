//! Ab-initio / fragment-assembly structure prediction — Rosetta-class.
//!
//! **Roadmap features 7-12.** When no template exists, structure must
//! be predicted from the sequence alone. The classical answer — the
//! one Rosetta `AbinitioRelax` made famous — is **fragment
//! assembly**: instead of moving individual atoms, the protocol
//! repeatedly swaps in short backbone *fragments* (3- and 9-residue
//! pieces of real protein geometry) and keeps the moves that lower a
//! knowledge-based energy.
//!
//! The pipeline:
//!
//! 1. **Fragment library** ([`fragments`]) — short backbone fragments
//!    keyed by a sequence window and the predicted secondary
//!    structure.
//! 2. **Secondary-structure prediction** ([`ss`]) — a classical
//!    GOR / Chou-Fasman predictor (no neural network).
//! 3. **Knowledge-based scoring** ([`score`]) — distance-dependent
//!    statistical potentials plus contact and torsion terms, the
//!    Rosetta-centroid-style energy.
//! 4. **Monte-Carlo fragment assembly** ([`assemble`]) — simulated
//!    annealing with fragment-insertion moves.
//! 5. **Coarse-to-fine protocol** ([`protocol`]) — centroid assembly
//!    then an all-atom pass.
//! 6. **Decoy selection** ([`decoy`]) — generate many models,
//!    cluster them, and pick a representative.
//!
//! **Honest accuracy note.** Classical fragment assembly samples
//! plausible folds and, for small single-domain proteins, sometimes
//! finds the right one — but it has no learnt co-evolutionary signal
//! and is *not* AlphaFold-accurate. It is a real, historically
//! important method; treat its models as low-resolution hypotheses.

pub mod assemble;
pub mod decoy;
pub mod dope;
pub mod fragments;
pub mod protocol;
pub mod score;
pub mod ss;

pub use assemble::{fragment_assembly, AssemblyOptions, AssemblyResult, AssemblyScorer};
pub use decoy::{cluster_decoys, select_model, DecoyCluster};
pub use dope::{dope_energy, dope_score, dope_table, DopePair, DopeScore, DopeWeights};
pub use fragments::{
    build_fragment_library, fragment_class_count, fragment_class_labels, Fragment, FragmentLibrary,
};
pub use protocol::{coarse_to_fine, coarse_to_fine_with_refine, ProtocolResult};
pub use score::{score_model, KnowledgeScore, ScoreWeights};
pub use ss::{predict_secondary_structure, SecondaryStructure, SsPrediction};
