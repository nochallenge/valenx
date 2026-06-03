//! Multiple-sequence alignment (Clustal / MUSCLE / MAFFT-class).
//!
//! The MSA pipeline:
//!
//! 1. [`guidetree`] — a pairwise [`guidetree::distance_matrix`] and a
//!    [`guidetree::GuideTree`] by [`guidetree::upgma`] or
//!    [`guidetree::neighbor_joining`].
//! 2. [`profile`] — the [`profile::Profile`] / PSSM column model and
//!    profile-profile / profile-sequence alignment.
//! 3. [`progressive`] — [`progressive::align`]: the guide-tree-driven
//!    progressive merge into an [`progressive::Msa`].
//! 4. [`mod@refine`] — [`refine::refine()`]: MUSCLE-class iterative
//!    refinement (partition, realign, keep-if-improved).
//! 5. [`analysis`] — column conservation, Shannon entropy and the
//!    consensus sequence.

pub mod analysis;
pub mod guidetree;
pub mod profile;
pub mod progressive;
pub mod refine;

pub use analysis::{
    column_conservation, column_entropy, consensus, conservation_profile, entropy_profile,
    ConsensusOptions,
};
pub use guidetree::{
    distance_matrix, neighbor_joining, upgma, DistanceMatrix, GuideTree, TreeNode,
};
pub use profile::{align_profile_sequence, align_profiles, Profile, ProfileAlignment};
pub use progressive::{align, Msa};
pub use refine::{refine, RefineParams};
