//! Physics-based protein design (non-ML).
//!
//! **Roadmap features 18-22.** Protein *design* is the inverse of
//! prediction: given a backbone, find a sequence that folds into it.
//! The neural approach (ProteinMPNN, RFdiffusion) is adapter-only
//! here; this module implements the **classical, physics-based
//! design** that Rosetta `fixbb` made standard — search amino-acid
//! sequence and rotamer space to minimise an energy on a fixed
//! backbone.
//!
//! - **Fixed-backbone design** ([`fixbb`]) — the Rosetta `fixbb`
//!   problem: choose the sequence of lowest energy for a given
//!   backbone.
//! - **Combinatorial design search** ([`search`]) — simulated
//!   annealing over the joint (residue → amino-acid, rotamer) space.
//! - **Design scoring** ([`score`]) — the statistical / physics
//!   potential the design search minimises.
//! - **Motif grafting** ([`graft`]) — graft a functional loop / motif
//!   onto a scaffold.
//! - **Interface design** ([`interface`]) — design a protein-protein
//!   binding interface.
//!
//! **Honest note.** Classical fixed-backbone design with a knowledge-
//! based score is a real, historically productive method — many
//! designed proteins predate the deep-learning era. It is *not*
//! ProteinMPNN: it has no learnt sequence model, and its success
//! rate on a hard backbone is lower. Use the ProteinMPNN /
//! RFdiffusion adapters when you need the network.

pub mod fixbb;
pub mod graft;
pub mod interface;
pub mod score;
pub mod search;

pub use fixbb::{design_fixed_backbone, FixbbResult};
pub use graft::{graft_motif, GraftResult};
pub use interface::{design_interface, InterfaceDesign};
pub use score::{design_score, DesignScore, DesignScoreWeights};
pub use search::{combinatorial_design, DesignSearchResult, ResiduePalette};
