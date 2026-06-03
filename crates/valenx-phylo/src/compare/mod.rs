//! Tree support, comparison and manipulation.
//!
//! - [`bootstrap`] — non-parametric bootstrap support: resample
//!   alignment columns, rebuild a tree per replicate, and map the
//!   replicate frequencies of each clade onto a reference tree.
//! - [`distance`] — topological distances between two trees: the
//!   Robinson-Foulds (symmetric bipartition) distance and the quartet
//!   distance.
//! - [`consensus`] — majority-rule and strict consensus of a tree set.
//! - [`manipulate`] — structural editing: reroot (outgroup / midpoint),
//!   ladderise, prune taxa and extract a subtree.

pub mod bootstrap;
pub mod consensus;
pub mod distance;
pub mod manipulate;

pub use bootstrap::{bootstrap_support, BootstrapResult};
pub use consensus::{consensus_tree, ConsensusKind};
pub use distance::{quartet_distance, robinson_foulds, RfResult};
pub use manipulate::{ladderize, midpoint_root, prune_taxa, reroot, subtree};
