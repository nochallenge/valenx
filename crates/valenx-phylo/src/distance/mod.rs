//! Distance-based phylogenetic inference.
//!
//! The classic two-step pipeline: turn an alignment into a pairwise
//! [`DistMatrix`] of evolutionary distances, then cluster that matrix
//! into a tree.
//!
//! - [`matrix`] — distance estimation from a multiple alignment:
//!   p-distance and the Jukes-Cantor, Kimura-2-parameter and
//!   Tamura-Nei model corrections.
//! - [`cluster`] — the tree-building algorithms: UPGMA, WPGMA,
//!   neighbor-joining and BIONJ.
//!
//! UPGMA / WPGMA produce *ultrametric*, rooted trees (they assume a
//! molecular clock); NJ / BIONJ produce *additive*, unrooted trees and
//! do not.

pub mod cluster;
pub mod matrix;

pub use cluster::{bionj, neighbor_joining, upgma, wpgma};
pub use matrix::{distance_matrix, DistMatrix, DistanceModel};
