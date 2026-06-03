//! Structure comparison and analysis.
//!
//! Tools for comparing structures and folding harder structure
//! classes:
//!
//! - [`distance`] — base-pair distance and the Zhang-Shasha
//!   tree-edit distance between two structures.
//! - [`align`] — RNAforester-class structure (tree) alignment.
//! - [`consensus`] — RNAalifold-class covariation-aware consensus
//!   folding of an alignment.
//! - [`pseudoknot`] — folding of restricted (H-type) pseudoknots.
//! - [`pknots_rg`] — pknotsRG-class pseudoknot folding (H-type and
//!   kissing-hairpin).

pub mod align;
pub mod consensus;
pub mod distance;
pub mod pknots_rg;
pub mod pseudoknot;
