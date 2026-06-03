//! Maximum-parsimony phylogenetics.
//!
//! Parsimony scores a tree by the minimum number of character-state
//! changes it implies for an alignment; the most-parsimonious tree is
//! the one needing the fewest.
//!
//! - [`fitch`] — the **small-parsimony** problem (Fitch 1971): given a
//!   *fixed* tree, find the minimum change count for one character and
//!   the ancestral state *sets* at every internal node.
//! - [`sankoff`] — the **weighted** generalisation (Sankoff 1975):
//!   different state-to-state changes carry different costs, supplied
//!   as a cost matrix.
//! - [`search`] — the **large-parsimony** problem: search tree space
//!   (NNI + SPR rearrangements, hill-climbing) for a low-score
//!   topology. Large-parsimony is NP-hard, so this is a heuristic.
//!
//! Characters are `u8` state indices; the alignment is a slice of
//! equal-length rows (gaps treated as a missing / wildcard state).

pub mod fitch;
pub mod sankoff;
pub mod search;

pub use fitch::{fitch_parsimony, FitchResult};
pub use sankoff::{sankoff_parsimony, CostMatrix, SankoffResult};
pub use search::{parsimony_search, ParsimonySearch, SearchReport};
