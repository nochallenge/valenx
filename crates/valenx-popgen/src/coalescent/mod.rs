//! Backward (coalescent) simulation.
//!
//! This module group is the `msprime` / `ms` / `discoal` side of
//! `valenx-popgen`: instead of running a population forward, it traces
//! the genealogy of a sample *backward* in time.
//!
//! - [`kingman`] — the [`coalescent()`] simulator (Kingman 1982),
//!   piecewise [`PopHistory`] support and the
//!   [`structured_coalescent`] for multi-deme samples. Output is a
//!   [`valenx_phylo::Tree`].
//! - [`arg`] — the coalescent *with recombination*, an ancestral
//!   recombination graph, recorded directly into a [`TreeSequence`].
//! - [`tree_sequence`] — the succinct `tskit`-class [`TreeSequence`]
//!   (node / edge / site / mutation tables) with local-tree
//!   extraction.
//! - [`overlay`] — dropping infinite-sites mutations onto a genealogy
//!   or tree sequence to obtain a [`crate::infer::GenotypeMatrix`].

pub mod arg;
pub mod kingman;
pub mod overlay;
pub mod tree_sequence;

pub use arg::{simulate_arg, ArgParams, RecombinationMap};
pub use kingman::{coalescent, structured_coalescent, PopHistory};
pub use overlay::{overlay_mutations, overlay_on_tree};
pub use tree_sequence::{Edge, TreeSequence, TsMutation, TsNode, TsSite};
