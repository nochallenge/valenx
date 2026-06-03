//! Tree file I/O — readers and writers for the standard phylogenetic
//! exchange formats.
//!
//! - [`newick`] — the ubiquitous parenthesised
//!   `((A:0.1,B:0.2):0.3,C:0.5);` format (read + write).
//! - [`nexus`] — the NEXUS block format with `TAXA`, `TREES` and a
//!   minimal `DATA` block (read + write).
//! - [`xml`] — PhyloXML and NeXML writers.
//!
//! Newick is the lingua franca; NEXUS and the XML formats are layered
//! on top of the same Newick string builder where they embed a tree.

pub mod newick;
pub mod nexus;
pub mod xml;

pub use newick::{read_newick, write_newick};
pub use nexus::{read_nexus, write_nexus, NexusFile};
pub use xml::{write_nexml, write_phyloxml};
