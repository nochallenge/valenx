//! Hidden Markov models for sequence alignment and search
//! (HMMER-class core).
//!
//! Three model families, all with log-space dynamic programming:
//!
//! - [`pairhmm`] — a [`pairhmm::PairHmm`], the three-state (M/X/Y)
//!   model that generates two sequences jointly; Viterbi gives the
//!   most probable alignment, Forward the relatedness likelihood.
//! - [`profilehmm`] — a [`profilehmm::ProfileHmm`], the Plan7-style
//!   match/insert/delete profile built from an MSA; Viterbi / Forward
//!   score a query against the family.
//! - [`pssm`] — a [`pssm::Pssm`] position-specific scoring matrix /
//!   position-weight matrix with motif scanning.

pub mod pairhmm;
pub mod profilehmm;
pub mod pssm;

pub use pairhmm::PairHmm;
pub use profilehmm::{HmmNode, ProfileHmm};
pub use pssm::{MotifHit, Pssm};
