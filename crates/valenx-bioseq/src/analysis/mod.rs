//! Sequence analysis: composition (GC / skew / entropy), k-mer
//! statistics, melting temperature, molecular weight, ProtParam-class
//! protein properties, and primer-grade DNA thermodynamics (ΔG, ΔH,
//! ΔS, hairpin / dimer scoring).

pub mod composition;
pub mod kmer;
pub mod protparam;
pub mod thermo;
pub mod tm;
pub mod weight;
