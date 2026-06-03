//! Format-specific readers + writers for the canonical types.
//! Each submodule covers one file format. Pull only the modules
//! your adapter needs to keep build times tight.

pub mod dcd;
pub mod fasta;
pub mod fastq;
pub mod mmcif;
pub mod pdb;
pub mod pdbqt;
pub mod sam;
pub mod vcf;
