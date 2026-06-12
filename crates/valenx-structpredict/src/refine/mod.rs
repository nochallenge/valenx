//! Model refinement & quality assessment.
//!
//! **Roadmap features 13-17.** A freshly-built model — whether from
//! homology or ab-initio — needs refinement: its sidechains may
//! clash, its backbone may be strained, and its quality needs an
//! honest score. This module group does that:
//!
//! - **Rotamer repacking** ([`repack`]) — resolve sidechain clashes
//!   by a combinatorial rotamer search (dead-end elimination plus
//!   simulated annealing — the SCWRL-class algorithm).
//! - **Energy minimisation** ([`minimize`]) — relax the model by
//!   gradient minimisation, driven by the [`valenx_md`] force-field
//!   minimisers against a structure-energy potential.
//! - **Ramachandran refinement** ([`ramachandran`]) — pull strained
//!   backbone φ/ψ angles into the allowed Ramachandran regions.
//! - **Quality assessment** ([`quality`]) — a classical model-quality
//!   score: clash score, Ramachandran outliers, packing.
//! - **Superposition** ([`superpose`]) — superpose a model on a
//!   reference and report RMSD and GDT (reuses
//!   [`valenx_biostruct`]'s Kabsch).

pub mod mcrefine;
pub mod minimize;
pub mod quality;
pub mod ramachandran;
pub mod repack;
pub mod superpose;

pub use mcrefine::{mc_refine, McRefineOptions, McRefineResult};
pub use minimize::{relax_model, RelaxResult};
pub use quality::{assess_quality, QualityReport};
pub use ramachandran::{refine_ramachandran, RamachandranRefinement};
pub use repack::{repack_sidechains, RepackResult};
pub use superpose::{ca_rmsd_superposed, gdt_ts, model_vs_reference, SuperposeReport};
