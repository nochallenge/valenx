//! Comparative / homology modelling — Modeller-class.
//!
//! **Roadmap features 1-6.** Comparative modelling is the oldest and
//! still most reliable structure-prediction approach *when a
//! homologous template exists*: a protein of known structure that is
//! evolutionarily related to the target. The pipeline is:
//!
//! 1. **Template search** ([`search`]) — score candidate templates
//!    by sequence alignment to the query; pick the best.
//! 2. **Target-template alignment** ([`align`]) — a careful pairwise
//!    alignment that defines which target residue maps to which
//!    template residue.
//! 3. **Backbone transfer** ([`transfer`]) — copy the template's
//!    backbone coordinates for the aligned ("equivalenced") regions.
//! 4. **Loop modelling** ([`loops`]) — build the insertions and the
//!    gap regions the template does not cover, by cyclic-coordinate-
//!    descent (CCD) loop closure.
//! 5. **Sidechain placement** ([`sidechains`]) — place the target's
//!    sidechains with the rotamer library.
//! 6. **Spatial-restraint assembly** ([`restraints`]) — Modeller's
//!    core idea: derive distance / dihedral restraints from the
//!    template and refine the model to satisfy them.
//!
//! **Honest accuracy note.** A comparative model is as good as its
//! template and alignment. With a close template (> 50 % identity)
//! the backbone is reliable; in the "twilight zone" (< 25 %) it is
//! not. This is genuinely useful classical modelling — it is *not*
//! AlphaFold, which learns a fold even with no template at all.

pub mod align;
pub mod loops;
pub mod restraints;
pub mod search;
pub mod sidechains;
pub mod transfer;

pub use align::{target_template_alignment, TargetTemplateAlignment};
pub use loops::{close_loop, model_loops, LoopClosure};
pub use restraints::{derive_restraints, satisfy_restraints, SpatialRestraint};
pub use search::{search_templates, TemplateCandidate, TemplateHit};
pub use sidechains::place_sidechains;
pub use transfer::transfer_backbone;
