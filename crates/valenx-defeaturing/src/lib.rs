//! # valenx-defeaturing
//!
//! Auto-removal of small features (fillets, holes, engraved text,
//! slivers) from imported models. v1: tessellation-based — walks
//! the triangulated mesh, classifies candidate triangles against
//! per-defeature size thresholds, and emits a cleaned mesh.
//!
//! Phase 43 of the FreeCAD-parity roadmap. FreeCAD `Defeaturing`
//! community workbench analogue.
//!
//! # Surface
//!
//! - [`Defeature`] — 4 variants: FilletRemove, HoleRemove,
//!   TextRemove, SliverRemove. Each carries its own size threshold.
//! - [`apply::apply`] — apply a sequence of defeatures and return
//!   the cleaned [`valenx_cad::Solid::Mesh`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod apply;
pub mod defeature;
pub mod error;

pub use apply::{apply, DEFAULT_TOLERANCE_MM};
pub use defeature::Defeature;
pub use error::{DefeatureError, ErrorCategory};
