//! # valenx-fields
//!
//! Canonical results types for Valenx — `Units`, `TimeKey`, `Field`,
//! `ScalarRecord`, `Provenance`, `Results`. Specified by
//! [RFC 0004](../rfcs/0004-results-and-fields.md).
//!
//! Every adapter's `collect()` returns a [`Results`]; every downstream
//! consumer (viewer, plot editor, table view, export pipeline, report
//! generator) reads from it. This is the single shared vocabulary for
//! simulation output.

#![forbid(unsafe_code)]
// `missing_docs` will be re-enabled to warn as the pre-alpha code
// matures and every public item gets a comment. Gating it now adds
// ~hundreds of warnings on purely mechanical doc stubs (SI unit
// constants, field accessors) without improving the code.
#![allow(missing_docs)]

pub mod artifact;
pub mod catalog;
pub mod colormap;
pub mod field;
pub mod frd;
pub mod integrals;
pub mod interp;
pub mod provenance;
pub mod pvd;
pub mod results;
pub mod scalar;
pub mod stress;
pub mod time;
pub mod units;
pub mod vtk_dispatch;
pub mod vtk_legacy;
pub mod vtu;

// Re-exports for convenience.
pub use artifact::{Artifact, ArtifactKind};
pub use catalog::{FieldCatalog, ScalarCatalog};
pub use field::{Field, FieldKind, Location, RegionRef};
pub use provenance::{Provenance, ProvenanceRef, Sha256Hex};
pub use results::{ResultMeta, Results};
pub use scalar::{ScalarRecord, ScalarSource};
pub use time::TimeKey;
pub use units::Units;
