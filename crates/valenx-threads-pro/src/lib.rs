//! # valenx-threads-pro
//!
//! Extended thread tables + helical thread geometry.  Phase 13
//! attached thread metadata to Hole / Bolt / Nut features but did NOT
//! sweep a real 3D thread profile; this crate fills that gap and adds
//! 9 more thread families on top of the Phase 13 four.
//!
//! Phase 48 of the FreeCAD-parity roadmap.  FreeCAD `Threaded
//! Profiles` community workbench equivalent.
//!
//! # Surface
//!
//! - [`standard::ThreadStandardPro`] (15 families) + [`standard::ProfileShape`]
//!   (V / Acme / Trapezoidal / Buttress).
//! - [`spec::ThreadSpecPro`] — one row in any table.
//! - [`tables`] — full ISO metric / metric fine / metric extra-fine /
//!   UNC / UNF / UNEF / BSPP / BSPT / NPT / NPTF / NPS / Acme /
//!   ISO Trapezoidal / Whitworth-BSW tables.
//! - [`thread::helix_polyline`] — sample a true helix.
//! - [`thread::profile_solid`] — sweep the profile to a mesh-backed
//!   Solid for preview.
//!
//! Phase 13's [`valenx_feature_tree::threads::ThreadStandard`] is
//! re-exported so existing Hole / Bolt code can still call into the
//! v1 four-family tags without an extra import.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod spec;
pub mod standard;
pub mod tables;
pub mod thread;

pub use error::{ErrorCategory, ThreadsProError};
pub use spec::ThreadSpecPro;
pub use standard::{ProfileShape, ThreadStandardPro};
pub use tables::{
    acme_table, bspp_table, bspt_table, metric_extra_fine_table, metric_fine_table,
    metric_full_table, nps_table, npt_table, nptf_table, trapezoidal_iso_table, unc_table,
    unef_table, unf_table, whitworth_bsw_table,
};
pub use thread::{helix_polyline, profile_solid};

// Re-export the Phase 13 family tag for downstream callers that still
// work with Bolt / Nut / Hole's ThreadSpec.
pub use valenx_feature_tree::threads::ThreadStandard;
