//! # valenx-print-bed
//!
//! 3D-printer multi-part bed layout, optimal-orientation search,
//! tree-style support generation, and slicer-bundle export.
//!
//! Phase 51 of the FreeCAD-parity roadmap.  FreeCAD `Print Bed
//! Layout` community workbench equivalent.
//!
//! # Surface
//!
//! - [`printer::Printer`] + [`printer::Part`] + [`printer::BedType`]
//!   + [`printer::BedMaterial`].
//! - [`nest::auto_pack`] — first-fit decreasing 2D bin pack.
//! - [`orient::optimal`] — pick the best face-down orientation.
//! - [`support::generate`] — tree-style supports under overhangs.
//! - [`gcode::export_layout`] — STL-bundle export to a directory.
//! - [`panel::PrintBedPanelState`] — UI state envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod gcode;
pub mod nest;
pub mod orient;
pub mod panel;
pub mod printer;
pub mod support;

pub use error::{ErrorCategory, PrintBedError};
pub use gcode::{export_layout, SlicerSettings};
pub use nest::auto_pack;
pub use orient::{optimal, Criterion};
pub use panel::PrintBedPanelState;
pub use printer::{BedMaterial, BedType, Part, Printer};
pub use support::generate as support_generate;
