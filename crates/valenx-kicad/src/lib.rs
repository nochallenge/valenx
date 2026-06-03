//! # valenx-kicad
//!
//! KiCad PCB round-trip — parametric [`KicadBoard`] description,
//! `.kicad_pcb` S-expression reader (minimal subset: outline +
//! drills + components), board-to-solid tessellator, and STEP
//! assembly export hook (delegates to Phase 6 assembly machinery).
//!
//! Phase 42 of the FreeCAD-parity roadmap. KicadStepUp community
//! workbench analogue.
//!
//! # Surface
//!
//! - [`KicadBoard`] / [`Pad`] / [`Component`] / [`PadShape`] — data
//!   model.
//! - [`tessellate::pcb_to_solid`] — extrude outline + drill-hole
//!   visual proxy.
//! - [`parse::import_kicad_pcb`] / [`parse::from_str`] — `.kicad_pcb`
//!   S-expression reader.
//! - [`export::build_assembly`] — bridge to
//!   [`valenx_assembly::Assembly`] using user-supplied 3D models.
//! - [`export::export_step_with_components`] — STEP-assembly emit
//!   (returns `KicadError::NotImplemented` until truck-stepio
//!   gains assembly support).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod board;
pub mod error;
pub mod export;
pub mod parse;
pub mod tessellate;

pub use board::{Component, KicadBoard, Pad, PadShape};
pub use error::{ErrorCategory, KicadError};
pub use export::{build_assembly, export_step_with_components};
pub use parse::{from_str, import_kicad_pcb};
pub use tessellate::pcb_to_solid;
