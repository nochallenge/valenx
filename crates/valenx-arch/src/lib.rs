//! # valenx-arch
//!
//! Arch / BIM workbench — walls, slabs, columns, beams, windows,
//! doors, stairs, roofs, spaces, plus IFC4 export, schedule (BOM-like)
//! reporting, and a BCF (BIM Collaboration Format) issue tracker.
//!
//! Phase 15 of the FreeCAD-parity roadmap.
//!
//! # Architecture
//!
//! An [`ArchDocument`] owns a flat list of `(id, ArchEntity)` pairs
//! plus a `next_id` counter. Each [`ArchEntity`] variant carries its
//! own parametric description (start/end points, height, thickness,
//! profile, etc.) and knows how to tessellate itself into a
//! [`valenx_cad::Solid`] for viewport rendering.
//!
//! Each entity exposes:
//! - a `tessellate()` method returning a [`valenx_cad::Solid`] for
//!   single-entity preview, and
//! - a contribution to the document-level fused mesh produced by
//!   [`ArchDocument::tessellate_all`].
//!
//! ## IFC export
//!
//! [`ifc::writer::write_document`] emits a minimal-but-real ISO-10303-21
//! IFC4 file with the canonical `IfcProject` / `IfcSite` / `IfcBuilding`
//! / `IfcBuildingStorey` hierarchy plus one IFC entity per wall, slab,
//! column, beam, window, door, and space.
//!
//! v1 limitations are documented honestly below.
//! Notably, the writer covers a tiny subset of the IFC4 schema's
//! ~1500 entity types — enough to round-trip through validators and
//! to give downstream tools (Revit, ArchiCAD, BlenderBIM, IFC.js) a
//! readable file, but not enough for production-grade BIM hand-off.
//!
//! ## BCF stub
//!
//! [`bcf`] exposes an in-memory [`Bcf`] model and a directory-form
//! writer. A true BCF file is a ZIP envelope; v1 emits the XML files
//! into a directory and documents the path forward (Phase 15.5: pull
//! in a `zip` workspace dep and zip the directory).
//!
//! # Example
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_arch::{ArchDocument, ArchEntity, WallParams};
//!
//! let mut doc = ArchDocument::new("House");
//! let wall = WallParams {
//!     start: Vector3::new(0.0, 0.0, 0.0),
//!     end: Vector3::new(5.0, 0.0, 0.0),
//!     height: 2.7,
//!     thickness: 0.2,
//!     material: "Concrete".into(),
//! };
//! let _id = doc.add_entity(ArchEntity::Wall(wall));
//! assert_eq!(doc.count(), 1);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bcf;
pub mod document;
pub mod entity;
pub mod error;
pub mod ifc;
pub mod mep;
pub mod persist;
pub mod schedule;
pub mod structural;

// Per-entity modules.
pub mod beam;
pub mod column;
pub mod door;
pub mod opening;
pub mod roof;
pub mod slab;
pub mod space;
pub mod stair;
pub mod wall;
pub mod window;

pub use bcf::{Bcf, BcfIssue, BcfStatus, BcfViewpoint};
pub use beam::{BeamParams, BeamSection};
pub use column::{ColumnParams, ColumnSection};
pub use document::ArchDocument;
pub use door::{DoorParams, DoorStyle, Side};
pub use entity::{ArchEntity, ArchEntityKind};
pub use error::{ArchError, ErrorCategory};
pub use mep::{
    CableSegmentParams, ConduitSegmentParams, DuctShape, DuctSegmentParams, EquipmentKind,
    FlowDirection, MepEquipmentParams, PipeSegmentParams,
};
pub use persist::ArchFile;
pub use roof::{RoofParams, RoofType};
pub use schedule::{Schedule, ScheduleEntry};
pub use slab::SlabParams;
pub use space::SpaceParams;
pub use stair::StairParams;
pub use structural::{
    asd_load_combination, export_structural_model, lrfd_factored_load, StructuralElement,
    StructuralLoad, StructuralMaterial, StructuralMember, StructuralModel, StructuralModelOptions,
    StructuralNode, StructuralSection, StructuralSupport, SupportKind,
};
pub use wall::WallParams;
pub use window::{WindowParams, WindowStyle};
