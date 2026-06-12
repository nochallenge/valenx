//! # valenx-occt-exchange
//!
//! Phases 101-130 — OpenCASCADE (OCCT) **data
//! exchange** feature parity for Valenx.
//!
//! OCCT's `DataExchange` package implements ~30 separate importers /
//! exporters across the STEP family (AP203 / AP214 / AP242),
//! IGES 5.3, the proprietary BREP formats (ACIS .sat, Parasolid X_T,
//! Siemens JT), and the mesh-interchange family (OBJ / PLY / STL /
//! glTF / VRML / X3D / COLLADA). This crate provides a Rust-native
//! API surface for 30 of those, one module per format-or-direction
//! so each can advance from "scaffold" → "v1 implementation" →
//! "production parity" independently.
//!
//! ## v1 strategy
//!
//! Each of the 30 functions either:
//!
//! 1. Maps cleanly onto an existing `valenx-step-iges` or
//!    `valenx-mesh` capability, in which case it implements honestly
//!    (parameter validation + delegation + typed errors); or
//! 2. Returns [`OcctExchangeError::NotYetImplemented`] with rustdoc
//!    describing what the real OCCT API does and which Phase `N.5`
//!    follow-up will deliver the deep implementation.
//!
//! Both kinds carry a public function signature so downstream crates
//! (the toolbox UI, the import/export dispatcher, integration tests)
//! can be written against the final API today and have the stubs
//! fill in later without churn.
//!
//! ## Feature catalogue (Phases 101-130)
//!
//! ### STEP family (Phases 101-108)
//! - [`step_ap203_writer()`], [`step_ap203_reader()`],
//!   [`step_ap214_writer()`], [`step_ap214_reader()`],
//!   [`step_ap242_full_writer()`], [`step_ap242_full_reader()`],
//!   [`step_ap203_assembly_writer()`], [`step_color_attributes_writer()`].
//!
//! ### IGES family (Phases 109-113)
//! - [`iges_5_3_writer()`], [`iges_5_3_reader()`],
//!   [`iges_trimmed_surface_writer()`],
//!   [`iges_trimmed_surface_reader()`], [`iges_color_attributes()`].
//!
//! ### Proprietary BREP (Phases 114-119)
//! - [`acis_sat_writer()`], [`acis_sat_reader()`],
//!   [`parasolid_xt_writer()`], [`parasolid_xt_reader()`],
//!   [`jt_writer()`], [`jt_reader()`].
//!
//! ### Mesh formats (Phases 120-127)
//! - [`obj_writer_extended()`], [`obj_reader_extended()`],
//!   [`ply_writer_extended()`], [`ply_reader_extended()`],
//!   [`stl_writer_extended()`], [`stl_reader_extended()`],
//!   [`gltf2_writer()`], [`gltf2_reader()`].
//!
//! ### Misc exchange (Phases 128-130)
//! - [`vrml_writer()`], [`x3d_writer()`], [`collada_writer()`].
//!
//! ## Error model
//!
//! All public APIs return [`Result<_, OcctExchangeError>`]. See
//! [`OcctExchangeError::code`] / [`OcctExchangeError::category`] for
//! the stable taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;

// STEP family (Phases 101-108)
pub mod step_ap203_assembly_writer;
pub mod step_ap203_reader;
pub mod step_ap203_writer;
pub mod step_ap214_reader;
pub mod step_ap214_writer;
pub mod step_ap242_full_reader;
pub mod step_ap242_full_writer;
pub mod step_color_attributes_writer;

// IGES family (Phases 109-113)
pub mod iges_5_3_reader;
pub mod iges_5_3_writer;
pub mod iges_color_attributes;
pub mod iges_trimmed_surface_reader;
pub mod iges_trimmed_surface_writer;

// Proprietary BREP (Phases 114-119)
pub mod acis_sat_reader;
pub mod acis_sat_writer;
pub mod jt_reader;
pub mod jt_writer;
pub mod parasolid_xt_reader;
pub mod parasolid_xt_writer;

// Mesh formats (Phases 120-127)
pub mod gltf2_reader;
pub mod gltf2_writer;
pub mod obj_reader_extended;
pub mod obj_writer_extended;
pub mod ply_reader_extended;
pub mod ply_writer_extended;
pub mod stl_reader_extended;
pub mod stl_writer_extended;

// Misc exchange (Phases 128-130)
pub mod collada_writer;
pub mod vrml_writer;
pub mod x3d_writer;

pub use error::{ErrorCategory, OcctExchangeError};

// Re-export the entry points so callers can `use valenx_occt_exchange::step_ap203_writer;`
// instead of the full module-path mouthful.
pub use acis_sat_reader::acis_sat_reader;
pub use acis_sat_writer::acis_sat_writer;
pub use collada_writer::collada_writer;
pub use gltf2_reader::gltf2_reader;
pub use gltf2_writer::gltf2_writer;
pub use iges_5_3_reader::iges_5_3_reader;
pub use iges_5_3_writer::iges_5_3_writer;
pub use iges_color_attributes::iges_color_attributes;
pub use iges_trimmed_surface_reader::iges_trimmed_surface_reader;
pub use iges_trimmed_surface_writer::iges_trimmed_surface_writer;
pub use jt_reader::{jt_reader, read_jt_model, JtModel, JtNode, JtTocEntry};
pub use jt_writer::jt_writer;
pub use obj_reader_extended::obj_reader_extended;
pub use obj_writer_extended::obj_writer_extended;
pub use parasolid_xt_reader::parasolid_xt_reader;
pub use parasolid_xt_writer::parasolid_xt_writer;
pub use ply_reader_extended::ply_reader_extended;
pub use ply_writer_extended::ply_writer_extended;
pub use step_ap203_assembly_writer::step_ap203_assembly_writer;
pub use step_ap203_reader::step_ap203_reader;
pub use step_ap203_writer::step_ap203_writer;
pub use step_ap214_reader::step_ap214_reader;
pub use step_ap214_writer::step_ap214_writer;
pub use step_ap242_full_reader::step_ap242_full_reader;
pub use step_ap242_full_writer::step_ap242_full_writer;
pub use step_color_attributes_writer::step_color_attributes_writer;
pub use stl_reader_extended::stl_reader_extended;
pub use stl_writer_extended::stl_writer_extended;
pub use vrml_writer::vrml_writer;
pub use x3d_writer::x3d_writer;
