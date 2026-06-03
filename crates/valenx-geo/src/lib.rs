//! # valenx-geo
//!
//! Canonical geometry types. CAD adapters produce [`Geometry`];
//! meshers consume it. Nobody downstream sees OpenCASCADE or
//! FreeCAD's internal BRep — only these types.
//!
//! Defined by [ARCHITECTURE.md § 4](../ARCHITECTURE.md).

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

pub mod bounding_box;
pub mod geometry;
pub mod polyline;
pub mod source;

pub use bounding_box::BoundingBox;
pub use geometry::{BRepHandle, Geometry};
pub use source::SourceFormat;
