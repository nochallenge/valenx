//! The opaque `Geometry` type — a handle to a BRep plus enough
//! metadata for the viewer and mesher to decide what to do with it.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::bounding_box::BoundingBox;
use crate::source::SourceFormat;

/// Opaque reference to an in-memory BRep. The actual BRep lives
/// behind an adapter (OpenCASCADE, fornjot, truck, …); callers move
/// the handle around but don't inspect its internals.
///
/// `backend` is a short owned string like `"opencascade"`, `"fornjot"`,
/// `"truck"`. Owned (not `&'static str`) so the type serializes cleanly.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BRepHandle {
    /// Backend the handle belongs to.
    pub backend: String,
    /// Stable opaque ID within the backend (e.g. a hash of the BRep
    /// contents or a handle into a backend-owned table).
    pub id: String,
}

/// Canonical geometry exchanged between CAD adapters and meshers.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Geometry {
    /// Project-local stable ID (used by `project.toml` entries).
    pub id: String,
    /// Optional on-disk source, relative to the project root.
    pub source_path: Option<PathBuf>,
    /// Format of the original source file, if any.
    pub source_format: Option<SourceFormat>,
    /// Axis-aligned bounding box, in model units (metres by
    /// default).
    pub bounds: BoundingBox,
    /// Handle to the BRep — absent for pure STL / OBJ imports that
    /// have no parametric representation.
    pub brep: Option<BRepHandle>,
    /// Number of top-level solids / shells in this geometry.
    pub solid_count: u32,
    /// Number of faces, if known.
    pub face_count: Option<u32>,
}

impl Geometry {
    /// Minimal constructor for when only a bounding box is known
    /// (e.g. a bare STL import).
    pub fn stl(id: impl Into<String>, bounds: BoundingBox) -> Self {
        Self {
            id: id.into(),
            source_path: None,
            source_format: Some(SourceFormat::Stl),
            bounds,
            brep: None,
            solid_count: 1,
            face_count: None,
        }
    }
}
