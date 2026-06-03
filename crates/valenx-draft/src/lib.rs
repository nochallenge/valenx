//! # valenx-draft
//!
//! Draft workbench — 2D drawing entities (lines, polylines, arcs,
//! circles, rectangles, regular polygons, dimensions, text labels)
//! placed on a [`WorkingPlane`] in 3D space. Not parametric; meant
//! for construction lines, references, and annotations on top of the
//! Sketcher / Part Design workbenches.
//!
//! Phase 4 of the FreeCAD-parity roadmap.
//!
//! # Example
//!
//! ```
//! use valenx_draft::{DraftDocument, DraftEntity, WorkingPlane};
//!
//! let mut doc = DraftDocument::new(WorkingPlane::from_xy());
//! doc.add_entity(DraftEntity::Line { start: [0.0, 0.0], end: [10.0, 0.0] });
//! doc.add_entity(DraftEntity::Circle { center: [5.0, 5.0], radius: 2.0 });
//! assert_eq!(doc.entity_count(), 2);
//! ```
//!
//! ## Persistence
//!
//! Documents serialize to RON via [`persist::DraftFile::write_to`] /
//! [`persist::DraftFile::read_from`].
//!
//! ## Snapping
//!
//! [`snap`] exposes helpers for endpoint, midpoint, nearest-candidate,
//! and grid snapping. These are pure-data helpers — the UI layer uses
//! them to drive visual snap markers and to round cursor coordinates
//! before persisting a click.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod document;
pub mod entity;
pub mod error;
pub mod persist;
pub mod plane;
pub mod snap;

pub use document::DraftDocument;
pub use entity::DraftEntity;
pub use error::{DraftError, ErrorCategory};
pub use persist::DraftFile;
pub use plane::WorkingPlane;
