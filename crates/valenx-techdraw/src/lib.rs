//! # valenx-techdraw
//!
//! TechDraw workbench — generate 2D engineering drawings from 3D
//! solids. Front / top / right / iso projection views, hidden-line
//! removal, section cuts, dimensions (linear / angular / radial /
//! diameter), bills of materials, sheet templates (A4–A0 landscape)
//! with a title block, and export to SVG / PDF / DXF.
//!
//! Phase 5 of the FreeCAD-parity roadmap.
//!
//! # Architecture
//!
//! A [`Drawing`] owns a [`Sheet`] (paper size + title-block fields), a
//! list of [`View`]s (each a projection of a source solid onto a 2D
//! plane), and a list of [`Dimension`]s overlaid on the views.
//!
//! Each [`View`] carries a [`ViewKind`] (Front / Top / iso / custom
//! camera), placement on the sheet, scale, and pre-computed visible /
//! hidden edge segments. Generation walks the source solid's BRep
//! edges (or mesh-triangle edges for mesh-backed solids), projects
//! them through the view's camera matrix, then runs an approximate
//! z-buffer hidden-line pass to classify each edge.
//!
//! Section views additionally cut the solid by a plane (delegating
//! to [`valenx_mesh::cut::intersect_plane_triangles`]) and hatch the
//! resulting cross-section with 45° parallel lines.
//!
//! Export targets are implemented in [`export`]:
//! - SVG via a hand-rolled writer (pure XML text — no transitive
//!   dependency churn).
//! - PDF: a minimal PDF 1.4 writer that emits one page per sheet with
//!   line segments and basic text.
//! - DXF: a minimal AutoCAD R12 DXF writer (HEADER + ENTITIES with
//!   LINE + TEXT records).
//!
//! Hand-rolling all three formats avoids pulling heavy unrelated
//! crates into the workspace (every export target is a few hundred
//! lines of text/binary serialization — well within Phase 5 scope).
//!
//! # Example
//!
//! ```no_run
//! use valenx_techdraw::{Drawing, Sheet, View, ViewKind};
//!
//! let mut d = Drawing::new(Sheet::a4_landscape("Bracket", "A. Engineer", "A"));
//! d.add_view(View::new(ViewKind::Front, 1.0, [80.0, 100.0]));
//! d.add_view(View::new(ViewKind::Top,   1.0, [80.0,  30.0]));
//! d.add_view(View::new(ViewKind::Isometric, 1.0, [220.0, 100.0]));
//! assert_eq!(d.views.len(), 3);
//! ```
//!
//! ## Persistence
//!
//! [`persist::TechDrawFile`] wraps a [`Drawing`] with a format version
//! and round-trips via RON (`from_ron` / `to_ron`, `read_from` /
//! `write_to`).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod balloon;
pub mod bom;
pub mod broken_view;
pub mod detail_view;
pub mod dim_chain;
pub mod dimension;
pub mod document;
pub mod error;
pub mod export;
pub mod gdt;
pub mod hatch_lib;
pub mod hlr;
pub mod leader;
pub mod parametric_view;
pub mod persist;
pub mod projection;
pub mod projection_group;
pub mod revision_block;
pub mod section;
pub mod sheet;
pub mod surface_finish;
pub mod view;
pub mod weld;

pub use balloon::{Balloon, BalloonStyle};
pub use bom::{Bom, BomItem};
pub use broken_view::{BreakAxis, BreakRegion, BreakStyle, BrokenEdges};
pub use detail_view::DetailView;
pub use dim_chain::{DimChain, DimChainKind};
pub use dimension::Dimension;
pub use document::Drawing;
pub use error::{ErrorCategory, TechDrawError};
pub use gdt::{Datum, DatumRef, GdtSymbol, GeometricCharacteristic, MaterialCondition};
pub use hatch_lib::HatchPattern;
pub use leader::{ArrowKind, Leader};
pub use parametric_view::ParametricView;
pub use persist::TechDrawFile;
pub use projection_group::{Projection, ProjectionGroup};
pub use revision_block::{RevisionBlock, RevisionEntry};
pub use sheet::{Sheet, SheetSize, SheetTemplate};
pub use surface_finish::{LayPattern, SurfaceFinish, SurfaceProcess};
pub use view::{View, ViewKind};
pub use weld::{WeldPosition, WeldSymbol, WeldType};
