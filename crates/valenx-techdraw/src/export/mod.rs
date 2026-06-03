//! Export drawings to SVG, PDF, and DXF.
//!
//! Hand-rolled writers — Phase 5 avoids pulling heavy external crates
//! (`svg`, `printpdf`, `dxf`) because each one drags in a tangle of
//! transitive dependencies that bumps workspace build time and
//! conflicts with our pinned versions. Every format we ship is
//! within a few hundred lines of plain text / binary serialization,
//! so we keep them in-tree where they're easy to maintain.
//!
//! Each submodule exposes a `write(drawing, path)` entry point:
//! - [`svg::write`] — SVG 1.1 plain text.
//! - [`pdf::write`] — minimal PDF 1.4, one page per sheet.
//! - [`dxf::write`] — AutoCAD R12 ASCII DXF.
//!
//! All three accept the same [`crate::Drawing`] and produce a file
//! readable by the standard tools in their respective ecosystems
//! (browsers / Inkscape for SVG, any PDF viewer, AutoCAD / FreeCAD /
//! LibreCAD for DXF).

pub mod dxf;
pub mod pdf;
pub mod svg;
