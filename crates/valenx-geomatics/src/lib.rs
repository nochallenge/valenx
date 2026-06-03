//! # valenx-geomatics
//!
//! Geomatics / GIS workbench — surveying + WGS84 / UTM coordinate
//! conversion + Digital Elevation Models. The FreeCAD `Geomatics`
//! community workbench equivalent.
//!
//! Phase 35 of the FreeCAD-parity roadmap.
//!
//! # Surface
//!
//! - [`LatLon`] + [`Utm`] — coordinate types.
//! - [`coord::wgs84_to_utm`] / [`coord::utm_to_wgs84`] — Kruger /
//!   Karney series approximation suitable for typical engineering
//!   work (sub-metre accuracy within the 6° UTM strip).
//! - [`Dem`] — regular grid Digital Elevation Model.
//! - [`dem::from_xyz_ascii`] — parser for `x y z` lines (whitespace
//!   separated). The grid spacing and dimensions are inferred from
//!   the data — non-regular inputs return `GeomaticsError`.
//! - [`dem::to_mesh`] — build a Tri3 surface mesh from the grid.
//! - [`dem::sample`] — bilinear interpolation at `(x, y)`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod coord;
pub mod dem;
pub mod error;

pub use coord::{LatLon, Utm, Hemisphere};
pub use dem::Dem;
pub use error::{ErrorCategory, GeomaticsError};
