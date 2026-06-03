//! # valenx-viz
//!
//! The 3D viewport. Built on `wgpu`. Handles the camera model, picking
//! and selection, overlays (axes, ground, ViewCube), render styles
//! (flat / shaded / hidden-line / X-ray / section), and result
//! visualizations (contour / vector / iso-surface / streamline).
//!
//! Design specified in [DESIGN.md § 8](../DESIGN.md#8-viewport-design).
//!
//! ## What's shipped
//!
//! - [`stl`] — self-contained ASCII + binary STL loader producing a
//!   triangle soup with vertex positions and per-face normals. No
//!   wgpu dep, no GPU context required — callers can hand the output
//!   to a renderer or run unit tests on geometry properties.
//! - [`camera`] — the orbit-camera model (target + azimuth + elevation
//!   + distance) used by the ViewCube and the scroll-to-zoom gesture.
//!
//! The rendering passes (shaded / hidden-line / section) land with
//! the eframe shell in Phase 1.

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

pub mod camera;
pub mod projection;
pub mod scene;
pub mod stl;

pub use camera::{OrbitCamera, ProjectionMode, ViewDirection};
pub use projection::{project_point, project_triangle, ScreenPoint};
pub use scene::{
    gizmo_axis_screen_dirs, gizmo_view_for_face, grid_lod_params, intersect_ground_y0,
    nice_grid_spacing, ray_from_screen, snap_ground_point, GizmoFace, Ray,
};
pub use stl::{StlError, StlFormat, StlTriangle, TriangleMesh};
