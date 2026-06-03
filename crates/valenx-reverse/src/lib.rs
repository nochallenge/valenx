//! # valenx-reverse
//!
//! Reverse Engineering workbench — turn raw scan data (point clouds)
//! into solid models. This is the FreeCAD `ReverseEngineering`
//! community workbench equivalent.
//!
//! Phase 26 of the FreeCAD-parity roadmap.
//!
//! # Pipeline
//!
//! 1. [`pointcloud::PointCloud`] — load points (`from_ply`,
//!    `from_xyz`) optionally with normals.
//! 2. [`pointcloud::estimate_normals`] — PCA-based per-point normal
//!    fit using k nearest neighbours. Adds normals when the loader
//!    didn't supply them.
//! 3. [`pointcloud::triangulate`] — k-NN local-Delaunay-style mesh
//!    reconstruction. v1 is "ball-pivot lite": for each point,
//!    connect it to its k nearest neighbours then keep only the
//!    triangles where the third vertex is mutual.
//! 4. [`reverse::cloud_to_brep`] — orchestrates 2 and 3 then calls
//!    `valenx_mesh_to_brep::brep_from_mesh` to produce a
//!    [`valenx_cad::Solid::Mesh`].
//!
//! ## v1 limitations
//!
//! - Triangulation is a deliberately simple k-NN mutual-pair
//!   reconstruction — good enough for cleanly-sampled industrial scans
//!   with low noise; not competitive with screened Poisson
//!   reconstruction. The follow-up (Phase 26.5) will swap in a real
//!   surface-reconstruction algorithm.
//! - PLY parser handles ASCII `ply / format ascii 1.0` files with
//!   `element vertex N` headers and `property float x/y/z` plus
//!   optional `nx/ny/nz` columns. Binary PLY and the colour /
//!   confidence channels are deferred.
//! - The final BRep is `valenx_cad::Solid::Mesh` — same v1 caveat as
//!   Phase 23 (mesh-backed; booleans + further BRep ops on it return
//!   `MeshBackedSolid`).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod persist;
pub mod pointcloud;
pub mod reverse;

pub use error::{ErrorCategory, ReverseError};
pub use persist::ReverseFile;
pub use pointcloud::{estimate_normals, from_ply, from_xyz, triangulate, PointCloud};
pub use reverse::cloud_to_brep;
