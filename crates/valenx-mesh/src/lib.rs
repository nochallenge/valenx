//! # valenx-mesh
//!
//! Canonical mesh types. Meshers produce [`Mesh`]; solvers consume
//! it. Independent of any specific mesher's file format.
//!
//! Defined by [ARCHITECTURE.md § 4](../ARCHITECTURE.md).

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale
// Surface future `&str` byte-offset slicing in clippy review — this
// crate parses untrusted text mesh formats (OBJ/PLY/etc.), where
// non-char-boundary slices panic. WARN (not deny): most existing slices
// are safe ASCII; this only flags NEW ones.
#![warn(clippy::string_slice)]

pub mod adjacency;
pub mod boolean;
pub mod cut;
pub mod decimate;
pub mod element;
pub mod format;
pub mod mesh;
pub mod quality;
pub mod region;
pub mod remesh;
pub mod repair;
pub mod smooth;
pub mod stats;
pub mod stl_write;
pub mod transform;

pub use adjacency::{
    build_edge_adjacency, build_face_adjacency, edge_nodes_for, face_nodes_for, BoundaryEdge,
    BoundaryFace, EdgeAdjacency, ElementEdgeRef, ElementFaceRef, FaceAdjacency, InteriorEdge,
    InteriorFace,
};
pub use cut::{intersect_plane, slice, LineSegment};
pub use decimate::quadric_error_decimate;
pub use element::{ElementBlock, ElementType};
pub use mesh::Mesh;
pub use quality::{
    aspect_ratio, aspect_ratio_histogram, equiangle_skewness, min_orthogonality,
    report as quality_report, signed_size, skewness_histogram, AspectRatioHistogram, QualityReport,
    SkewnessHistogram, DEFAULT_AR_BUCKETS, DEFAULT_SKEW_BUCKETS,
};
pub use region::{BoundaryGroup, Region};
pub use remesh::{collapse_short_edges, flip_to_improve_valence, isotropic, split_long_edges};
pub use repair::{boundary_loops, fill_holes, is_manifold, self_intersections};
pub use smooth::{edge_length_stats, laplacian, taubin, vertex_neighbors};
pub use stats::MeshStats;
pub use stl_write::write_stl_binary;
pub use transform::{mirror, rotate_axis, scale_per_axis, scale_uniform, translate, Axis, Plane};
