//! # valenx-meshpart
//!
//! MeshPart workbench — finer-grained mesh ↔ BRep round-trip
//! utilities than Phase 23. The FreeCAD `MeshPart` workbench
//! equivalent.
//!
//! Phase 32 of the FreeCAD-parity roadmap.
//!
//! # Surface
//!
//! - [`brep_to_polyhedron`] — explicit polyhedral tessellation with
//!   adjustable chord-error tolerance (thin shim over
//!   `valenx_cad::solid_to_mesh` that returns the canonical
//!   [`valenx_mesh::Mesh`]).
//! - [`polyhedron_to_brep`] — naïve mesh → BRep promoter: each
//!   triangle becomes one planar face; the result is wrapped as a
//!   [`valenx_cad::Solid::Mesh`] (true BRep sewing of fitted faces is
//!   still gated on truck-modeling APIs — Phase 32.5).
//! - [`segment_by_normal`] — group triangles by similar face normal
//!   (angle threshold in degrees).
//! - [`extract_boundary_loop`] — boundary polyline of a connected
//!   triangle group.
//! - [`flatten_boundary`] — project a 3D loop onto a 2D plane for
//!   sketch generation.
//! - [`triangulate_polygon`] — ear-clipping fill for a 2D polygon.
//! - [`merge_meshes`] — concatenate vertex arrays + connectivity with
//!   the right offsets.
//! - [`split_mesh_by_planes`] — cut along an axis-aligned plane list,
//!   return the regions.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod boundary;
pub mod convert;
pub mod error;
pub mod merge;
pub mod segment;
pub mod split;
pub mod triangulate;

pub use boundary::{extract_boundary_loop, flatten_boundary};
pub use convert::{brep_to_polyhedron, polyhedron_to_brep};
pub use error::{ErrorCategory, MeshPartError};
pub use merge::merge_meshes;
pub use segment::{segment_by_normal, TriangleGroup};
pub use split::split_mesh_by_planes;
pub use triangulate::triangulate_polygon;

use valenx_mesh::element::ElementBlock;

/// Validate that every Tri3 connectivity entry indexes an existing node.
///
/// `valenx_mesh::Mesh` exposes public `nodes`/`element_blocks` and enforces no
/// in-bounds invariant, so a corrupt or hand-built mesh can carry a
/// connectivity index past `nodes`. The mesh-part algorithms index
/// `mesh.nodes[conn[..]]` directly, so they call this first to fail with a
/// typed error rather than panic out of bounds.
pub(crate) fn check_connectivity(
    block: &ElementBlock,
    n_nodes: usize,
) -> Result<(), MeshPartError> {
    for &idx in &block.connectivity {
        if idx as usize >= n_nodes {
            return Err(MeshPartError::BadParameter {
                name: "connectivity",
                reason: format!("node index {idx} >= node count {n_nodes}"),
            });
        }
    }
    Ok(())
}
