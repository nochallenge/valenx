//! External-geometry references — sketch entities that project
//! geometry from another part of the feature tree into the sketcher.
//!
//! Phase 12D. An external edge / vertex / face from another solid is
//! projected onto the sketch plane and treated as a *fixed* primitive
//! (the solver freezes its variables; see
//! [`crate::sketch::Sketch::is_var_frozen`]).
//!
//! The actual projection happens in
//! [`crate::sketch::Sketch::resolve_externals`], which walks the
//! provided feature tree, finds each source entity, and rewrites the
//! frozen sketch primitives' variables to match the projected
//! coordinates.

use serde::{Deserialize, Serialize};

/// Opaque feature id — index into a caller-provided feature tree.
/// Mirrors the Phase 2 feature-tree shape; we don't depend on a
/// concrete tree type so this crate stays self-contained.
pub type FeatureId = u64;

/// Reference to a geometric element of another feature.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum ExternalGeomRef {
    /// Edge of a source feature's geometry.
    Edge {
        /// Identifier of the source feature in the parent tree.
        source_feature: FeatureId,
        /// Index of the edge within that feature's geometry.
        edge_index: usize,
    },
    /// Vertex of a source feature.
    Vertex {
        /// Identifier of the source feature in the parent tree.
        source_feature: FeatureId,
        /// Index of the vertex within that feature's geometry.
        vertex_index: usize,
    },
    /// Face of a source feature — projected as a polyline of edges.
    Face {
        /// Identifier of the source feature in the parent tree.
        source_feature: FeatureId,
        /// Index of the face within that feature's geometry.
        face_index: usize,
    },
}

/// Trait the host provides so [`crate::sketch::Sketch::resolve_externals`]
/// can look up source geometry without binding to a concrete tree.
///
/// Each lookup returns endpoints in world coordinates; the sketch is
/// assumed to be on the world XY plane (Phase 1's working assumption,
/// reused here) so projection is `z=0` after subtracting whatever
/// plane offset the host has.
pub trait FeatureTreeLookup {
    /// Return the world-space endpoints of an edge.
    fn edge_endpoints(&self, ext: &ExternalGeomRef) -> Option<((f64, f64), (f64, f64))>;
    /// Return the world-space coordinates of a vertex.
    fn vertex_xy(&self, ext: &ExternalGeomRef) -> Option<(f64, f64)>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn external_geom_ref_serializes_round_trip() {
        let e = ExternalGeomRef::Edge {
            source_feature: 7,
            edge_index: 2,
        };
        let ron = ron::ser::to_string(&e).unwrap();
        let back: ExternalGeomRef = ron::from_str(&ron).unwrap();
        assert_eq!(e, back);
    }
}
