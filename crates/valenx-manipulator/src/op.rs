//! Direct-edit operations enum.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Push/pull operations on a mesh-backed solid.
///
/// v1 indexes are into the tessellated mesh: `face_idx` is a
/// triangle-block index, `edge_idx` and `vertex_idx` are unique
/// edge / vertex ids built from the triangle list.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ManipulateOp {
    /// Move every vertex of the selected triangle face by `delta`.
    MoveFace {
        /// Triangle face index in the mesh.
        face_idx: usize,
        /// World-space translation.
        delta: Vector3<f64>,
    },
    /// Rotate the selected triangle face around an arbitrary axis
    /// through the face centroid.
    RotateFace {
        /// Triangle face index.
        face_idx: usize,
        /// Rotation axis (will be normalised).
        axis: Vector3<f64>,
        /// Angle in degrees.
        angle_deg: f64,
    },
    /// Move both vertices of the selected unique mesh edge by
    /// `delta`.
    MoveEdge {
        /// Unique edge id.
        edge_idx: usize,
        /// World-space translation.
        delta: Vector3<f64>,
    },
    /// Move a single vertex by `delta`.
    MoveVertex {
        /// Vertex id (node index in the mesh).
        vertex_idx: usize,
        /// World-space translation.
        delta: Vector3<f64>,
    },
    /// Extrude the selected triangle face along its normal by
    /// `distance`. v1 duplicates the face vertices and adds side
    /// walls connecting old to new vertices.
    ExtrudeFace {
        /// Triangle face index.
        face_idx: usize,
        /// Distance along the outward normal.
        distance: f64,
    },
    /// Offset the selected triangle face along its normal by
    /// `distance` (no vertex duplication).
    OffsetFace {
        /// Triangle face index.
        face_idx: usize,
        /// Distance along the outward normal.
        distance: f64,
    },
}

impl ManipulateOp {
    /// Short kebab-case identifier (for error / log messages).
    pub fn kind(&self) -> &'static str {
        match self {
            Self::MoveFace { .. } => "move-face",
            Self::RotateFace { .. } => "rotate-face",
            Self::MoveEdge { .. } => "move-edge",
            Self::MoveVertex { .. } => "move-vertex",
            Self::ExtrudeFace { .. } => "extrude-face",
            Self::OffsetFace { .. } => "offset-face",
        }
    }
}
