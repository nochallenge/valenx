//! Triangle-mesh primitive — vertices + triangle index list.
//!
//! libigl operates on `V` (V x 3 matrix of vertices) + `F` (F x 3
//! matrix of triangle indices). We mirror that with a simple owned
//! struct.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Triangle mesh.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TriMesh {
    /// Vertices.
    pub vertices: Vec<Vector3<f64>>,
    /// Triangles as 3 vertex indices each.
    pub triangles: Vec<[usize; 3]>,
}

impl TriMesh {
    /// Empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Vertex count.
    pub fn n_verts(&self) -> usize {
        self.vertices.len()
    }

    /// Triangle count.
    pub fn n_tris(&self) -> usize {
        self.triangles.len()
    }

    /// Centroid of all vertices.
    pub fn centroid(&self) -> Vector3<f64> {
        if self.vertices.is_empty() {
            return Vector3::zeros();
        }
        let s: Vector3<f64> = self.vertices.iter().sum();
        s / self.vertices.len() as f64
    }

    /// Adjacency — set of neighbour vertex ids per vertex. Used by
    /// the discrete Laplacian + heat-geodesics.
    pub fn vertex_one_ring(&self) -> Vec<std::collections::BTreeSet<usize>> {
        let mut out = vec![std::collections::BTreeSet::<usize>::new(); self.vertices.len()];
        for tri in &self.triangles {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                out[a].insert(b);
                out[b].insert(a);
            }
        }
        out
    }
}
