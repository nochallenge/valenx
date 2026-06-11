//! Polygon mesh primitive used by the subdivision schemes.
//!
//! Catmull-Clark needs general polygons (it produces quads); Loop
//! needs triangles. We model both with a single [`SubdivMesh`] —
//! faces are `Vec<usize>` so a triangle is length 3, a quad length
//! 4, etc.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Polygon mesh — vertices + variable-arity faces.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SubdivMesh {
    /// Vertex positions.
    pub vertices: Vec<Vector3<f64>>,
    /// Faces, each a vertex-index loop in CCW order.
    pub faces: Vec<Vec<usize>>,
}

impl SubdivMesh {
    /// Empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Vertex count.
    pub fn n_verts(&self) -> usize {
        self.vertices.len()
    }

    /// Face count.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Check that every face references only in-range vertex indices.
    ///
    /// `SubdivMesh` exposes public fields, so a caller (or a corrupt import)
    /// can build a face indexing past `vertices`; the subdivision schemes would
    /// then index out of bounds. Returns
    /// [`crate::error::SubdivError::IndexOutOfRange`] for the first offender.
    pub fn validate(&self) -> Result<(), crate::error::SubdivError> {
        let limit = self.vertices.len();
        for face in &self.faces {
            for &v in face {
                if v >= limit {
                    return Err(crate::error::SubdivError::IndexOutOfRange {
                        kind: "vertex",
                        idx: v,
                        limit,
                    });
                }
            }
        }
        Ok(())
    }

    /// Unit-cube polygon mesh — six quads. Useful starting point
    /// for tests; the truck-modeling kernel is not pulled in just
    /// for this.
    pub fn unit_cube() -> Self {
        let v = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let f = vec![
            vec![0, 3, 2, 1], // bottom (-Z)
            vec![4, 5, 6, 7], // top    (+Z)
            vec![0, 1, 5, 4], // front  (-Y)
            vec![1, 2, 6, 5], // right  (+X)
            vec![2, 3, 7, 6], // back   (+Y)
            vec![3, 0, 4, 7], // left   (-X)
        ];
        Self { vertices: v, faces: f }
    }

    /// Unit-tetrahedron triangle mesh — four triangles. Used by
    /// the Loop scheme tests.
    pub fn tetrahedron() -> Self {
        let v = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 1.0, 0.0),
            Vector3::new(0.5, 0.5, 1.0),
        ];
        let f = vec![
            vec![0, 1, 2],
            vec![0, 2, 3],
            vec![0, 3, 1],
            vec![1, 3, 2],
        ];
        Self { vertices: v, faces: f }
    }

    /// Undirected edge set with consistent ordering (low, high).
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut out = std::collections::BTreeSet::<(usize, usize)>::new();
        for face in &self.faces {
            for k in 0..face.len() {
                let a = face[k];
                let b = face[(k + 1) % face.len()];
                let key = if a < b { (a, b) } else { (b, a) };
                out.insert(key);
            }
        }
        out.into_iter().collect()
    }

    /// Per-vertex one-ring neighbours.
    pub fn vertex_one_ring(&self) -> Vec<std::collections::BTreeSet<usize>> {
        let mut out = vec![std::collections::BTreeSet::<usize>::new(); self.vertices.len()];
        for face in &self.faces {
            for k in 0..face.len() {
                let a = face[k];
                let b = face[(k + 1) % face.len()];
                out[a].insert(b);
                out[b].insert(a);
            }
        }
        out
    }
}
