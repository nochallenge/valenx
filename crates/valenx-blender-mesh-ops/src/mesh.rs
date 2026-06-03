//! Polygon mesh primitive shared by all Blender-style ops.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Polygon mesh.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Mesh {
    /// Vertex positions.
    pub vertices: Vec<Vector3<f64>>,
    /// Faces — vertex-index loops in CCW order.
    pub faces: Vec<Vec<usize>>,
}

impl Mesh {
    /// Empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// Unit cube as six quads.
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
            vec![0, 3, 2, 1], // bottom
            vec![4, 5, 6, 7], // top
            vec![0, 1, 5, 4], // front
            vec![1, 2, 6, 5], // right
            vec![2, 3, 7, 6], // back
            vec![3, 0, 4, 7], // left
        ];
        Self { vertices: v, faces: f }
    }

    /// Vertex count.
    pub fn n_verts(&self) -> usize {
        self.vertices.len()
    }

    /// Face count.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Undirected edges in sorted (low, high) order.
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut set = std::collections::BTreeSet::<(usize, usize)>::new();
        for f in &self.faces {
            for k in 0..f.len() {
                let a = f[k];
                let b = f[(k + 1) % f.len()];
                let key = if a < b { (a, b) } else { (b, a) };
                set.insert(key);
            }
        }
        set.into_iter().collect()
    }
}

/// Newell face normal.
pub(crate) fn face_normal(verts: &[Vector3<f64>], face: &[usize]) -> Vector3<f64> {
    let mut n: Vector3<f64> = Vector3::zeros();
    let k = face.len();
    for i in 0..k {
        let p = verts[face[i]];
        let q = verts[face[(i + 1) % k]];
        n.x += (p.y - q.y) * (p.z + q.z);
        n.y += (p.z - q.z) * (p.x + q.x);
        n.z += (p.x - q.x) * (p.y + q.y);
    }
    let len = n.norm();
    if len < 1e-12 {
        Vector3::<f64>::z()
    } else {
        n / len
    }
}

/// Face centroid.
pub(crate) fn face_centroid(verts: &[Vector3<f64>], face: &[usize]) -> Vector3<f64> {
    let s: Vector3<f64> = face.iter().map(|&i| verts[i]).sum();
    s / face.len() as f64
}
