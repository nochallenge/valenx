//! UV-aware quadric construction.
//!
//! Extends each vertex quadric with a 2x2 UV-space term so an edge
//! collapse that would stretch the parametrisation pays a higher cost.
//!
//! The matrix shape kept around is the canonical 4x4 spatial Q (set
//! to zero here since the spatial part is owned by valenx-mesh) plus a
//! 2x2 UV quadric, packed into a flat 6-tuple. v1 uses this for the
//! per-vertex UV-stretch penalty exported via [`uv_stretch_weight`]
//! that the [`crate::decimate_pro`] driver multiplies into the
//! curvature weight.

use nalgebra::{Matrix2, Vector2};

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

/// Wrapper around the per-vertex 2x2 UV quadric used by the UV-aware
/// driver. Stored as a packed `[a, b, c]` upper triangle of a 2x2
/// symmetric matrix so we don't drag nalgebra-serialize.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct QuadricMatrix {
    /// Upper-triangular packed `[uu, uv, vv]` of the 2x2 symmetric
    /// UV-space quadric (units: UV-distance squared).
    pub uv: [f64; 3],
}

impl QuadricMatrix {
    /// Convert back into a dense 2x2 matrix.
    pub fn to_matrix(self) -> Matrix2<f64> {
        Matrix2::new(self.uv[0], self.uv[1], self.uv[1], self.uv[2])
    }

    /// Apply the UV quadric to the supplied UV point `p` and return
    /// the scalar quadric form `pT Q p`.
    pub fn quadric_form(self, p: Vector2<f64>) -> f64 {
        let m = self.to_matrix();
        (p.transpose() * m * p)[(0, 0)]
    }
}

/// Build a per-vertex UV-aware quadric. Each triangle contributes a
/// rank-1 outer product of its UV-stretch direction to the quadrics of
/// its three incident vertices, scaled by the triangle's spatial area
/// (so big tris dominate small ones — matching Hoppe's "appearance-
/// preserving" weighting).
pub fn uv_aware_quadric(mesh: &Mesh, uvs: &[[f64; 2]]) -> Vec<QuadricMatrix> {
    let n = mesh.nodes.len();
    let mut q = vec![QuadricMatrix::default(); n];
    if uvs.len() != n {
        return q; // caller surfaces the size-mismatch error.
    }
    for block in &mesh.element_blocks {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        for tri in block.connectivity.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            let (i, j, k) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if i >= n || j >= n || k >= n {
                continue;
            }
            let area = 0.5
                * (mesh.nodes[j] - mesh.nodes[i])
                    .cross(&(mesh.nodes[k] - mesh.nodes[i]))
                    .norm();
            // Edge vectors in UV space, scaled by area.
            let euv1 = Vector2::new(uvs[j][0] - uvs[i][0], uvs[j][1] - uvs[i][1]);
            let euv2 = Vector2::new(uvs[k][0] - uvs[i][0], uvs[k][1] - uvs[i][1]);
            add_outer(&mut q[i], euv1, area);
            add_outer(&mut q[i], euv2, area);
            add_outer(&mut q[j], euv1, area);
            add_outer(&mut q[k], euv2, area);
        }
    }
    q
}

/// Scalar weight a curvature/feature/uv-aware driver multiplies into
/// the QEM cost so high-stretch verts resist collapse. Returns
/// `‖Q_v‖_F` (Frobenius norm of the UV quadric).
pub fn uv_stretch_weight(qm: QuadricMatrix) -> f64 {
    let m = qm.to_matrix();
    let s: f64 = m.iter().map(|x| x * x).sum();
    s.sqrt()
}

fn add_outer(q: &mut QuadricMatrix, v: Vector2<f64>, scale: f64) {
    q.uv[0] += scale * v.x * v.x;
    q.uv[1] += scale * v.x * v.y;
    q.uv[2] += scale * v.y * v.y;
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_mesh::element::ElementBlock;

    #[test]
    fn empty_mesh_returns_no_quadrics() {
        let mesh = Mesh::new("e");
        let q = uv_aware_quadric(&mesh, &[]);
        assert!(q.is_empty());
    }

    #[test]
    fn uv_quadric_is_zero_on_degenerate_uvs() {
        let mut m = Mesh::new("d");
        m.nodes.extend_from_slice(&[
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ]);
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(b);
        m.recompute_stats();
        let uvs = vec![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
        let q = uv_aware_quadric(&m, &uvs);
        for entry in q {
            assert_eq!(entry.uv, [0.0, 0.0, 0.0]);
        }
    }

    #[test]
    fn uv_stretch_weight_grows_with_uv_span() {
        let mut m = Mesh::new("s");
        m.nodes.extend_from_slice(&[
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ]);
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(b);
        m.recompute_stats();
        let uvs_small = vec![[0.0, 0.0], [0.01, 0.0], [0.0, 0.01]];
        let uvs_large = vec![[0.0, 0.0], [10.0, 0.0], [0.0, 10.0]];
        let q1 = uv_aware_quadric(&m, &uvs_small);
        let q2 = uv_aware_quadric(&m, &uvs_large);
        assert!(uv_stretch_weight(q2[0]) > uv_stretch_weight(q1[0]));
    }
}
