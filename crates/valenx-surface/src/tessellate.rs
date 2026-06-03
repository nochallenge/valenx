//! Tessellation — sample NURBS curves + surfaces into discrete
//! polylines + triangle meshes for the viewport.
//!
//! The mesh produced by [`surface`] is a Tri3 element block on a
//! `nu × nv` rectangular grid: every quad cell `(i, j)` becomes two
//! triangles (`(i, j) (i+1, j) (i+1, j+1)` and
//! `(i, j) (i+1, j+1) (i, j+1)`), so the total element count is
//! `2 (nu - 1) (nv - 1)` and the total node count is `nu * nv`.

use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};

use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;

/// Sample a NURBS curve into `n_segments + 1` evenly-spaced points
/// across its valid parameter range. Returns the points in order.
///
/// `n_segments` must be `>= 1`; smaller values are clamped up.
pub fn curve(curve: &NurbsCurve, n_segments: usize) -> Vec<Vector3<f64>> {
    let n_segments = n_segments.max(1);
    let (u_min, u_max) = curve.parameter_range();
    let mut out = Vec::with_capacity(n_segments + 1);
    for i in 0..=n_segments {
        let t = i as f64 / n_segments as f64;
        let u = u_min + t * (u_max - u_min);
        out.push(curve.evaluate(u));
    }
    out
}

/// Sample a NURBS surface at a `nu × nv` parameter grid and emit a
/// Tri3 [`valenx_mesh::Mesh`] suitable for the viewport.
///
/// Each parameter sample becomes one mesh node; each interior quad
/// cell is split into two triangles via the diagonal `(i,j)—(i+1,j+1)`.
/// `nu, nv >= 2` (smaller values are clamped up to 2 to produce a
/// degenerate-but-valid 1×1 quad).
pub fn surface(surface: &NurbsSurface, nu: usize, nv: usize) -> valenx_mesh::Mesh {
    let nu = nu.max(2);
    let nv = nv.max(2);
    let (u_min, u_max) = surface.u_range();
    let (v_min, v_max) = surface.v_range();

    let mut nodes = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let s = i as f64 / (nu - 1) as f64;
        let u = u_min + s * (u_max - u_min);
        for j in 0..nv {
            let t = j as f64 / (nv - 1) as f64;
            let v = v_min + t * (v_max - v_min);
            nodes.push(surface.evaluate(u, v));
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = Vec::with_capacity(2 * (nu - 1) * (nv - 1) * 3);
    let idx = |i: usize, j: usize| -> u32 { (i * nv + j) as u32 };
    for i in 0..(nu - 1) {
        for j in 0..(nv - 1) {
            let a = idx(i, j);
            let b = idx(i + 1, j);
            let c = idx(i + 1, j + 1);
            let d = idx(i, j + 1);
            // Triangle 1: a-b-c
            block.connectivity.push(a);
            block.connectivity.push(b);
            block.connectivity.push(c);
            // Triangle 2: a-c-d
            block.connectivity.push(a);
            block.connectivity.push(c);
            block.connectivity.push(d);
        }
    }

    let mut mesh = valenx_mesh::Mesh::new("nurbs-surface");
    mesh.nodes = nodes;
    mesh.element_blocks = vec![block];
    mesh.recompute_stats();
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic_bezier_curve(cps: [Vector3<f64>; 4]) -> NurbsCurve {
        NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            cps.to_vec(),
            vec![1.0; 4],
        )
        .unwrap()
    }

    fn planar_unit_square_surface() -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        Vector3::new(u, v, 0.0)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    #[test]
    fn curve_sampling_includes_endpoints() {
        let c = cubic_bezier_curve([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        let pts = curve(&c, 10);
        assert_eq!(pts.len(), 11);
        assert!((pts[0] - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-10);
        assert!((pts[10] - Vector3::new(3.0, 0.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn curve_clamps_zero_segments() {
        let c = cubic_bezier_curve([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        let pts = curve(&c, 0);
        // Clamped up to 1 segment → 2 points.
        assert_eq!(pts.len(), 2);
    }

    #[test]
    fn surface_tessellation_is_planar_with_correct_counts() {
        let s = planar_unit_square_surface();
        let m = surface(&s, 5, 7);
        // 5 * 7 = 35 nodes.
        assert_eq!(m.nodes.len(), 35);
        assert_eq!(m.stats.node_count, 35);
        // 2 * (5-1) * (7-1) = 48 triangles.
        assert_eq!(m.total_elements(), 48);
        assert_eq!(m.stats.element_count, 48);
        // Every node has z = 0 because the surface is planar.
        for n in &m.nodes {
            assert!(n.z.abs() < 1e-10, "node has z = {}", n.z);
        }
        // Connectivity has 3 indices per triangle.
        let block = &m.element_blocks[0];
        assert_eq!(block.element_type, ElementType::Tri3);
        assert_eq!(block.connectivity.len(), 48 * 3);
    }
}
