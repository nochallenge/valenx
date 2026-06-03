//! Triangle mesh ↔ triangle mesh collision.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::aabb::{intersect, Aabb};

/// Description of a detected collision between two meshes.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CollisionInfo {
    /// Approximate contact point (centroid of the first intersecting
    /// triangle pair).
    pub contact_point: Vector3<f64>,
    /// Triangle index in mesh A (linear across Tri3 blocks).
    pub tri_a: usize,
    /// Triangle index in mesh B.
    pub tri_b: usize,
}

/// AABB prune + per-pair triangle test. Returns the first detected
/// collision (or `None`).
pub fn collide(m1: &Mesh, m2: &Mesh) -> Option<CollisionInfo> {
    let bb1 = bbox_of(m1)?;
    let bb2 = bbox_of(m2)?;
    if !intersect(&bb1, &bb2) {
        return None;
    }

    let tris1 = collect_tris(m1);
    let tris2 = collect_tris(m2);
    for (i, t1) in tris1.iter().enumerate() {
        let bb_t1 = Aabb::from_points(t1.iter());
        if !intersect(&bb_t1, &bb2) {
            continue;
        }
        for (j, t2) in tris2.iter().enumerate() {
            let bb_t2 = Aabb::from_points(t2.iter());
            if !intersect(&bb_t1, &bb_t2) {
                continue;
            }
            if tri_tri_test(t1, t2) {
                let cp = (t1[0] + t1[1] + t1[2] + t2[0] + t2[1] + t2[2]) / 6.0;
                return Some(CollisionInfo {
                    contact_point: cp,
                    tri_a: i,
                    tri_b: j,
                });
            }
        }
    }
    None
}

/// Bounding box over every node in the mesh.
pub fn bbox_of(m: &Mesh) -> Option<Aabb> {
    if m.nodes.is_empty() {
        return None;
    }
    Some(Aabb::from_points(m.nodes.iter()))
}

fn collect_tris(m: &Mesh) -> Vec<[Vector3<f64>; 3]> {
    let mut out = Vec::new();
    for block in &m.element_blocks {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        for chunk in block.connectivity.chunks(3) {
            let a = m.nodes[chunk[0] as usize];
            let b = m.nodes[chunk[1] as usize];
            let c = m.nodes[chunk[2] as usize];
            out.push([a, b, c]);
        }
    }
    out
}

/// Möller separating-axis triangle-triangle test (overlap version,
/// without coplanar case). Approximates by testing all 9 edge-cross
/// axes + the two face normals.
fn tri_tri_test(t1: &[Vector3<f64>; 3], t2: &[Vector3<f64>; 3]) -> bool {
    let n1 = (t1[1] - t1[0]).cross(&(t1[2] - t1[0]));
    let n2 = (t2[1] - t2[0]).cross(&(t2[2] - t2[0]));
    if !overlap_on_axis(t1, t2, &n1) {
        return false;
    }
    if !overlap_on_axis(t1, t2, &n2) {
        return false;
    }
    let edges1 = [t1[1] - t1[0], t1[2] - t1[1], t1[0] - t1[2]];
    let edges2 = [t2[1] - t2[0], t2[2] - t2[1], t2[0] - t2[2]];
    for ea in &edges1 {
        for eb in &edges2 {
            let axis = ea.cross(eb);
            if axis.norm_squared() < 1e-12 {
                continue;
            }
            if !overlap_on_axis(t1, t2, &axis) {
                return false;
            }
        }
    }
    true
}

fn overlap_on_axis(
    t1: &[Vector3<f64>; 3],
    t2: &[Vector3<f64>; 3],
    axis: &Vector3<f64>,
) -> bool {
    let (a_lo, a_hi) = project(t1, axis);
    let (b_lo, b_hi) = project(t2, axis);
    !(a_hi < b_lo || b_hi < a_lo)
}

fn project(t: &[Vector3<f64>; 3], axis: &Vector3<f64>) -> (f64, f64) {
    let p0 = t[0].dot(axis);
    let p1 = t[1].dot(axis);
    let p2 = t[2].dot(axis);
    let lo = p0.min(p1).min(p2);
    let hi = p0.max(p1).max(p2);
    (lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_mesh::element::ElementBlock;

    fn single_tri(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Mesh {
        let mut m = Mesh::new("tri");
        m.nodes.push(a);
        m.nodes.push(b);
        m.nodes.push(c);
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(blk);
        m.recompute_stats();
        m
    }

    #[test]
    fn coincident_tris_collide() {
        let a = single_tri(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        );
        let b = single_tri(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        );
        assert!(collide(&a, &b).is_some());
    }

    #[test]
    fn disjoint_tris_do_not_collide() {
        let a = single_tri(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        );
        let b = single_tri(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(5.0, 1.0, 0.0),
        );
        assert!(collide(&a, &b).is_none());
    }

    #[test]
    fn crossing_tris_collide() {
        // Two tris that cross in the +x axis.
        let a = single_tri(
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        );
        let b = single_tri(
            Vector3::new(0.0, -1.0, 0.5),
            Vector3::new(0.0, 1.0, 0.5),
            Vector3::new(0.0, 0.0, -0.5),
        );
        assert!(collide(&a, &b).is_some());
    }
}
