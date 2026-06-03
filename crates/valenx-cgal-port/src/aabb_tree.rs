//! AABB tree over 3D triangles — analogous to `valenx-opencamlib`'s
//! tree but exposes the **ray-intersection-with-triangle** API CGAL
//! users expect.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Stable id into the triangle list.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct TriangleId(pub usize);

/// A single triangle.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Triangle3 {
    /// Three vertices.
    pub v: [Vector3<f64>; 3],
}

/// AABB tree node.
#[derive(Clone, Debug)]
pub enum Node {
    /// Internal node — two child indices.
    Internal {
        /// AABB min.
        min: Vector3<f64>,
        /// AABB max.
        max: Vector3<f64>,
        /// Left child index.
        left: usize,
        /// Right child index.
        right: usize,
    },
    /// Leaf — list of triangle ids.
    Leaf {
        /// AABB min.
        min: Vector3<f64>,
        /// AABB max.
        max: Vector3<f64>,
        /// Triangle indices.
        triangles: Vec<TriangleId>,
    },
}

impl Node {
    /// AABB of this node.
    pub fn aabb(&self) -> (Vector3<f64>, Vector3<f64>) {
        match self {
            Node::Internal { min, max, .. } | Node::Leaf { min, max, .. } => (*min, *max),
        }
    }
}

/// AABB tree.
pub struct AabbTree {
    /// Triangle list (owned).
    pub triangles: Vec<Triangle3>,
    /// Flat node list.
    pub nodes: Vec<Node>,
}

impl AabbTree {
    /// Build from an owned triangle list.
    pub fn new(triangles: Vec<Triangle3>) -> Self {
        let mut nodes = Vec::new();
        if triangles.is_empty() {
            return Self { triangles, nodes };
        }
        let indices: Vec<usize> = (0..triangles.len()).collect();
        build(&triangles, indices, &mut nodes);
        Self { triangles, nodes }
    }

    /// All triangles whose AABB intersects the ray `(origin, direction)`.
    pub fn intersect_ray(
        &self,
        origin: Vector3<f64>,
        direction: Vector3<f64>,
    ) -> Vec<TriangleId> {
        let mut out = Vec::new();
        if self.nodes.is_empty() {
            return out;
        }
        let inv = Vector3::new(
            inv(direction.x),
            inv(direction.y),
            inv(direction.z),
        );
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let (min, max) = self.nodes[i].aabb();
            if !ray_aabb(min, max, origin, inv) {
                continue;
            }
            match &self.nodes[i] {
                Node::Internal { left, right, .. } => {
                    stack.push(*left);
                    stack.push(*right);
                }
                Node::Leaf { triangles, .. } => {
                    out.extend_from_slice(triangles);
                }
            }
        }
        out
    }
}

fn inv(v: f64) -> f64 {
    if v.abs() > 1e-18 {
        1.0 / v
    } else {
        f64::INFINITY
    }
}

fn build(tris: &[Triangle3], mut indices: Vec<usize>, nodes: &mut Vec<Node>) -> usize {
    let aabb = tris_aabb(tris, &indices);
    if indices.len() <= 4 {
        let idx = nodes.len();
        nodes.push(Node::Leaf {
            min: aabb.0,
            max: aabb.1,
            triangles: indices.into_iter().map(TriangleId).collect(),
        });
        return idx;
    }
    let extent = aabb.1 - aabb.0;
    let axis = if extent.x >= extent.y && extent.x >= extent.z {
        0
    } else if extent.y >= extent.z {
        1
    } else {
        2
    };
    indices.sort_by(|a, b| {
        let ca = (tris[*a].v[0][axis] + tris[*a].v[1][axis] + tris[*a].v[2][axis]) / 3.0;
        let cb = (tris[*b].v[0][axis] + tris[*b].v[1][axis] + tris[*b].v[2][axis]) / 3.0;
        ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
    });
    let mid = indices.len() / 2;
    let right = indices.split_off(mid);
    let placeholder = nodes.len();
    nodes.push(Node::Leaf {
        min: aabb.0,
        max: aabb.1,
        triangles: Vec::new(),
    });
    let left_idx = build(tris, indices, nodes);
    let right_idx = build(tris, right, nodes);
    nodes[placeholder] = Node::Internal {
        min: aabb.0,
        max: aabb.1,
        left: left_idx,
        right: right_idx,
    };
    placeholder
}

fn tris_aabb(tris: &[Triangle3], indices: &[usize]) -> (Vector3<f64>, Vector3<f64>) {
    let mut lo = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &i in indices {
        for p in &tris[i].v {
            lo = Vector3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vector3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    (lo, hi)
}

fn ray_aabb(
    min: Vector3<f64>,
    max: Vector3<f64>,
    origin: Vector3<f64>,
    inv_dir: Vector3<f64>,
) -> bool {
    let t1 = (min.x - origin.x) * inv_dir.x;
    let t2 = (max.x - origin.x) * inv_dir.x;
    let t3 = (min.y - origin.y) * inv_dir.y;
    let t4 = (max.y - origin.y) * inv_dir.y;
    let t5 = (min.z - origin.z) * inv_dir.z;
    let t6 = (max.z - origin.z) * inv_dir.z;
    let tmin = t1.min(t2).max(t3.min(t4)).max(t5.min(t6));
    let tmax = t1.max(t2).min(t3.max(t4)).min(t5.max(t6));
    tmax >= tmin.max(0.0)
}
