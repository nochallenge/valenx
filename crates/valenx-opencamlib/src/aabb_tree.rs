//! AABB tree over a flat triangle list.
//!
//! Recursive median split on the longest axis. Stops splitting once a
//! leaf holds ≤ `LEAF_SIZE` triangles. Used by [`crate::cutter::DropCutter`]
//! and [`crate::cutter::PushCutter`] to skip triangles that can't
//! possibly affect the current query.

use nalgebra::Vector3;

use crate::triangle::Triangle;

/// Triangles per leaf — small enough for cache-friendly scan.
pub const LEAF_SIZE: usize = 8;

/// Node in the AABB tree — internal or leaf.
#[derive(Clone, Debug)]
pub enum Node {
    /// Internal node with left + right child indices.
    Internal {
        /// Tight AABB enclosing both children.
        min: Vector3<f64>,
        /// Max corner of the AABB.
        max: Vector3<f64>,
        /// Left subtree node index.
        left: usize,
        /// Right subtree node index.
        right: usize,
    },
    /// Leaf node — holds triangle indices into the flat list.
    Leaf {
        /// Min AABB corner.
        min: Vector3<f64>,
        /// Max AABB corner.
        max: Vector3<f64>,
        /// Indices into the triangle list passed to [`AabbTree::new`].
        triangles: Vec<usize>,
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
#[derive(Clone, Debug)]
pub struct AabbTree {
    /// Flat list of nodes — `nodes[0]` is the root if non-empty.
    pub nodes: Vec<Node>,
}

impl AabbTree {
    /// Build the tree from a triangle list.
    pub fn new(tris: &[Triangle]) -> Self {
        let mut nodes = Vec::new();
        if tris.is_empty() {
            return Self { nodes };
        }
        let indices: Vec<usize> = (0..tris.len()).collect();
        build(tris, indices, &mut nodes);
        Self { nodes }
    }

    /// Query: indices of every triangle whose AABB intersects the
    /// ray `(origin, direction)` — returned in arbitrary order.
    /// Direction NOT required to be unit length.
    pub fn intersect_ray(
        &self,
        origin: Vector3<f64>,
        direction: Vector3<f64>,
    ) -> Vec<usize> {
        let mut out = Vec::new();
        if self.nodes.is_empty() {
            return out;
        }
        let inv_dir = Vector3::new(
            if direction.x.abs() > 1e-18 {
                1.0 / direction.x
            } else {
                f64::INFINITY
            },
            if direction.y.abs() > 1e-18 {
                1.0 / direction.y
            } else {
                f64::INFINITY
            },
            if direction.z.abs() > 1e-18 {
                1.0 / direction.z
            } else {
                f64::INFINITY
            },
        );
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            match &self.nodes[i] {
                Node::Internal {
                    min,
                    max,
                    left,
                    right,
                } => {
                    if ray_aabb(*min, *max, origin, inv_dir) {
                        stack.push(*left);
                        stack.push(*right);
                    }
                }
                Node::Leaf {
                    min,
                    max,
                    triangles,
                } => {
                    if ray_aabb(*min, *max, origin, inv_dir) {
                        out.extend_from_slice(triangles);
                    }
                }
            }
        }
        out
    }

    /// Query: indices of every triangle whose XY footprint contains
    /// `(x, y)` — used by [`crate::cutter::DropCutter::drop`] to enumerate
    /// candidate triangles for the Z evaluation.
    pub fn xy_query(&self, x: f64, y: f64) -> Vec<usize> {
        let mut out = Vec::new();
        if self.nodes.is_empty() {
            return out;
        }
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let (min, max) = self.nodes[i].aabb();
            if x < min.x - 1e-12 || x > max.x + 1e-12 || y < min.y - 1e-12 || y > max.y + 1e-12 {
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

fn build(tris: &[Triangle], mut indices: Vec<usize>, nodes: &mut Vec<Node>) -> usize {
    let aabb = tris_aabb(tris, &indices);
    if indices.len() <= LEAF_SIZE {
        let idx = nodes.len();
        nodes.push(Node::Leaf {
            min: aabb.0,
            max: aabb.1,
            triangles: indices,
        });
        return idx;
    }
    // Split along longest axis.
    let extent = aabb.1 - aabb.0;
    let axis = if extent.x >= extent.y && extent.x >= extent.z {
        0
    } else if extent.y >= extent.z {
        1
    } else {
        2
    };
    // Sort by centroid on the split axis.
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

fn tris_aabb(tris: &[Triangle], indices: &[usize]) -> (Vector3<f64>, Vector3<f64>) {
    let mut lo = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &i in indices {
        let (l, h) = tris[i].aabb();
        lo = Vector3::new(lo.x.min(l.x), lo.y.min(l.y), lo.z.min(l.z));
        hi = Vector3::new(hi.x.max(h.x), hi.y.max(h.y), hi.z.max(h.z));
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
