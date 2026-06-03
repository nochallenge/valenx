//! Octree spatial index for fast surface queries.
//!
//! Stores triangle indices in cubic buckets. Build cost is `O(n log n)`,
//! query cost is `O(log n + k)` where `k` is the number of triangles
//! actually returned. Used by [`crate::cutter::AdaptiveDropCutter`] to
//! decide whether a region needs subdivision.

use nalgebra::Vector3;

use crate::triangle::Triangle;

/// Max triangles per leaf before subdivision.
pub const MAX_PER_LEAF: usize = 16;
/// Max recursion depth (safety cap).
pub const MAX_DEPTH: u32 = 12;

/// Octree node.
#[derive(Clone, Debug)]
pub enum OctNode {
    /// Leaf with triangle indices.
    Leaf {
        /// Min AABB corner.
        min: Vector3<f64>,
        /// Max AABB corner.
        max: Vector3<f64>,
        /// Indices into the triangle list.
        triangles: Vec<usize>,
    },
    /// Internal node with 8 child indices in `[XYZ]` octant order.
    Internal {
        /// Min AABB corner.
        min: Vector3<f64>,
        /// Max AABB corner.
        max: Vector3<f64>,
        /// 8 child node indices.
        children: [usize; 8],
    },
}

impl OctNode {
    /// AABB of this node.
    pub fn aabb(&self) -> (Vector3<f64>, Vector3<f64>) {
        match self {
            OctNode::Leaf { min, max, .. } | OctNode::Internal { min, max, .. } => (*min, *max),
        }
    }
}

/// Octree.
#[derive(Clone, Debug)]
pub struct Octree {
    /// All nodes; `nodes[0]` is the root when present.
    pub nodes: Vec<OctNode>,
}

impl Octree {
    /// Build the octree from a triangle list.
    pub fn new(tris: &[Triangle]) -> Self {
        let mut nodes = Vec::new();
        if tris.is_empty() {
            return Self { nodes };
        }
        let mut lo = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for t in tris {
            let (l, h) = t.aabb();
            lo = Vector3::new(lo.x.min(l.x), lo.y.min(l.y), lo.z.min(l.z));
            hi = Vector3::new(hi.x.max(h.x), hi.y.max(h.y), hi.z.max(h.z));
        }
        let indices: Vec<usize> = (0..tris.len()).collect();
        build(tris, lo, hi, indices, 0, &mut nodes);
        Self { nodes }
    }

    /// Triangles overlapping the XY column at `(x, y)`.
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
                OctNode::Internal { children, .. } => {
                    stack.extend_from_slice(children);
                }
                OctNode::Leaf { triangles, .. } => {
                    out.extend_from_slice(triangles);
                }
            }
        }
        out
    }
}

fn build(
    tris: &[Triangle],
    min: Vector3<f64>,
    max: Vector3<f64>,
    indices: Vec<usize>,
    depth: u32,
    nodes: &mut Vec<OctNode>,
) -> usize {
    if indices.len() <= MAX_PER_LEAF || depth >= MAX_DEPTH {
        let idx = nodes.len();
        nodes.push(OctNode::Leaf {
            min,
            max,
            triangles: indices,
        });
        return idx;
    }
    let mid = (min + max) * 0.5;
    // Bucket every triangle by which child octant its centroid lies in.
    let mut buckets: [Vec<usize>; 8] = Default::default();
    for &i in &indices {
        let c = (tris[i].v[0] + tris[i].v[1] + tris[i].v[2]) / 3.0;
        let oct = (if c.x > mid.x { 1 } else { 0 })
            | (if c.y > mid.y { 2 } else { 0 })
            | (if c.z > mid.z { 4 } else { 0 });
        buckets[oct].push(i);
    }
    let placeholder = nodes.len();
    nodes.push(OctNode::Leaf {
        min,
        max,
        triangles: Vec::new(),
    });
    let mut children = [0usize; 8];
    for (k, bucket) in buckets.into_iter().enumerate() {
        let cmin = Vector3::new(
            if k & 1 != 0 { mid.x } else { min.x },
            if k & 2 != 0 { mid.y } else { min.y },
            if k & 4 != 0 { mid.z } else { min.z },
        );
        let cmax = Vector3::new(
            if k & 1 != 0 { max.x } else { mid.x },
            if k & 2 != 0 { max.y } else { mid.y },
            if k & 4 != 0 { max.z } else { mid.z },
        );
        children[k] = build(tris, cmin, cmax, bucket, depth + 1, nodes);
    }
    nodes[placeholder] = OctNode::Internal { min, max, children };
    placeholder
}
