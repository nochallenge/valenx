//! Triangle primitive + helpers used by every cutter algorithm.
//!
//! OpenCamLib works on flat triangle soups. We mirror that to keep the
//! Drop / Push / Edge cutter code small — callers convert from
//! [`valenx_mesh::Mesh`] into [`Triangle`] vectors at the boundary.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Single triangle in 3D world space.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Triangle {
    /// Three vertices in CCW order (when viewed from outside).
    pub v: [Vector3<f64>; 3],
}

impl Triangle {
    /// Construct from three points.
    pub fn new(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Self {
        Self { v: [a, b, c] }
    }

    /// Axis-aligned bounding box `(min, max)`.
    pub fn aabb(&self) -> (Vector3<f64>, Vector3<f64>) {
        let mut lo = self.v[0];
        let mut hi = self.v[0];
        for p in &self.v[1..] {
            lo.x = lo.x.min(p.x);
            lo.y = lo.y.min(p.y);
            lo.z = lo.z.min(p.z);
            hi.x = hi.x.max(p.x);
            hi.y = hi.y.max(p.y);
            hi.z = hi.z.max(p.z);
        }
        (lo, hi)
    }

    /// Triangle normal (NOT unit-length).
    pub fn normal(&self) -> Vector3<f64> {
        (self.v[1] - self.v[0]).cross(&(self.v[2] - self.v[0]))
    }

    /// Z value at `(x, y)` via barycentric interpolation. Returns
    /// `None` if `(x, y)` is outside the triangle's XY projection or
    /// the triangle is degenerate.
    pub fn z_at_xy(&self, x: f64, y: f64) -> Option<f64> {
        let (a, b, c) = (self.v[0], self.v[1], self.v[2]);
        let denom = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
        if denom.abs() < 1e-18 {
            return None;
        }
        let l1 = ((b.y - c.y) * (x - c.x) + (c.x - b.x) * (y - c.y)) / denom;
        let l2 = ((c.y - a.y) * (x - c.x) + (a.x - c.x) * (y - c.y)) / denom;
        let l3 = 1.0 - l1 - l2;
        if l1 < -1e-9 || l2 < -1e-9 || l3 < -1e-9 {
            return None;
        }
        Some(l1 * a.z + l2 * b.z + l3 * c.z)
    }
}

/// Convert a [`valenx_mesh::Mesh`] of Tri3 elements into a flat
/// triangle list. Non-Tri3 blocks are skipped.
pub fn from_valenx_mesh(m: &valenx_mesh::Mesh) -> Vec<Triangle> {
    let mut out = Vec::new();
    for block in &m.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = m.nodes[tri[0] as usize];
            let b = m.nodes[tri[1] as usize];
            let c = m.nodes[tri[2] as usize];
            out.push(Triangle::new(a, b, c));
        }
    }
    out
}
