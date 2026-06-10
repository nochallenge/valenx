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

    /// Triangle area `½·‖(v₁−v₀)×(v₂−v₀)‖` — half the magnitude of the non-unit
    /// [`normal`](Self::normal). A degenerate (collinear) triangle returns `0.0`.
    pub fn area(&self) -> f64 {
        self.normal().norm() * 0.5
    }

    /// Triangle perimeter — the sum of the three edge lengths
    /// `‖v₁−v₀‖ + ‖v₂−v₁‖ + ‖v₀−v₂‖`. A length, distinct from [`area`](Self::area).
    pub fn perimeter(&self) -> f64 {
        (self.v[1] - self.v[0]).norm()
            + (self.v[2] - self.v[1]).norm()
            + (self.v[0] - self.v[2]).norm()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn area_via_cross_product_magnitude() {
        // Legs 3 and 4 at a right angle in the XY plane → area = ½·3·4 = 6.
        let t = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
        );
        assert!((t.area() - 6.0).abs() < 1e-12);
        // Same legs in the YZ plane (out of the XY plane) → same area.
        let yz = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 4.0),
        );
        assert!((yz.area() - 6.0).abs() < 1e-12);
        // Collinear (degenerate) triangle → 0.
        let flat = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
        );
        assert_eq!(flat.area(), 0.0);
    }

    #[test]
    fn perimeter_sums_the_three_edges() {
        // 3-4-5 right triangle → perimeter 3 + 5 + 4 = 12.
        let t = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
        );
        assert!((t.perimeter() - 12.0).abs() < 1e-12);
        // Equilateral side 2 → perimeter 6; each edge < the perimeter.
        let eq = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(1.0, 3.0_f64.sqrt(), 0.0),
        );
        assert!((eq.perimeter() - 6.0).abs() < 1e-10);
        // Collinear (0,0,0)/(1,0,0)/(2,0,0) → edges 1 + 1 + 2 = 4.
        let flat = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
        );
        assert!((flat.perimeter() - 4.0).abs() < 1e-12);
    }
}
