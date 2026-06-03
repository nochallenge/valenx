//! 4/5-axis types and transforms (Phase 17E).
//!
//! Parallel types to the 3-axis [`crate::toolpath::Toolpath`] /
//! [`crate::toolpath::Move`] — keeping the 3-axis surface untouched
//! preserves Phase 10 callers verbatim.
//!
//! ## Convention
//!
//! - **`a`** — rotation about the world X axis (degrees).
//! - **`b`** — rotation about the world Y axis (degrees).
//!
//! Some controllers use different axis conventions (e.g. Sinumerik
//! uses `C` for the Z-axis rotary); the postprocessor is responsible
//! for relabelling at G-code emit time.

use nalgebra::{Matrix4, Vector3};
use serde::{Deserialize, Serialize};

use crate::toolpath::MoveKind;

/// A single 5-axis tool motion. Extends [`crate::toolpath::Move`]
/// with `a` / `b` rotary positions (degrees).
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Move5Ax {
    /// Motion kind.
    pub kind: MoveKind,
    /// Absolute target position (mm) in stock-local coordinates.
    pub position: Vector3<f64>,
    /// A-axis target (degrees).
    pub a_deg: f64,
    /// B-axis target (degrees).
    pub b_deg: f64,
    /// Feed (mm/min). Ignored for `Rapid`.
    pub feed: f64,
}

impl Move5Ax {
    /// Convenience constructor.
    pub fn new(kind: MoveKind, position: Vector3<f64>, a_deg: f64, b_deg: f64, feed: f64) -> Self {
        Self {
            kind,
            position,
            a_deg,
            b_deg,
            feed,
        }
    }
}

/// 5-axis toolpath. Parallel to [`crate::toolpath::Toolpath`] for
/// the 3-axis case.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Toolpath5Ax {
    /// All 5-axis moves in execution order.
    pub moves: Vec<Move5Ax>,
}

impl Toolpath5Ax {
    /// Empty toolpath.
    pub fn new() -> Self {
        Self { moves: Vec::new() }
    }

    /// Push a single move.
    pub fn push(&mut self, m: Move5Ax) {
        self.moves.push(m);
    }

    /// `true` when no moves are queued.
    pub fn is_empty(&self) -> bool {
        self.moves.is_empty()
    }

    /// Number of moves.
    pub fn len(&self) -> usize {
        self.moves.len()
    }
}

/// Build a 4×4 homogeneous rotation matrix for the given A/B (deg)
/// pose. Useful for transforming 3-axis output into a fixed-axis
/// indexed orientation.
pub fn rotation_matrix(a_deg: f64, b_deg: f64) -> Matrix4<f64> {
    let a = a_deg.to_radians();
    let b = b_deg.to_radians();
    let (sa, ca) = a.sin_cos();
    let (sb, cb) = b.sin_cos();
    // Rx(a) — rotate around X.
    let rx = Matrix4::new(
        1.0, 0.0, 0.0, 0.0, 0.0, ca, -sa, 0.0, 0.0, sa, ca, 0.0, 0.0, 0.0, 0.0, 1.0,
    );
    // Ry(b) — rotate around Y.
    let ry = Matrix4::new(
        cb, 0.0, sb, 0.0, 0.0, 1.0, 0.0, 0.0, -sb, 0.0, cb, 0.0, 0.0, 0.0, 0.0, 1.0,
    );
    // Apply Y first, then X (matches Fanuc 30i Z⇒B⇒A chain).
    rx * ry
}

/// Transform a homogeneous point `(x, y, z, 1)` by `m`.
pub fn transform_point(m: &Matrix4<f64>, p: Vector3<f64>) -> Vector3<f64> {
    let h = m * nalgebra::Vector4::new(p.x, p.y, p.z, 1.0);
    Vector3::new(h.x, h.y, h.z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn move5ax_construct() {
        let m = Move5Ax::new(
            MoveKind::Cut,
            Vector3::new(1.0, 2.0, 3.0),
            30.0,
            45.0,
            500.0,
        );
        assert_eq!(m.a_deg, 30.0);
        assert_eq!(m.b_deg, 45.0);
        assert_eq!(m.position.x, 1.0);
    }

    #[test]
    fn rotation_zero_is_identity() {
        let m = rotation_matrix(0.0, 0.0);
        let p = Vector3::new(1.0, 2.0, 3.0);
        let q = transform_point(&m, p);
        assert!((q - p).norm() < 1e-9);
    }

    #[test]
    fn rotation_90deg_a_swaps_yz() {
        let m = rotation_matrix(90.0, 0.0);
        let q = transform_point(&m, Vector3::new(0.0, 1.0, 0.0));
        assert!(q.x.abs() < 1e-9);
        assert!(q.y.abs() < 1e-9);
        assert!((q.z - 1.0).abs() < 1e-9);
    }

    #[test]
    fn toolpath5ax_basics() {
        let mut tp = Toolpath5Ax::new();
        assert!(tp.is_empty());
        tp.push(Move5Ax::new(
            MoveKind::Rapid,
            Vector3::zeros(),
            0.0,
            0.0,
            0.0,
        ));
        assert_eq!(tp.len(), 1);
    }
}
