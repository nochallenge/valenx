//! Axis-aligned bounding box in model units (assumed metres unless
//! the containing `Geometry` specifies otherwise).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Axis-aligned bounding box.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct BoundingBox {
    pub min: Vector3<f64>,
    pub max: Vector3<f64>,
}

impl BoundingBox {
    /// Construct from two corners; normalised so `min` <= `max`
    /// componentwise.
    pub fn new(a: Vector3<f64>, b: Vector3<f64>) -> Self {
        let min = Vector3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z));
        let max = Vector3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z));
        Self { min, max }
    }

    /// Empty box at the origin. Useful as a sentinel.
    pub fn empty() -> Self {
        Self {
            min: Vector3::zeros(),
            max: Vector3::zeros(),
        }
    }

    /// Midpoint of the box (mean of `min` and `max`).
    pub fn center(&self) -> Vector3<f64> {
        (self.min + self.max) * 0.5
    }

    /// `max - min` — the per-axis extents.
    pub fn size(&self) -> Vector3<f64> {
        self.max - self.min
    }

    /// Length of the main diagonal `‖size()‖₂`.
    pub fn diagonal(&self) -> f64 {
        self.size().norm()
    }

    /// Volume `Δx · Δy · Δz`. Returns 0 for a degenerate box.
    pub fn volume(&self) -> f64 {
        let s = self.size();
        s.x * s.y * s.z
    }

    /// Extend to contain `other`. Returns the expanded box.
    pub fn union(mut self, other: &BoundingBox) -> Self {
        self.min = Vector3::new(
            self.min.x.min(other.min.x),
            self.min.y.min(other.min.y),
            self.min.z.min(other.min.z),
        );
        self.max = Vector3::new(
            self.max.x.max(other.max.x),
            self.max.y.max(other.max.y),
            self.max.z.max(other.max.z),
        );
        self
    }

    /// `true` if `point` lies within the closed box (inclusive of
    /// the boundary).
    pub fn contains(&self, point: &Vector3<f64>) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalisation() {
        let bb = BoundingBox::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 2.0));
        assert_eq!(bb.min, Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(bb.max, Vector3::new(1.0, 1.0, 2.0));
    }

    #[test]
    fn union() {
        let a = BoundingBox::new(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let b = BoundingBox::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::zeros());
        let u = a.union(&b);
        assert_eq!(u.min, Vector3::new(-1.0, -1.0, -1.0));
        assert_eq!(u.max, Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn contains_point() {
        let bb = BoundingBox::new(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        assert!(bb.contains(&Vector3::new(0.5, 0.5, 0.5)));
        assert!(!bb.contains(&Vector3::new(1.1, 0.0, 0.0)));
    }
}
