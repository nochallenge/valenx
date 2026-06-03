//! Axis-aligned bounding-box utilities.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// 3D axis-aligned bounding box.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Aabb {
    /// Minimum corner.
    pub min: Vector3<f64>,
    /// Maximum corner.
    pub max: Vector3<f64>,
}

impl Aabb {
    /// Empty (inverted-bounds) AABB.
    pub fn empty() -> Self {
        Self {
            min: Vector3::repeat(f64::INFINITY),
            max: Vector3::repeat(f64::NEG_INFINITY),
        }
    }

    /// AABB from a set of points.
    pub fn from_points<'a>(pts: impl Iterator<Item = &'a Vector3<f64>>) -> Self {
        let mut bb = Self::empty();
        for p in pts {
            bb.include(p);
        }
        bb
    }

    /// Expand to include `p`.
    pub fn include(&mut self, p: &Vector3<f64>) {
        self.min = self.min.zip_map(p, f64::min);
        self.max = self.max.zip_map(p, f64::max);
    }

    /// Centre + half-extents.
    pub fn center(&self) -> Vector3<f64> {
        0.5 * (self.min + self.max)
    }

    /// Returns `true` if the bounds form a valid non-empty box.
    pub fn is_valid(&self) -> bool {
        self.min.x <= self.max.x && self.min.y <= self.max.y && self.min.z <= self.max.z
    }
}

/// Fast overlap test — returns `true` when the two boxes share any
/// volume (touching counts as overlap).
pub fn intersect(a: &Aabb, b: &Aabb) -> bool {
    if !a.is_valid() || !b.is_valid() {
        return false;
    }
    a.min.x <= b.max.x
        && a.max.x >= b.min.x
        && a.min.y <= b.max.y
        && a.max.y >= b.min.y
        && a.min.z <= b.max.z
        && a.max.z >= b.min.z
}

/// Squared L2 separation between the two boxes. Returns 0 when
/// they overlap.
pub fn distance_squared(a: &Aabb, b: &Aabb) -> f64 {
    if intersect(a, b) {
        return 0.0;
    }
    let dx = (b.min.x - a.max.x).max(a.min.x - b.max.x).max(0.0);
    let dy = (b.min.y - a.max.y).max(a.min.y - b.max.y).max(0.0);
    let dz = (b.min.z - a.max.z).max(a.min.z - b.max.z).max(0.0);
    dx * dx + dy * dy + dz * dz
}

/// L2 separation between the two boxes. Returns 0 when they overlap.
pub fn distance(a: &Aabb, b: &Aabb) -> f64 {
    distance_squared(a, b).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersect_overlapping_boxes() {
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 10.0, 10.0),
        };
        let b = Aabb {
            min: Vector3::new(5.0, 5.0, 5.0),
            max: Vector3::new(15.0, 15.0, 15.0),
        };
        assert!(intersect(&a, &b));
        assert!((distance(&a, &b) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn intersect_disjoint_boxes() {
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Vector3::new(2.0, 0.0, 0.0),
            max: Vector3::new(3.0, 1.0, 1.0),
        };
        assert!(!intersect(&a, &b));
        assert!((distance(&a, &b) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn touching_counts_as_overlap() {
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Vector3::new(1.0, 0.0, 0.0),
            max: Vector3::new(2.0, 1.0, 1.0),
        };
        assert!(intersect(&a, &b));
    }
}
