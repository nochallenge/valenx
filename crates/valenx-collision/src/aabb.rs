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

    /// Volume of the box — the product of the per-axis extents, each clamped to `≥ 0`. Returns
    /// `0.0` for an inverted or empty box (any axis where `max < min`).
    pub fn volume(&self) -> f64 {
        let dx = (self.max.x - self.min.x).max(0.0);
        let dy = (self.max.y - self.min.y).max(0.0);
        let dz = (self.max.z - self.min.z).max(0.0);
        dx * dy * dz
    }

    /// Surface area of the box — `2·(dx·dy + dy·dz + dz·dx)` with each extent clamped to `≥ 0`.
    /// Returns `0.0` for an inverted or empty box (any axis where `max < min`).
    pub fn surface_area(&self) -> f64 {
        let dx = (self.max.x - self.min.x).max(0.0);
        let dy = (self.max.y - self.min.y).max(0.0);
        let dz = (self.max.z - self.min.z).max(0.0);
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Space-diagonal length of the box — `√(dx²+dy²+dz²)` with each extent clamped to `≥ 0`.
    /// Returns `0.0` for an inverted or empty box (any axis where `max < min`).
    pub fn diagonal(&self) -> f64 {
        let dx = (self.max.x - self.min.x).max(0.0);
        let dy = (self.max.y - self.min.y).max(0.0);
        let dz = (self.max.z - self.min.z).max(0.0);
        (dx * dx + dy * dy + dz * dz).sqrt()
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

    #[test]
    fn volume_of_box() {
        // 10×20×30 box → 6000.
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 20.0, 30.0),
        };
        assert!((a.volume() - 6000.0).abs() < 1e-9);
        // Unit cube → 1.
        let unit = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        assert!((unit.volume() - 1.0).abs() < 1e-12);
        // Inverted box → 0.
        let inv = Aabb {
            min: Vector3::new(10.0, 10.0, 10.0),
            max: Vector3::new(5.0, 5.0, 5.0),
        };
        assert_eq!(inv.volume(), 0.0);
        // Flat (zero-height) box → 0.
        let flat = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 10.0, 0.0),
        };
        assert_eq!(flat.volume(), 0.0);
    }

    #[test]
    fn surface_area_of_box() {
        // 10×20×30 → 2·(200+600+300) = 2200.
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 20.0, 30.0),
        };
        assert!((a.surface_area() - 2200.0).abs() < 1e-9);
        // Unit cube → 6.
        let unit = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        assert!((unit.surface_area() - 6.0).abs() < 1e-12);
        // Inverted → 0.
        let inv = Aabb {
            min: Vector3::new(10.0, 10.0, 10.0),
            max: Vector3::new(5.0, 5.0, 5.0),
        };
        assert_eq!(inv.surface_area(), 0.0);
        // Flat plate 10×10×0 → 2·100 = 200 (two faces).
        let flat = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 10.0, 0.0),
        };
        assert!((flat.surface_area() - 200.0).abs() < 1e-9);
    }

    #[test]
    fn diagonal_of_box() {
        // 10×20×30 → √1400.
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 20.0, 30.0),
        };
        assert!((a.diagonal() - 1400.0_f64.sqrt()).abs() < 1e-9);
        // Unit cube → √3.
        let unit = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        assert!((unit.diagonal() - 3.0_f64.sqrt()).abs() < 1e-12);
        // Pythagorean 3×4×0 → exactly 5.
        let pyth = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(3.0, 4.0, 0.0),
        };
        assert!((pyth.diagonal() - 5.0).abs() < 1e-12);
        // Inverted → 0.
        let inv = Aabb {
            min: Vector3::new(10.0, 10.0, 10.0),
            max: Vector3::new(5.0, 5.0, 5.0),
        };
        assert_eq!(inv.diagonal(), 0.0);
    }
}
