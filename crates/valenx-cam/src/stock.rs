//! Stock: the rectangular block the tool starts from and cuts down
//! into.
//!
//! v1 ships rectangular stock only — origin (corner) + extent in
//! each axis. The top face is at `origin.z + size.z`; ops above the
//! stock should rapid up to `top_z() + safe_z_clearance`.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;

/// A rectangular block of material. Origin is the corner with the
/// smallest XYZ; size is the extent along each axis (must be > 0).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Stock {
    /// World-space corner at minimum XYZ.
    pub origin: Vector3<f64>,
    /// Extent along each axis (must each be > 0).
    pub size: Vector3<f64>,
    /// Free-form material descriptor (e.g. `"6061-T6 aluminum"`).
    pub material: String,
}

impl Default for Stock {
    fn default() -> Self {
        Self {
            origin: Vector3::zeros(),
            size: Vector3::new(100.0, 100.0, 10.0),
            material: "aluminum".into(),
        }
    }
}

impl Stock {
    /// Construct a validated stock. Returns
    /// [`CamError::EmptyStock`] if any extent is ≤ 0.
    pub fn new(
        origin: Vector3<f64>,
        size: Vector3<f64>,
        material: impl Into<String>,
    ) -> Result<Self, CamError> {
        if !(size.x > 0.0 && size.y > 0.0 && size.z > 0.0) {
            return Err(CamError::EmptyStock);
        }
        Ok(Self {
            origin,
            size,
            material: material.into(),
        })
    }

    /// Z coordinate of the top face — the safe entry plane for ops.
    pub fn top_z(&self) -> f64 {
        self.origin.z + self.size.z
    }

    /// Z coordinate of the bottom face.
    pub fn bottom_z(&self) -> f64 {
        self.origin.z
    }

    /// Min / max corners of the axis-aligned bounding box.
    pub fn aabb(&self) -> (Vector3<f64>, Vector3<f64>) {
        (self.origin, self.origin + self.size)
    }

    /// Returns `true` if `xy` lies inside the stock's XY footprint
    /// (Z is ignored). Used by ops that need to stay over the stock.
    pub fn contains_xy(&self, x: f64, y: f64) -> bool {
        let (min, max) = self.aabb();
        x >= min.x && x <= max.x && y >= min.y && y <= max.y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_100x100x10() {
        let s = Stock::default();
        assert_eq!(s.size, Vector3::new(100.0, 100.0, 10.0));
        assert!((s.top_z() - 10.0).abs() < 1e-12);
        assert!((s.bottom_z() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn empty_extent_rejected() {
        let bad = Stock::new(Vector3::zeros(), Vector3::new(10.0, 0.0, 5.0), "");
        assert!(matches!(bad, Err(CamError::EmptyStock)));
    }

    #[test]
    fn aabb_is_origin_plus_size() {
        let s = Stock::new(
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(10.0, 20.0, 30.0),
            "x",
        )
        .unwrap();
        let (min, max) = s.aabb();
        assert_eq!(min, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(max, Vector3::new(11.0, 22.0, 33.0));
    }

    #[test]
    fn contains_xy_works() {
        let s = Stock::default();
        assert!(s.contains_xy(50.0, 50.0));
        assert!(!s.contains_xy(150.0, 50.0));
    }
}
