//! Plane3d primitive — point + normal in 3D.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// A 3D plane defined by `normal . (x - origin) = 0`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Plane3d {
    /// Origin point on the plane.
    pub origin: Vector3<f64>,
    /// Unit normal.
    pub normal: Vector3<f64>,
}

impl Plane3d {
    /// XY plane at `z = z0`.
    pub fn xy_at(z0: f64) -> Self {
        Self {
            origin: Vector3::new(0.0, 0.0, z0),
            normal: Vector3::z(),
        }
    }

    /// XZ plane at `y = y0`.
    pub fn xz_at(y0: f64) -> Self {
        Self {
            origin: Vector3::new(0.0, y0, 0.0),
            normal: Vector3::y(),
        }
    }

    /// YZ plane at `x = x0`.
    pub fn yz_at(x0: f64) -> Self {
        Self {
            origin: Vector3::new(x0, 0.0, 0.0),
            normal: Vector3::x(),
        }
    }
}
