//! Camera description.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Pinhole camera.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Camera {
    /// Eye position in world space.
    pub position: Vector3<f64>,
    /// Point the camera is looking at.
    pub target: Vector3<f64>,
    /// Approximate world up — used to derive the camera right axis.
    pub up: Vector3<f64>,
    /// Vertical field of view in radians.
    pub fov_v_rad: f64,
    /// Output image width in pixels.
    pub image_width: u32,
    /// Output image height in pixels.
    pub image_height: u32,
}

impl Default for Camera {
    /// Sensible default: 1920×1080, looking at the origin from
    /// (5, 5, 5) with 60° vertical fov.
    fn default() -> Self {
        Self {
            position: Vector3::new(5.0, 5.0, 5.0),
            target: Vector3::zeros(),
            up: Vector3::z(),
            fov_v_rad: 60f64.to_radians(),
            image_width: 1920,
            image_height: 1080,
        }
    }
}

impl Camera {
    /// Aspect ratio (width / height).
    pub fn aspect(&self) -> f64 {
        self.image_width as f64 / self.image_height.max(1) as f64
    }
}
