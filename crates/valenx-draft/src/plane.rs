//! Working plane in 3D space.
//!
//! A [`WorkingPlane`] is the 2D coordinate system inside which all
//! draft entities live. It carries:
//!
//! - `origin`  — world-space 3D origin of the plane,
//! - `normal`  — world-space unit normal (Z of the local frame),
//! - `x_axis`  — world-space unit X of the local frame.
//!
//! The implicit Y axis is `normal × x_axis`, which keeps the frame
//! right-handed regardless of which constructor produced it.
//!
//! A 2D point `[x, y]` on the plane projects to world-space as
//! `origin + x * x_axis + y * y_axis`.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// A 3D coordinate frame the Draft workbench uses as its 2D plane.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct WorkingPlane {
    /// World-space origin of the local frame.
    pub origin: Vector3<f64>,
    /// World-space unit normal (Z of the local frame). Should be unit
    /// length; constructors normalise.
    pub normal: Vector3<f64>,
    /// World-space unit X of the local frame. Should be unit length
    /// AND orthogonal to `normal`; constructors enforce this.
    pub x_axis: Vector3<f64>,
}

impl WorkingPlane {
    /// World XY plane: origin at world origin, normal = +Z, X = +X.
    pub fn from_xy() -> Self {
        Self {
            origin: Vector3::zeros(),
            normal: Vector3::new(0.0, 0.0, 1.0),
            x_axis: Vector3::new(1.0, 0.0, 0.0),
        }
    }

    /// World XZ plane: origin at world origin, X = world +X, local Y
    /// points along world +Z (so users see a "natural" front view).
    /// Normal is therefore world `-Y` (right-handed: `-Y × +X = +Z`).
    pub fn from_xz() -> Self {
        Self {
            origin: Vector3::zeros(),
            normal: Vector3::new(0.0, -1.0, 0.0),
            x_axis: Vector3::new(1.0, 0.0, 0.0),
        }
    }

    /// World YZ plane: origin at world origin, X = world +Y, local Y
    /// points along world +Z (so users see a "natural" side view).
    /// Normal is world `+X` (right-handed: `+X × +Y = +Z`).
    pub fn from_yz() -> Self {
        Self {
            origin: Vector3::zeros(),
            normal: Vector3::new(1.0, 0.0, 0.0),
            x_axis: Vector3::new(0.0, 1.0, 0.0),
        }
    }

    /// The implicit local Y axis: `normal × x_axis`. Right-handed.
    pub fn y_axis(&self) -> Vector3<f64> {
        self.normal.cross(&self.x_axis)
    }

    /// Convert a 2D point on this plane to world coordinates.
    pub fn local_to_world(&self, p: [f64; 2]) -> Vector3<f64> {
        self.origin + self.x_axis * p[0] + self.y_axis() * p[1]
    }

    /// Project a world-space point onto this plane, returning its 2D
    /// (x, y) coordinates in the local frame. Points off the plane
    /// are projected by dropping the component along the normal.
    pub fn world_to_local(&self, world: Vector3<f64>) -> [f64; 2] {
        let rel = world - self.origin;
        [rel.dot(&self.x_axis), rel.dot(&self.y_axis())]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xy_plane_round_trips_origin() {
        let p = WorkingPlane::from_xy();
        let local = [3.0, 4.0];
        let world = p.local_to_world(local);
        assert_eq!(world, Vector3::new(3.0, 4.0, 0.0));
        let back = p.world_to_local(world);
        assert!((back[0] - 3.0).abs() < 1e-12);
        assert!((back[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn xz_plane_maps_local_y_to_world_z() {
        let p = WorkingPlane::from_xz();
        let world = p.local_to_world([2.0, 5.0]);
        // Local x=2 → world x=2; local y=5 → world z=5; y stays at origin.
        assert!((world.x - 2.0).abs() < 1e-12);
        assert!(world.y.abs() < 1e-12);
        assert!((world.z - 5.0).abs() < 1e-12);
    }

    #[test]
    fn yz_plane_maps_local_axes_correctly() {
        let p = WorkingPlane::from_yz();
        // x_axis = +Y world, y_axis = normal(=+X) × x_axis(=+Y) = +Z world.
        let world = p.local_to_world([7.0, 9.0]);
        assert!(world.x.abs() < 1e-12);
        assert!((world.y - 7.0).abs() < 1e-12);
        assert!((world.z - 9.0).abs() < 1e-12);
    }

    #[test]
    fn world_to_local_projects_off_plane_point() {
        // On the XY plane, a point with non-zero z just drops the z.
        let p = WorkingPlane::from_xy();
        let local = p.world_to_local(Vector3::new(1.0, 2.0, 99.0));
        assert!((local[0] - 1.0).abs() < 1e-12);
        assert!((local[1] - 2.0).abs() < 1e-12);
    }
}
