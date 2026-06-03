//! Phase 190 — `AIS_Animation_Object::SetLocation` — animate an
//! object's pose between two transforms.
//!
//! ## What OCCT does
//!
//! `AIS_Animation_Object` carries a start + end `gp_Trsf` (4×4 rigid
//! transform). On each `Update(time)` it interpolates the translation
//! component linearly and the rotation component as a slerp on the
//! corresponding unit quaternion, then composes them back into a
//! `gp_Trsf` it applies via `SetLocalTransformation` to the AIS
//! interactive object.
//!
//! ## v1 status
//!
//! **Honest v1.** Returns a 4×4 [`nalgebra::Matrix4`] computed from
//! lerp'd translation + slerp'd rotation. Scale is not supported (OCCT
//! `gp_Trsf` is rigid-only; AIS doesn't expose scale either). The
//! caller is responsible for plumbing the resulting matrix into their
//! per-object model matrix — Valenx's wgpu pipeline already accepts a
//! per-object MVP uniform so this drops in cleanly.

use nalgebra::{Matrix3, Matrix4, UnitQuaternion, Vector3};

use crate::error::OcctVizError;

/// Rigid pose: rotation (quaternion) + translation.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pose {
    /// Rotation as a unit quaternion `(w, x, y, z)`.
    pub rotation: UnitQuaternion<f32>,
    /// World-space translation.
    pub translation: Vector3<f32>,
}

impl Pose {
    /// Identity pose — no rotation, zero translation.
    pub fn identity() -> Self {
        Self {
            rotation: UnitQuaternion::identity(),
            translation: Vector3::zeros(),
        }
    }
}

/// Interpolate `start → end` at parameter `t ∈ [0, 1]`. Translation
/// lerps; rotation slerps (shortest arc).
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `t` is not finite.
pub fn view_animation_object_motion(
    start: &Pose,
    end: &Pose,
    t: f32,
) -> Result<Matrix4<f32>, OcctVizError> {
    if !t.is_finite() {
        return Err(OcctVizError::bad_input("t", "must be finite"));
    }
    let t = t.clamp(0.0, 1.0);

    let r = start.rotation.slerp(&end.rotation, t);
    let p = start.translation + (end.translation - start.translation) * t;

    // Assemble 4×4 from quaternion + translation.
    let m3: Matrix3<f32> = r.to_rotation_matrix().into_inner();
    let mut out = Matrix4::identity();
    for i in 0..3 {
        for j in 0..3 {
            out[(i, j)] = m3[(i, j)];
        }
    }
    out[(0, 3)] = p.x;
    out[(1, 3)] = p.y;
    out[(2, 3)] = p.z;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn rejects_nan_t() {
        let a = Pose::identity();
        let b = Pose::identity();
        let err = view_animation_object_motion(&a, &b, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn t_zero_returns_identity_for_identity_start() {
        let a = Pose::identity();
        let b = Pose {
            translation: Vector3::new(10.0, 0.0, 0.0),
            ..Pose::identity()
        };
        let m = view_animation_object_motion(&a, &b, 0.0).unwrap();
        // Translation column = (0, 0, 0).
        assert!((m[(0, 3)]).abs() < 1e-4);
        assert!((m[(1, 3)]).abs() < 1e-4);
        assert!((m[(2, 3)]).abs() < 1e-4);
    }

    #[test]
    fn t_one_reaches_end_translation() {
        let a = Pose::identity();
        let b = Pose {
            translation: Vector3::new(10.0, 5.0, -3.0),
            ..Pose::identity()
        };
        let m = view_animation_object_motion(&a, &b, 1.0).unwrap();
        assert!((m[(0, 3)] - 10.0).abs() < 1e-4);
        assert!((m[(1, 3)] - 5.0).abs() < 1e-4);
        assert!((m[(2, 3)] - (-3.0)).abs() < 1e-4);
    }

    #[test]
    fn midpoint_translation_is_halfway() {
        let a = Pose::identity();
        let b = Pose {
            translation: Vector3::new(10.0, 0.0, 0.0),
            ..Pose::identity()
        };
        let m = view_animation_object_motion(&a, &b, 0.5).unwrap();
        assert!((m[(0, 3)] - 5.0).abs() < 1e-4);
    }

    #[test]
    fn rotation_slerps_through_intermediate() {
        let a = Pose::identity();
        let b = Pose {
            rotation: UnitQuaternion::from_axis_angle(&Vector3::z_axis(), PI / 2.0),
            ..Pose::identity()
        };
        let m = view_animation_object_motion(&a, &b, 0.5).unwrap();
        // After 45° z-rotation, the (0,0) entry should be cos(45°) ≈ 0.707.
        assert!((m[(0, 0)] - (PI / 4.0).cos()).abs() < 1e-3);
    }
}
