//! Parts — a [`Part`] wraps a [`valenx_cad::Solid`] with a rigid-body
//! transform ([`PartTransform`]) and a `fixed` flag that tells the
//! solver "do not move this part" (the assembly equivalent of a
//! grounded body in multi-body dynamics).

use nalgebra::{UnitQuaternion, Vector3};
use serde::{Deserialize, Serialize};

/// Rigid-body transform composed of a translation and an orientation
/// (unit quaternion).
///
/// The quaternion convention matches `nalgebra::UnitQuaternion`:
/// `apply_point(p) = orientation * p + translation`. Pure rotations
/// (`apply_vector`) drop the translation; the orientation is applied
/// alone.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PartTransform {
    /// World-space translation applied after the rotation.
    pub translation: Vector3<f64>,
    /// World-space orientation as a unit quaternion. Right-multiplied
    /// onto local-frame vectors.
    pub orientation: UnitQuaternion<f64>,
}

impl PartTransform {
    /// Identity transform — no rotation, no translation. The default
    /// for any new [`Part`].
    pub fn identity() -> Self {
        Self {
            translation: Vector3::zeros(),
            orientation: UnitQuaternion::identity(),
        }
    }

    /// Apply the full transform to a local-frame point:
    /// `world = R * local + t`.
    pub fn apply_point(&self, p: Vector3<f64>) -> Vector3<f64> {
        self.orientation * p + self.translation
    }

    /// Apply the rotation only to a local-frame vector (no translation):
    /// `world = R * local`. Use this for direction vectors.
    pub fn apply_vector(&self, v: Vector3<f64>) -> Vector3<f64> {
        self.orientation * v
    }

    /// Return the inverse transform. Useful for "express a world-space
    /// point in the local frame of this part" computations.
    pub fn inverse(&self) -> Self {
        let inv_rot = self.orientation.inverse();
        Self {
            translation: -(inv_rot * self.translation),
            orientation: inv_rot,
        }
    }
}

impl Default for PartTransform {
    fn default() -> Self {
        Self::identity()
    }
}

/// Axis-aligned bounding box in world space.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct AABB {
    /// Lower corner.
    pub min: Vector3<f64>,
    /// Upper corner.
    pub max: Vector3<f64>,
}

impl AABB {
    /// Build an AABB from a min + max corner. Both inputs assumed to
    /// already form a well-ordered bounding box (no swap performed).
    pub fn new(min: Vector3<f64>, max: Vector3<f64>) -> Self {
        Self { min, max }
    }

    /// Empty / "no points" sentinel — `min` set to +infinity, `max` to
    /// -infinity so any subsequent `expand_point` produces a sensible
    /// box.
    pub fn empty() -> Self {
        Self {
            min: Vector3::repeat(f64::INFINITY),
            max: Vector3::repeat(f64::NEG_INFINITY),
        }
    }

    /// Expand the box to include `p`.
    pub fn expand_point(&mut self, p: Vector3<f64>) {
        self.min = self.min.zip_map(&p, f64::min);
        self.max = self.max.zip_map(&p, f64::max);
    }
}

/// One part in an assembly — a solid + its world-space transform +
/// the "fixed" anchor flag.
///
/// Parts carry their own id so mates / joints can reference them
/// without holding `&Part` references (which would force borrows that
/// the solver can't satisfy).
///
/// `solid` is intentionally not `Serialize` (truck BRep isn't
/// serializable). Persistence ([`crate::persist`]) stores a sentinel
/// marker for the solid and expects callers to re-attach the geometry
/// after loading — same pattern as Phase 5's TechDraw, which projects
/// from the currently-loaded solid in the app.
#[derive(Clone, Debug)]
pub struct Part {
    /// Stable id assigned by [`crate::Assembly::add_part`]. Unique
    /// within one assembly; not reused after delete.
    pub id: usize,
    /// Display name (shown in the scene tree).
    pub name: String,
    /// The underlying solid geometry.
    pub solid: valenx_cad::Solid,
    /// Current world-space pose.
    pub transform: PartTransform,
    /// When `true` the solver excludes this part from the pose vector
    /// (its DOFs are removed). "Ground" / "world frame" parts go here.
    pub fixed: bool,
}

impl Part {
    /// Build a new part at the identity transform, not fixed.
    pub fn new(id: usize, name: impl Into<String>, solid: valenx_cad::Solid) -> Self {
        Self {
            id,
            name: name.into(),
            solid,
            transform: PartTransform::identity(),
            fixed: false,
        }
    }

    /// World-space axis-aligned bounding box of the part — tessellates
    /// the solid (at a fairly coarse tolerance) and sweeps every
    /// vertex through the current transform.
    ///
    /// Returns the empty AABB if tessellation fails (mesh-backed solid
    /// with zero triangles, BRep tess errors). Callers that need a
    /// guaranteed-valid AABB should pre-validate the solid.
    pub fn bounding_box(&self) -> AABB {
        let mesh = match valenx_cad::solid_to_mesh(&self.solid, valenx_cad::DEFAULT_TESS_TOLERANCE)
        {
            Ok(m) => m,
            Err(_) => return AABB::empty(),
        };
        let mut aabb = AABB::empty();
        for n in &mesh.nodes {
            aabb.expand_point(self.transform.apply_point(*n));
        }
        aabb
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn identity_transform_is_a_no_op() {
        let t = PartTransform::identity();
        let p = Vector3::new(3.0, 4.0, 5.0);
        assert!((t.apply_point(p) - p).norm() < 1e-12);
        assert!((t.apply_vector(p) - p).norm() < 1e-12);
    }

    #[test]
    fn translation_only_moves_points_but_not_vectors() {
        let t = PartTransform {
            translation: Vector3::new(1.0, 2.0, 3.0),
            orientation: UnitQuaternion::identity(),
        };
        assert!((t.apply_point(Vector3::zeros()) - Vector3::new(1.0, 2.0, 3.0)).norm() < 1e-12);
        // Vectors (directions) are translation-invariant.
        let v = Vector3::new(5.0, 0.0, 0.0);
        assert!((t.apply_vector(v) - v).norm() < 1e-12);
    }

    #[test]
    fn rotation_around_z_sends_x_axis_to_y_axis() {
        let t = PartTransform {
            translation: Vector3::zeros(),
            orientation: UnitQuaternion::from_axis_angle(&Vector3::z_axis(), PI / 2.0),
        };
        let v = Vector3::new(1.0, 0.0, 0.0);
        let r = t.apply_vector(v);
        assert!((r - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-9, "got {r:?}");
    }

    #[test]
    fn inverse_round_trips() {
        let t = PartTransform {
            translation: Vector3::new(3.0, -2.0, 5.0),
            orientation: UnitQuaternion::from_axis_angle(
                &nalgebra::Unit::new_normalize(Vector3::new(1.0, 1.0, 0.0)),
                0.7,
            ),
        };
        let inv = t.inverse();
        let p = Vector3::new(2.0, 4.0, -1.0);
        let r = inv.apply_point(t.apply_point(p));
        assert!((r - p).norm() < 1e-9, "round-trip drift: {r:?}");
    }

    #[test]
    fn part_with_unit_cube_has_unit_aabb_at_origin() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let p = Part::new(0, "cube", cube);
        let aabb = p.bounding_box();
        assert!((aabb.min - Vector3::new(0.0, 0.0, 0.0)).norm() < 0.01);
        assert!((aabb.max - Vector3::new(1.0, 1.0, 1.0)).norm() < 0.01);
    }

    #[test]
    fn part_aabb_picks_up_translation() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let mut p = Part::new(0, "cube", cube);
        p.transform.translation = Vector3::new(10.0, 0.0, 0.0);
        let aabb = p.bounding_box();
        assert!(aabb.min.x >= 9.9 && aabb.max.x <= 11.1);
    }
}
