//! Joint kinematics — translate a [`crate::Joint`]'s `parameter` into
//! a relative transform between its two parts.
//!
//! ## Semantics
//!
//! Each joint defines an axis or pivot in part_a's *local frame*.
//! Apply means: compute part_b's transform = part_a's transform ∘
//! (parameter-dependent local transform). This is the "child follows
//! parent" pattern from a scene graph — part_a is the parent, part_b
//! is the child.
//!
//! ## Per-kind table
//!
//! | Joint        | parameter   | what changes                            |
//! |--------------|-------------|-----------------------------------------|
//! | Fixed        | unused      | part_b ← part_a                         |
//! | Revolute     | angle (rad) | rotate part_b about axis on part_a      |
//! | Prismatic    | dist (m)    | translate part_b along axis_dir on part_a |
//! | Cylindrical  | angle (rad) | rotate part_b about axis (slide left to user) |
//! | Spherical    | unused      | snap part_b's pivot to part_a's anchor  |
//! | Planar       | unused      | project part_b's transform onto the plane |
//!
//! Spherical / Planar leave the orientation of part_b alone — they
//! only constrain *position*. v1 is preview-grade: the rotational
//! DOFs are not parameterized via the slider. Phase 6.5 will add
//! per-DOF parameter vectors.

use nalgebra::{Translation3, UnitQuaternion, Vector3};

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::joint::{Joint, JointKind};
use crate::part::PartTransform;

/// Apply one joint's parameter to its `part_b`'s transform. The
/// `part_a` side is left untouched.
///
/// **Semantics (v1 — absolute deterministic pose):** part_b's
/// transform is computed as a function of part_a's transform and the
/// joint's `parameter` alone. b's prior pose is **discarded** (with
/// exceptions for Spherical / Planar — see below). This lets the UI's
/// "joint slider" repeatedly drive motion idempotently: at `parameter
/// = 0` the joint is in its neutral pose, and any value can be set
/// without keeping a rest state.
///
/// Per-kind:
///
/// - **Fixed** — `b_pose = a_pose`. Parameter unused.
/// - **Revolute** — `b_pose = T(o_w) ∘ R(axis_w, θ) ∘ T(-o_w) ∘ a_pose`
///   where `o_w` and `axis_w` are the axis origin / direction in
///   world space (computed from a's pose). Parameter = rotation angle.
/// - **Prismatic** — `b_translation = a_translation + d · axis_w`,
///   `b_orientation = a_orientation`. Parameter = slide distance.
/// - **Cylindrical** — same as Revolute (the slide DOF is left at zero
///   for v1; Phase 6.5 will dual-parameterize).
/// - **Spherical** — `b_translation = world position of pivot in a`,
///   `b_orientation` is **preserved** (Spherical is free-rotation).
///   Parameter unused.
/// - **Planar** — `b_translation` is projected onto the plane on a,
///   `b_orientation` is **preserved**. Parameter unused.
pub fn apply_joint(a: &mut Assembly, j: &Joint) -> Result<(), AssemblyError> {
    if j.suppressed {
        return Ok(());
    }
    let (id_a, id_b) = j.kind.parts();
    let xform_a = a.get_part(id_a)?.transform.clone();
    let xform_b = a.get_part(id_b)?.transform.clone();
    let new_b: PartTransform = match &j.kind {
        JointKind::Fixed { .. } => xform_a,
        JointKind::Revolute {
            axis_origin,
            axis_dir,
            ..
        } => apply_revolute(&xform_a, *axis_origin, *axis_dir, j.parameter)?,
        JointKind::Prismatic { axis_dir, .. } => apply_prismatic(&xform_a, *axis_dir, j.parameter)?,
        JointKind::Cylindrical {
            axis_origin,
            axis_dir,
            ..
        } => {
            // Cylindrical = Revolute with the slide-along-axis component
            // left at zero (Phase 6.5 will dual-parameterize).
            apply_revolute(&xform_a, *axis_origin, *axis_dir, j.parameter)?
        }
        JointKind::Spherical { point, .. } => apply_spherical(&xform_a, *point, &xform_b)?,
        JointKind::Planar {
            plane_origin,
            plane_normal,
            ..
        } => apply_planar(&xform_a, *plane_origin, *plane_normal, &xform_b)?,
    };
    a.get_part_mut(id_b)?.transform = new_b;
    Ok(())
}

/// Apply every (un-suppressed) joint in the assembly. Joints are
/// applied in insertion order — for a tree-structured assembly this
/// is fine; closed kinematic loops need the constraint solver
/// instead.
pub fn apply_all_joints(a: &mut Assembly) -> Result<(), AssemblyError> {
    let joints = a.joints.clone();
    for j in joints.iter() {
        apply_joint(a, j)?;
    }
    Ok(())
}

/// Build part_b's transform for a Revolute joint at angle `theta`.
///
/// At θ=0 b coincides with a (`b_pose = a_pose`). For θ ≠ 0, b's pose
/// is a's pose with an extra rotation about the axis `(axis_origin,
/// axis_dir)` expressed in a's local frame.
///
/// Math: world-space axis = a.apply_vector(axis_dir); world-space
/// origin = a.apply_point(axis_origin). The rotation Q(axis_w, θ) is
/// composed about that pivot:
///   b_pose = T(o_w) ∘ Q ∘ T(-o_w) ∘ a_pose
fn apply_revolute(
    xform_a: &PartTransform,
    axis_origin: Vector3<f64>,
    axis_dir: Vector3<f64>,
    theta: f64,
) -> Result<PartTransform, AssemblyError> {
    let axis_dir_len = axis_dir.norm();
    if !axis_dir_len.is_finite() || axis_dir_len < 1e-12 {
        return Err(AssemblyError::BadParameter {
            name: "axis_dir",
            reason: "must be a non-zero finite vector".into(),
        });
    }
    let axis_unit = nalgebra::Unit::new_unchecked(axis_dir / axis_dir_len);
    let axis_world = xform_a.apply_vector(axis_unit.into_inner());
    let axis_world_unit = nalgebra::Unit::new_normalize(axis_world);
    let origin_world = xform_a.apply_point(axis_origin);

    // Compose: T(o_w) ∘ Q ∘ T(-o_w) ∘ a_pose. Walking through:
    //   x ↦ R_a · x + t_a                      (apply a_pose)
    //   x ↦ R_a · x + t_a - o_w                (T(-o_w))
    //   x ↦ Q · (R_a · x + t_a - o_w)          (Q)
    //   x ↦ Q · (R_a · x + t_a - o_w) + o_w    (T(o_w))
    //   = (Q · R_a) · x + Q · (t_a - o_w) + o_w
    let q_rot = UnitQuaternion::from_axis_angle(&axis_world_unit, theta);
    let new_rot = q_rot * xform_a.orientation;
    let new_t = q_rot * (xform_a.translation - origin_world) + origin_world;
    Ok(PartTransform {
        translation: new_t,
        orientation: new_rot,
    })
}

/// Build part_b's transform for a Prismatic joint at distance `d`
/// along `axis_dir` (expressed in part_a's local frame). At d=0,
/// b's translation = a's translation; b's orientation = a's
/// orientation.
fn apply_prismatic(
    xform_a: &PartTransform,
    axis_dir: Vector3<f64>,
    d: f64,
) -> Result<PartTransform, AssemblyError> {
    let axis_dir_len = axis_dir.norm();
    if !axis_dir_len.is_finite() || axis_dir_len < 1e-12 {
        return Err(AssemblyError::BadParameter {
            name: "axis_dir",
            reason: "must be a non-zero finite vector".into(),
        });
    }
    let axis_unit = axis_dir / axis_dir_len;
    let axis_world = xform_a.apply_vector(axis_unit);
    let offset = Translation3::from(axis_world * d);
    Ok(PartTransform {
        translation: offset.vector + xform_a.translation,
        orientation: xform_a.orientation,
    })
}

/// Spherical joint — snap part_b's translation to the pivot point on
/// part_a, leaving part_b's orientation untouched.
fn apply_spherical(
    xform_a: &PartTransform,
    point: Vector3<f64>,
    xform_b: &PartTransform,
) -> Result<PartTransform, AssemblyError> {
    let pivot_world = xform_a.apply_point(point);
    Ok(PartTransform {
        translation: pivot_world,
        orientation: xform_b.orientation,
    })
}

/// Planar joint — project part_b's translation onto the plane
/// (origin + normal expressed in part_a's local frame), keeping
/// part_b's orientation unchanged.
fn apply_planar(
    xform_a: &PartTransform,
    plane_origin: Vector3<f64>,
    plane_normal: Vector3<f64>,
    xform_b: &PartTransform,
) -> Result<PartTransform, AssemblyError> {
    let n_len = plane_normal.norm();
    if !n_len.is_finite() || n_len < 1e-12 {
        return Err(AssemblyError::BadParameter {
            name: "plane_normal",
            reason: "must be a non-zero finite vector".into(),
        });
    }
    let n_unit = plane_normal / n_len;
    let plane_origin_w = xform_a.apply_point(plane_origin);
    let n_world = xform_a.apply_vector(n_unit).normalize();

    // Project xform_b.translation onto the plane defined by
    // (plane_origin_w, n_world).
    let delta = xform_b.translation - plane_origin_w;
    let perpendicular = n_world * delta.dot(&n_world);
    let projected = xform_b.translation - perpendicular;
    Ok(PartTransform {
        translation: projected,
        orientation: xform_b.orientation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::Part;
    use std::f64::consts::PI;

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    /// Task 24 — Revolute around Z-axis. Part a is fixed at origin.
    /// Part b is connected by a Revolute joint. Setting parameter =
    /// π/2 rotates b by 90° about Z relative to a.
    #[test]
    fn revolute_rotates_b_about_z_by_parameter() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));

        let mut joint = Joint::new(
            0,
            JointKind::Revolute {
                part_a: id_a,
                part_b: id_b,
                axis_origin: Vector3::zeros(),
                axis_dir: Vector3::z(),
            },
        );
        joint.parameter = PI / 2.0;
        a.add_joint(joint);

        apply_all_joints(&mut a).unwrap();

        // A point originally at (1, 0, 0) on b's local frame should
        // now be at (0, 1, 0) in world space.
        let pb = a.get_part(id_b).unwrap();
        let p_local = Vector3::new(1.0, 0.0, 0.0);
        let p_world = pb.transform.apply_point(p_local);
        assert!(
            (p_world - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-9,
            "got {p_world:?}"
        );
    }

    #[test]
    fn revolute_about_offset_axis_keeps_pivot_invariant() {
        // Absolute-pose semantic: at θ=π, b's *pose* equals a's pose
        // rotated 180° about the pivot. The pivot point itself stays
        // put in world coordinates (rotation invariant), so b's
        // world-space copy of the pivot point sits where it always did.
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        let pivot = Vector3::new(1.0, 0.0, 0.0);
        let mut joint = Joint::new(
            0,
            JointKind::Revolute {
                part_a: id_a,
                part_b: id_b,
                axis_origin: pivot,
                axis_dir: Vector3::z(),
            },
        );
        joint.parameter = PI;
        a.add_joint(joint);

        apply_all_joints(&mut a).unwrap();

        // b's local copy of the pivot (= the pivot expressed in b's local
        // frame, same value because b is rigid) should map back to the
        // same world location — that's the defining property of rotation
        // about an axis. b's *origin*, however, swings to the opposite
        // side: from (0,0,0) → (2,0,0).
        let pb = a.get_part(id_b).unwrap();
        let pivot_world = pb.transform.apply_point(pivot);
        assert!(
            (pivot_world - pivot).norm() < 1e-9,
            "pivot moved: {pivot_world:?}"
        );
        let origin_world = pb.transform.apply_point(Vector3::zeros());
        assert!(
            (origin_world - Vector3::new(2.0, 0.0, 0.0)).norm() < 1e-9,
            "origin at {origin_world:?}"
        );
    }

    /// Task 25 — Prismatic + chain. Three parts in a slider mechanism:
    /// a (fixed at origin), b (slides along X relative to a), c
    /// (slides along Y relative to b). Setting parameters drives them.
    #[test]
    fn prismatic_plus_cylindrical_chain() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        let id_c = a.add_part(unit_cube("c"));

        let mut j1 = Joint::new(
            0,
            JointKind::Prismatic {
                part_a: id_a,
                part_b: id_b,
                axis_dir: Vector3::x(),
            },
        );
        j1.parameter = 3.0;
        a.add_joint(j1);

        let mut j2 = Joint::new(
            0,
            JointKind::Prismatic {
                part_a: id_b,
                part_b: id_c,
                axis_dir: Vector3::y(),
            },
        );
        j2.parameter = 2.0;
        a.add_joint(j2);

        apply_all_joints(&mut a).unwrap();

        let pb = a.get_part(id_b).unwrap();
        let pc = a.get_part(id_c).unwrap();
        let b_origin = pb.transform.apply_point(Vector3::zeros());
        let c_origin = pc.transform.apply_point(Vector3::zeros());
        assert!(
            (b_origin - Vector3::new(3.0, 0.0, 0.0)).norm() < 1e-9,
            "b at {b_origin:?}"
        );
        assert!(
            (c_origin - Vector3::new(3.0, 2.0, 0.0)).norm() < 1e-9,
            "c at {c_origin:?}"
        );
    }

    #[test]
    fn fixed_joint_clones_a_transform_to_b() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.transform.translation = Vector3::new(5.0, 7.0, 9.0);
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        a.add_joint(Joint::new(
            0,
            JointKind::Fixed {
                part_a: id_a,
                part_b: id_b,
            },
        ));
        apply_all_joints(&mut a).unwrap();
        let pb = a.get_part(id_b).unwrap();
        assert!((pb.transform.translation - Vector3::new(5.0, 7.0, 9.0)).norm() < 1e-12);
    }

    #[test]
    fn spherical_pins_b_origin_to_pivot_preserves_orientation() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(10.0, 10.0, 10.0);
        // Give b a non-identity orientation to verify it's preserved.
        p1.transform.orientation = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.7);
        let id_b = a.add_part(p1);
        a.add_joint(Joint::new(
            0,
            JointKind::Spherical {
                part_a: id_a,
                part_b: id_b,
                point: Vector3::new(2.0, 0.0, 0.0),
            },
        ));
        apply_all_joints(&mut a).unwrap();
        let pb = a.get_part(id_b).unwrap();
        // b's translation should now be (2, 0, 0) in world (a's local
        // frame is the world frame since a is at the identity).
        assert!((pb.transform.translation - Vector3::new(2.0, 0.0, 0.0)).norm() < 1e-12);
        // Spherical preserves b's orientation (free rotation).
        let expected_q = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.7);
        let dot = pb
            .transform
            .orientation
            .coords
            .dot(&expected_q.coords)
            .abs();
        assert!((dot - 1.0).abs() < 1e-9, "orientation drifted");
    }

    #[test]
    fn planar_projects_b_onto_plane() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("b");
        // b at (1, 1, 5) — should project onto z=0 plane to (1, 1, 0).
        p1.transform.translation = Vector3::new(1.0, 1.0, 5.0);
        let id_b = a.add_part(p1);
        a.add_joint(Joint::new(
            0,
            JointKind::Planar {
                part_a: id_a,
                part_b: id_b,
                plane_origin: Vector3::zeros(),
                plane_normal: Vector3::z(),
            },
        ));
        apply_all_joints(&mut a).unwrap();
        let pb = a.get_part(id_b).unwrap();
        assert!((pb.transform.translation - Vector3::new(1.0, 1.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn revolute_with_zero_axis_errors() {
        let mut a = Assembly::new();
        let id_a = a.add_part(unit_cube("a"));
        let id_b = a.add_part(unit_cube("b"));
        let mut j = Joint::new(
            0,
            JointKind::Revolute {
                part_a: id_a,
                part_b: id_b,
                axis_origin: Vector3::zeros(),
                axis_dir: Vector3::zeros(),
            },
        );
        j.parameter = 0.5;
        let err = apply_joint(&mut a, &j).unwrap_err();
        assert_eq!(err.code(), "assembly.bad_parameter");
    }

    #[test]
    fn suppressed_joint_is_a_no_op() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        let mut j = Joint::new(
            0,
            JointKind::Prismatic {
                part_a: id_a,
                part_b: id_b,
                axis_dir: Vector3::x(),
            },
        );
        j.parameter = 100.0;
        j.suppressed = true;
        a.add_joint(j);
        apply_all_joints(&mut a).unwrap();
        // b's transform should still be identity since the joint was suppressed.
        let pb = a.get_part(id_b).unwrap();
        assert!(pb.transform.translation.norm() < 1e-12);
    }
}
