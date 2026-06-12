//! Assembly-constraint-aware gizmo — drag one part, re-solve the
//! mates / joints across the whole assembly.
//!
//! ## What OCCT (and a real CAD assembly UX) does
//!
//! A bare manipulator ([`TranslationGizmo`](crate::TranslationGizmo) /
//! [`RotationGizmo`](crate::RotationGizmo)) moves *one* object. In a
//! CAD **assembly**, the parts are tied together by mate / joint
//! constraints — a coincident mate, a distance mate, a revolute joint.
//! Grabbing one part and dragging it must therefore not just move that
//! part: the mated parts have to *follow*, re-solving the constraint
//! network so the assembly stays consistent. SolidWorks calls this
//! "move with triad + dynamic mate solve"; the constraint solver runs
//! live inside the drag loop.
//!
//! ## What this module is
//!
//! The cross-crate integration of Valenx's two existing systems:
//!
//! - [`valenx_assembly`] (Phase 6) — parts, mates, joints, and the
//!   Newton-Raphson [`valenx_assembly::solver`] 3-D constraint solver.
//! - The gizmos in this crate — [`crate::TranslationGizmo`] /
//!   [`crate::RotationGizmo`] — which turn a cursor drag into a
//!   rigid-body delta.
//!
//! [`apply_constraint_drag`] is the bridge: given an assembly, the id
//! of the dragged part, and the gizmo's drag delta, it (1) applies the
//! drag to the dragged part, (2) **pins** the dragged part so the
//! re-solve treats it as the anchor (the drag "wins"), (3) runs the
//! assembly constraint solver so every mated part moves to keep its
//! mates satisfied, and (4) restores the original `fixed` flags. The
//! result is the **constraint-propagating drag** — drag one part, the
//! whole assembly follows.
//!
//! [`AssemblyDragSession`] wraps that into a begin/update/end drag
//! state machine so an app can call `update` once per cursor move.
//!
//! ## Honest scope
//!
//! - The solve is a **full re-solve per drag update** — every call to
//!   [`apply_constraint_drag`] (or [`AssemblyDragSession::update`])
//!   runs the Newton-Raphson solver from the current pose. For the
//!   typical handful-of-parts assembly that is fast and correct; a
//!   large assembly would want an incremental / warm-started solve
//!   (the solver already warm-starts from the current pose, so a small
//!   drag step converges in few iterations, but the cost still scales
//!   with the constraint count).
//! - The **live in-app drag-loop hookup** — wiring this to the egui +
//!   wgpu viewport's per-frame cursor events and the gizmo pick — is
//!   the app-layer follow-on. This module is the crate-layer
//!   deliverable: the constraint-propagating drag as a pure,
//!   tested function the app calls.
//! - Mates drive the re-solve (the constraint solver consumes mates).
//!   Joints carry a *parameter* and are applied by
//!   [`valenx_assembly::kinematics`], not by the constraint solver —
//!   so a joint-only assembly re-solves to a no-op here; pair a joint
//!   with the matching mate, or drive the joint parameter via the
//!   kinematics applier, for joint-constrained motion.

use nalgebra::{UnitQuaternion, Vector3};

use valenx_assembly::solver::{solve, SolverConfig, SolverReport};
use valenx_assembly::{Assembly, AssemblyError};

use crate::error::OcctVizError;

/// A rigid-body drag delta produced by a gizmo — a translation plus a
/// rotation, to be composed onto the dragged part's transform.
///
/// A pure translation drag (from [`crate::TranslationGizmo`]) sets
/// `rotation` to the identity; a pure rotation drag (from
/// [`crate::RotationGizmo`]) sets `translation` to zero. A 6-DOF drag
/// sets both.
#[derive(Clone, Copy, Debug)]
pub struct DragDelta {
    /// World-space translation to add to the dragged part.
    pub translation: Vector3<f64>,
    /// Rotation to compose onto the dragged part's orientation.
    /// Identity for a translation-only drag.
    pub rotation: UnitQuaternion<f64>,
}

impl DragDelta {
    /// A translation-only drag delta (the [`crate::TranslationGizmo`]
    /// case).
    pub fn translation(delta: Vector3<f64>) -> DragDelta {
        DragDelta {
            translation: delta,
            rotation: UnitQuaternion::identity(),
        }
    }

    /// A rotation-only drag delta about `axis` (need not be unit) by
    /// `angle` radians — the [`crate::RotationGizmo`] case. A
    /// near-zero axis yields the identity rotation.
    pub fn rotation(axis: Vector3<f64>, angle: f64) -> DragDelta {
        let rot = match nalgebra::Unit::try_new(axis, 1e-12) {
            Some(unit) => UnitQuaternion::from_axis_angle(&unit, angle),
            None => UnitQuaternion::identity(),
        };
        DragDelta {
            translation: Vector3::zeros(),
            rotation: rot,
        }
    }

    /// A combined 6-DOF drag delta.
    pub fn rigid(translation: Vector3<f64>, rotation: UnitQuaternion<f64>) -> DragDelta {
        DragDelta {
            translation,
            rotation,
        }
    }

    /// The zero (no-op) drag delta.
    pub fn identity() -> DragDelta {
        DragDelta {
            translation: Vector3::zeros(),
            rotation: UnitQuaternion::identity(),
        }
    }
}

/// The outcome of a constraint-propagating drag.
#[derive(Clone, Debug)]
pub struct ConstraintDragResult {
    /// The constraint solver's report from the post-drag re-solve —
    /// status (converged / max-iterations), iteration count, and the
    /// final residual norm. A converged report means the assembly's
    /// mates are satisfied after the drag.
    pub solver_report: SolverReport,
    /// The id of the part that was dragged (echoed for convenience).
    pub dragged_part: usize,
}

impl ConstraintDragResult {
    /// True when the post-drag constraint solve converged — the
    /// assembly's mates are satisfied to the solver tolerance.
    pub fn converged(&self) -> bool {
        use valenx_assembly::solver::SolverStatus;
        matches!(self.solver_report.status, SolverStatus::Converged)
    }
}

/// Apply a gizmo drag to one part of an assembly and re-solve the
/// constraint network so the mated parts follow.
///
/// This is the constraint-propagating drag. It:
///
/// 1. Composes `delta` onto `dragged_part`'s transform (translation
///    added, rotation pre-multiplied onto the orientation).
/// 2. Temporarily marks `dragged_part` as `fixed` so the constraint
///    solver treats it as the anchor — the user's drag is authoritative
///    and the *other* parts move to satisfy their mates.
/// 3. Runs the [`valenx_assembly::solver`] Newton-Raphson solver,
///    which adjusts every non-fixed part's pose to drive the mate
///    residuals to zero.
/// 4. Restores every part's original `fixed` flag.
///
/// The solver warm-starts from the current pose, so a small drag step
/// converges quickly.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `dragged_part` is not a part id in
///   `assembly`, or `delta` carries a non-finite component.
/// - [`OcctVizError::BadInput`] (wrapping the solver's
///   [`AssemblyError`]) if the constraint solve fails structurally.
pub fn apply_constraint_drag(
    assembly: &mut Assembly,
    dragged_part: usize,
    delta: &DragDelta,
    config: SolverConfig,
) -> Result<ConstraintDragResult, OcctVizError> {
    // Validate the drag delta — a NaN must not silently corrupt poses.
    if !delta.translation.iter().all(|c| c.is_finite())
        || !delta.rotation.coords.iter().all(|c| c.is_finite())
    {
        return Err(OcctVizError::bad_input(
            "delta",
            "drag delta has a non-finite component",
        ));
    }

    // Snapshot every part's `fixed` flag so we can restore it after
    // the re-solve, and confirm the dragged part exists.
    let mut original_fixed: Vec<(usize, bool)> =
        assembly.parts.iter().map(|p| (p.id, p.fixed)).collect();
    if !original_fixed.iter().any(|&(id, _)| id == dragged_part) {
        return Err(OcctVizError::bad_input(
            "dragged_part",
            format!("no part with id {dragged_part} in the assembly"),
        ));
    }

    // (1) Apply the drag to the dragged part. Translation adds;
    // rotation pre-multiplies (a world-frame rotation of the part).
    {
        let part = assembly
            .get_part_mut(dragged_part)
            .map_err(assembly_err_to_bad_input)?;
        part.transform.translation += delta.translation;
        part.transform.orientation = delta.rotation * part.transform.orientation;
    }

    // (2) Pin the dragged part so the re-solve anchors on it.
    {
        let part = assembly
            .get_part_mut(dragged_part)
            .map_err(assembly_err_to_bad_input)?;
        part.fixed = true;
    }

    // (3) Re-solve the constraint network — mated parts move to keep
    // their mates satisfied around the dragged anchor.
    let solve_result = solve(assembly, config);

    // (4) Restore every original `fixed` flag, whatever the solve did.
    for part in &mut assembly.parts {
        if let Some(&(_, was_fixed)) = original_fixed.iter().find(|&&(id, _)| id == part.id) {
            part.fixed = was_fixed;
        }
    }
    original_fixed.clear();

    let solver_report = solve_result.map_err(assembly_err_to_bad_input)?;
    Ok(ConstraintDragResult {
        solver_report,
        dragged_part,
    })
}

/// A constraint-aware drag session — the begin / update / end state
/// machine an app drives from cursor events.
///
/// The session records the dragged part and the assembly pose at the
/// drag start. Each [`AssemblyDragSession::update`] receives the
/// **total** drag delta from the drag origin (exactly what
/// [`crate::TranslationGizmo::update_drag`] returns), restores the
/// assembly to the drag-start pose, re-applies the cumulative delta,
/// and re-solves — so the drag is always evaluated against the
/// original pose and dragging back to the origin returns the assembly
/// to where it started.
#[derive(Clone, Debug)]
pub struct AssemblyDragSession {
    /// The id of the part being dragged.
    dragged_part: usize,
    /// The assembly pose (per-part transform) captured at drag start.
    start_poses: Vec<(usize, PoseSnapshot)>,
    /// Solver configuration used for every re-solve.
    config: SolverConfig,
}

/// A captured rigid-body pose — the snapshot the drag session restores
/// to before re-applying the cumulative delta.
#[derive(Clone, Copy, Debug)]
struct PoseSnapshot {
    translation: Vector3<f64>,
    orientation: UnitQuaternion<f64>,
}

impl AssemblyDragSession {
    /// Begin a constraint-aware drag of `dragged_part` in `assembly`.
    ///
    /// Captures the current assembly pose. The session does not modify
    /// the assembly until the first [`AssemblyDragSession::update`].
    ///
    /// # Errors
    ///
    /// [`OcctVizError::BadInput`] if `dragged_part` is not a part id in
    /// `assembly`.
    pub fn begin(
        assembly: &Assembly,
        dragged_part: usize,
        config: SolverConfig,
    ) -> Result<AssemblyDragSession, OcctVizError> {
        if !assembly.parts.iter().any(|p| p.id == dragged_part) {
            return Err(OcctVizError::bad_input(
                "dragged_part",
                format!("no part with id {dragged_part} in the assembly"),
            ));
        }
        let start_poses = assembly
            .parts
            .iter()
            .map(|p| {
                (
                    p.id,
                    PoseSnapshot {
                        translation: p.transform.translation,
                        orientation: p.transform.orientation,
                    },
                )
            })
            .collect();
        Ok(AssemblyDragSession {
            dragged_part,
            start_poses,
            config,
        })
    }

    /// Update the drag with the **total** delta from the drag start.
    ///
    /// Restores every part to its drag-start pose, applies `total_delta`
    /// to the dragged part, and re-solves the constraint network. The
    /// returned [`ConstraintDragResult`] reports the solve outcome.
    ///
    /// Passing a delta back toward zero walks the assembly back toward
    /// its drag-start configuration.
    ///
    /// # Errors
    ///
    /// As [`apply_constraint_drag`].
    pub fn update(
        &self,
        assembly: &mut Assembly,
        total_delta: &DragDelta,
    ) -> Result<ConstraintDragResult, OcctVizError> {
        // Restore the drag-start pose so the cumulative delta is always
        // applied to the original configuration, not compounded.
        for part in &mut assembly.parts {
            if let Some(&(_, snap)) = self.start_poses.iter().find(|&&(id, _)| id == part.id) {
                part.transform.translation = snap.translation;
                part.transform.orientation = snap.orientation;
            }
        }
        apply_constraint_drag(assembly, self.dragged_part, total_delta, self.config)
    }

    /// End the drag. The assembly keeps its current (last-update'd)
    /// pose — ending is purely a session-lifetime marker, so this just
    /// consumes the session.
    ///
    /// Returns the id of the part that was dragged.
    pub fn end(self) -> usize {
        self.dragged_part
    }

    /// The id of the part this session is dragging.
    pub fn dragged_part(&self) -> usize {
        self.dragged_part
    }
}

/// Map a [`valenx_assembly::AssemblyError`] onto an
/// [`OcctVizError::BadInput`] — the assembly errors (unknown part /
/// mate / joint) are all caller-input problems from this crate's
/// vantage point.
fn assembly_err_to_bad_input(e: AssemblyError) -> OcctVizError {
    OcctVizError::bad_input("assembly", e.to_string())
}

/// Construct a constraint-aware drag session — the entry point an app
/// calls when a gizmo drag starts on an assembly part.
///
/// Thin wrapper over [`AssemblyDragSession::begin`] with the default
/// solver configuration; use [`AssemblyDragSession::begin`] directly to
/// pass a tuned [`SolverConfig`].
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if `dragged_part` is not a part id in
/// `assembly`.
pub fn transformation_assembly_gizmo(
    assembly: &Assembly,
    dragged_part: usize,
) -> Result<AssemblyDragSession, OcctVizError> {
    AssemblyDragSession::begin(assembly, dragged_part, SolverConfig::default())
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_assembly::{Mate, MateKind, Part};

    /// A unit-cube part.
    fn cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    /// Build the canonical two-part test assembly: a fixed anchor cube
    /// at the origin and a moving cube, joined by a Coincident mate
    /// pinning the moving cube's local `point_b` to the anchor's
    /// `point_a`. Returns `(assembly, anchor_id, moving_id)`.
    fn coincident_pair() -> (Assembly, usize, usize) {
        let mut a = Assembly::new();
        let mut anchor = cube("anchor");
        anchor.fixed = true;
        let anchor_id = a.add_part(anchor);
        let moving = cube("moving");
        let moving_id = a.add_part(moving);
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: anchor_id,
                point_a: Vector3::zeros(),
                part_b: moving_id,
                point_b: Vector3::zeros(),
            },
        ));
        (a, anchor_id, moving_id)
    }

    #[test]
    fn drag_delta_constructors_are_correct() {
        let t = DragDelta::translation(Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(t.translation, Vector3::new(1.0, 2.0, 3.0));
        assert!(
            t.rotation.angle().abs() < 1e-12,
            "translation drag has no rotation"
        );

        let r = DragDelta::rotation(Vector3::z(), std::f64::consts::FRAC_PI_2);
        assert!(
            r.translation.norm() < 1e-12,
            "rotation drag has no translation"
        );
        assert!((r.rotation.angle() - std::f64::consts::FRAC_PI_2).abs() < 1e-9);

        // A zero rotation axis degrades to the identity rotation.
        let degenerate = DragDelta::rotation(Vector3::zeros(), 1.0);
        assert!(degenerate.rotation.angle().abs() < 1e-12);

        let id = DragDelta::identity();
        assert!(id.translation.norm() < 1e-12);
        assert!(id.rotation.angle().abs() < 1e-12);
    }

    #[test]
    fn rejects_an_unknown_dragged_part() {
        let (mut a, _anchor, _moving) = coincident_pair();
        let err = apply_constraint_drag(
            &mut a,
            999,
            &DragDelta::translation(Vector3::x()),
            SolverConfig::default(),
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_a_non_finite_drag_delta() {
        let (mut a, _anchor, moving) = coincident_pair();
        let bad = DragDelta::translation(Vector3::new(f64::NAN, 0.0, 0.0));
        let err = apply_constraint_drag(&mut a, moving, &bad, SolverConfig::default()).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn dragging_the_anchor_pulls_the_mated_part_along() {
        // The headline test of constraint propagation. The two cubes
        // are pinned together at their origins. Drag the *anchor* part
        // by (5, 0, 0); the constraint solver must move the mated
        // part so its mate stays satisfied — the mated part follows.
        let (mut a, anchor, moving) = coincident_pair();
        // The moving part starts at the origin (mate already satisfied).
        let before = a
            .get_part(moving)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(before.norm() < 1e-9, "mate starts satisfied");

        let result = apply_constraint_drag(
            &mut a,
            anchor,
            &DragDelta::translation(Vector3::new(5.0, 0.0, 0.0)),
            SolverConfig::default(),
        )
        .unwrap();
        assert!(
            result.converged(),
            "the re-solve should converge: {result:?}"
        );

        // The anchor moved to x = 5.
        let anchor_pos = a
            .get_part(anchor)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (anchor_pos.x - 5.0).abs() < 1e-6,
            "anchor at {anchor_pos:?}"
        );

        // The mated part FOLLOWED — its mate point coincides with the
        // anchor's mate point again (both now at x ≈ 5).
        let moving_pos = a
            .get_part(moving)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (moving_pos - anchor_pos).norm() < 1e-5,
            "the mated part must follow the drag: anchor {anchor_pos:?}, moving {moving_pos:?}"
        );
    }

    #[test]
    fn the_original_fixed_flags_are_restored_after_the_drag() {
        // The anchor was `fixed` before the drag; the moving part was
        // not. After `apply_constraint_drag` both flags must be exactly
        // as they started — the temporary pin is internal only.
        let (mut a, anchor, moving) = coincident_pair();
        assert!(a.get_part(anchor).unwrap().fixed, "anchor starts fixed");
        assert!(!a.get_part(moving).unwrap().fixed, "moving starts free");

        apply_constraint_drag(
            &mut a,
            moving,
            &DragDelta::translation(Vector3::new(1.0, 0.0, 0.0)),
            SolverConfig::default(),
        )
        .unwrap();

        assert!(
            a.get_part(anchor).unwrap().fixed,
            "anchor's fixed flag must be restored"
        );
        assert!(
            !a.get_part(moving).unwrap().fixed,
            "moving's fixed flag must be restored (the pin was temporary)"
        );
    }

    #[test]
    fn dragging_a_distance_mated_part_keeps_the_distance() {
        // Two cubes joined by a Distance mate of 5 units. `apply_
        // constraint_drag` pins whatever part is dragged (the drag is
        // authoritative) and re-solves the *other* parts — so to
        // exercise distance-mate propagation the dragged part must be
        // the one with a free mated partner. Drag the anchor; the free
        // moving part re-solves to stay 5 units from the dragged anchor.
        let mut a = Assembly::new();
        let mut anchor = cube("anchor");
        anchor.fixed = true;
        let anchor_id = a.add_part(anchor);
        let mut moving = cube("moving");
        moving.transform.translation = Vector3::new(5.0, 0.0, 0.0);
        let moving_id = a.add_part(moving);
        a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: anchor_id,
                point_a: Vector3::zeros(),
                part_b: moving_id,
                point_b: Vector3::zeros(),
                target: 5.0,
            },
        ));

        // Drag the anchor out by (7, 0, 0) to x = 7 — the Distance
        // mate is now violated (the moving part is still at x = 5, only
        // 2 units away); the re-solve must move the free part back onto
        // the 5-unit sphere around the dragged anchor.
        let result = apply_constraint_drag(
            &mut a,
            anchor_id,
            &DragDelta::translation(Vector3::new(7.0, 0.0, 0.0)),
            SolverConfig::default(),
        )
        .unwrap();
        assert!(result.converged(), "{result:?}");
        let anchor_pos = a
            .get_part(anchor_id)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        let moving_pos = a
            .get_part(moving_id)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            ((moving_pos - anchor_pos).norm() - 5.0).abs() < 1e-5,
            "the Distance mate must hold after the drag: separation = {}",
            (moving_pos - anchor_pos).norm()
        );
    }

    #[test]
    fn drag_session_begin_rejects_unknown_part() {
        let (a, _anchor, _moving) = coincident_pair();
        let err = AssemblyDragSession::begin(&a, 777, SolverConfig::default()).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn drag_session_update_then_back_to_origin_restores_the_assembly() {
        // A drag session evaluates the *total* delta against the
        // drag-start pose. Dragging out, then back to a zero delta,
        // must return every part to its start pose.
        let (mut a, anchor, moving) = coincident_pair();
        let session = AssemblyDragSession::begin(&a, anchor, SolverConfig::default()).unwrap();
        let moving_start = a
            .get_part(moving)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());

        // Drag the anchor out by (4, 1, 0).
        session
            .update(&mut a, &DragDelta::translation(Vector3::new(4.0, 1.0, 0.0)))
            .unwrap();
        let moving_dragged = a
            .get_part(moving)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (moving_dragged - moving_start).norm() > 1.0,
            "the mated part should have moved during the drag"
        );

        // Drag back to a zero total delta — the assembly returns to
        // its drag-start configuration.
        session.update(&mut a, &DragDelta::identity()).unwrap();
        let moving_back = a
            .get_part(moving)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (moving_back - moving_start).norm() < 1e-5,
            "dragging back to zero must restore the assembly: start {moving_start:?}, back {moving_back:?}"
        );
        assert_eq!(session.dragged_part(), anchor);
    }

    #[test]
    fn a_rotation_drag_propagates_through_a_coincident_mate() {
        // Rotate the anchor part about Z. The Coincident mate pins the
        // two parts at their *local* origins, so the mated part must
        // re-solve to keep its mate point on the anchor's mate point.
        // (Both mate points are the local origin, which a rotation
        // about the part origin leaves at the world origin — so the
        // mated part stays satisfied; the test confirms the rotation
        // drag path runs and the solve converges.)
        let (mut a, anchor, _moving) = coincident_pair();
        let result = apply_constraint_drag(
            &mut a,
            anchor,
            &DragDelta::rotation(Vector3::z(), std::f64::consts::FRAC_PI_2),
            SolverConfig::default(),
        )
        .unwrap();
        assert!(result.converged(), "rotation-drag re-solve: {result:?}");
        // The anchor genuinely rotated 90° about Z.
        let rotated_x = a
            .get_part(anchor)
            .unwrap()
            .transform
            .apply_vector(Vector3::x());
        assert!(
            (rotated_x - Vector3::y()).norm() < 1e-9,
            "anchor should have rotated +X onto +Y, got {rotated_x:?}"
        );
    }

    #[test]
    fn transformation_assembly_gizmo_builds_a_session() {
        // The thin entry point builds a session with the default
        // solver config.
        let (a, anchor, _moving) = coincident_pair();
        let session = transformation_assembly_gizmo(&a, anchor).unwrap();
        assert_eq!(session.dragged_part(), anchor);
        // end() consumes the session and returns the dragged id.
        assert_eq!(session.end(), anchor);
    }

    #[test]
    fn a_three_part_chain_propagates_the_drag_down_the_chain() {
        // Three cubes a—b—c joined by two Coincident mates. Anchor `a`
        // fixed. Drag `a`; both `b` and `c` must follow so the chain
        // stays assembled.
        let mut asm = Assembly::new();
        let mut a = cube("a");
        a.fixed = true;
        let id_a = asm.add_part(a);
        let id_b = asm.add_part(cube("b"));
        let id_c = asm.add_part(cube("c"));
        asm.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        asm.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_b,
                point_a: Vector3::new(1.0, 0.0, 0.0),
                part_b: id_c,
                point_b: Vector3::zeros(),
            },
        ));

        let result = apply_constraint_drag(
            &mut asm,
            id_a,
            &DragDelta::translation(Vector3::new(0.0, 3.0, 0.0)),
            SolverConfig::default(),
        )
        .unwrap();
        assert!(result.converged(), "{result:?}");

        // a—b coincident: b's origin meets a's origin (now at y=3).
        let pa = asm
            .get_part(id_a)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        let pb = asm
            .get_part(id_b)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!((pa - pb).norm() < 1e-5, "b must follow a: {pa:?} vs {pb:?}");
        // b—c coincident: c's origin meets b's local (1,0,0).
        let pb1 = asm
            .get_part(id_b)
            .unwrap()
            .transform
            .apply_point(Vector3::new(1.0, 0.0, 0.0));
        let pc = asm
            .get_part(id_c)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (pb1 - pc).norm() < 1e-5,
            "c must follow b: {pb1:?} vs {pc:?}"
        );
    }
}
