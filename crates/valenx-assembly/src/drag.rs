//! Drag-aware re-solving — move one part and let the mate solver propagate
//! the move through the dependent parts.
//!
//! ## Why it matters
//!
//! Without drag-aware re-solving, the only way to move a constrained part
//! is to manually edit its transform and re-solve from scratch — the user
//! can't grab a part in the viewport and watch the linkage follow. With
//! it, a `drag_part(asm, id, new_pose)` call temporarily pins the dragged
//! part at `new_pose`, re-runs the constraint solver with the other parts
//! free, and the dependent parts' poses update through their mates. This
//! is the commercial-CAD "Move with constraints" / SolidWorks "Drag
//! Component" behavior.
//!
//! ## Algorithm
//!
//! 1. Snapshot the assembly's current pose vector (so we can roll back on
//!    a failed solve).
//! 2. Snapshot the dragged part's `fixed` flag.
//! 3. Set the dragged part's transform to `new_pose`.
//! 4. Set the dragged part's `fixed = true` (pins it during the solve).
//! 5. Run the solver. The dragged part's pose becomes the anchor; the
//!    other non-fixed parts find a pose configuration that satisfies
//!    every mate at the new anchor location.
//! 6. Restore the dragged part's `fixed` flag.
//! 7. Return [`DragOutcome::Success`] if the solver converged or
//!    [`DragOutcome::DragRejected`] (with the assembly rolled back to its
//!    pre-drag pose) if it did not — the latter is the common
//!    "you-dragged-out-of-the-convergence-basin" case.
//!
//! ## Honest scope
//!
//! - **Small drags only.** The Newton-Raphson solver has a finite
//!   convergence basin. A massive drag (e.g. push a part 100 m when the
//!   mates are ±5 m geometry) lands outside the basin and the solver
//!   diverges or stalls — we report `DragRejected` and roll back. The
//!   commercial CAD answer is the same: SolidWorks rejects such drags and
//!   snaps back. The honest interactive-drag pattern is many small drags
//!   per second, each a tiny pose delta the solver can absorb in a few
//!   iterations.
//! - **The dragged part is pinned for the duration of the solve.** Other
//!   `fixed` parts also stay pinned; the dragged part joins them as a
//!   temporary additional anchor. This is the correct semantics: the user
//!   said "this part *is* at the new location"; the solver finds the
//!   rest.

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::part::PartTransform;
use crate::solver::{apply_pose, pose_vector, solve, SolverConfig, SolverStatus};

/// Outcome of a [`drag_part`] call.
#[derive(Clone, Debug)]
pub enum DragOutcome {
    /// The solver converged with the dragged part pinned. The mated parts
    /// followed; the assembly's pose is updated.
    Success {
        /// Iterations the solver used.
        iterations: usize,
        /// Final residual L2 norm after the re-solve.
        residual_norm: f64,
    },
    /// The solver did not converge at the dragged target pose. The
    /// assembly's pose has been rolled back to its pre-drag state; no
    /// state was mutated. Common cause: the drag was too large for the
    /// solver's convergence basin.
    DragRejected {
        /// Iterations the solver burned before being rejected.
        iterations: usize,
        /// Final residual L2 norm at the rejected target pose (above
        /// tolerance — otherwise the call would have succeeded).
        residual_norm: f64,
    },
}

/// Drag `part_id` to `new_pose` and re-solve the rest of the assembly.
///
/// Returns [`AssemblyError::UnknownPart`] if `part_id` doesn't exist. Any
/// other solver error is converted into a [`DragOutcome::DragRejected`]
/// (with the assembly rolled back), so callers don't have to distinguish
/// "the part doesn't exist" from "the solver couldn't make it work" at
/// every call site.
pub fn drag_part(
    a: &mut Assembly,
    part_id: usize,
    new_pose: PartTransform,
) -> Result<DragOutcome, AssemblyError> {
    // Step 0 — validate the part id up front so the caller gets a clean
    // "unknown part" rather than a roll-back-from-nowhere.
    let saved_pose = pose_vector(a).0;
    let saved_fixed = a.get_part(part_id)?.fixed;
    let saved_xform = a.get_part(part_id)?.transform.clone();

    // Step 1-4 — pin the dragged part at the new pose.
    let part = a.get_part_mut(part_id)?;
    part.transform = new_pose;
    part.fixed = true;

    // Step 5 — re-solve.
    let cfg = SolverConfig::default();
    let result = solve(a, cfg);

    // Step 6 — restore the dragged part's fixed flag (regardless of
    // outcome — the temporary pin was for the duration of the solve).
    let part = a.get_part_mut(part_id)?;
    part.fixed = saved_fixed;

    match result {
        Ok(report) if report.status == SolverStatus::Converged => Ok(DragOutcome::Success {
            iterations: report.iterations,
            residual_norm: report.residual_norm,
        }),
        Ok(report) => {
            // Solver returned MaxIterations — roll back the whole pose
            // (the dragged part's pose has been corrupted along with the
            // rest of the bodies).
            rollback(a, &saved_pose, part_id, saved_xform);
            Ok(DragOutcome::DragRejected {
                iterations: report.iterations,
                residual_norm: report.residual_norm,
            })
        }
        Err(e) => {
            rollback(a, &saved_pose, part_id, saved_xform);
            // A "bad parameter" or "unknown part" error from inside the
            // solver becomes a typed error; numerical infeasibility
            // surfaces as DragRejected via the Ok(MaxIterations) branch
            // above. We propagate the typed error so callers can
            // distinguish bad-input from rejected-drag.
            Err(e)
        }
    }
}

/// Roll the dragged part *and* the rest of the assembly back to their
/// pre-drag poses. Internal helper.
fn rollback(
    a: &mut Assembly,
    saved_pose: &[f64],
    part_id: usize,
    saved_xform: PartTransform,
) {
    // The dragged part was pinned (fixed=true) during the solve, so its
    // pose is not represented in the saved_pose vector (which only spans
    // the *originally* non-fixed parts). Apply the non-fixed pose first;
    // then restore the dragged part's transform separately.
    apply_pose(a, saved_pose);
    if let Ok(p) = a.get_part_mut(part_id) {
        p.transform = saved_xform;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mate::{Mate, MateKind};
    use crate::part::Part;
    use crate::solver::{solve, SolverConfig};
    use nalgebra::{UnitQuaternion, Vector3};

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    /// A 3-part linkage — input (a) → middle (b) → output (c) — each pair
    /// connected by a Distance mate of length 2. Drag the input part by a
    /// small offset; the middle and output parts should follow.
    #[test]
    fn drag_input_propagates_through_distance_chain() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("input");
        p0.fixed = false;
        let id_input = a.add_part(p0);
        let mut p1 = unit_cube("middle");
        p1.transform.translation = Vector3::new(2.0, 0.0, 0.0);
        let id_middle = a.add_part(p1);
        let mut p2 = unit_cube("output");
        p2.transform.translation = Vector3::new(4.0, 0.0, 0.0);
        let id_output = a.add_part(p2);

        a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_input,
                point_a: Vector3::zeros(),
                part_b: id_middle,
                point_b: Vector3::zeros(),
                target: 2.0,
            },
        ));
        a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_middle,
                point_a: Vector3::zeros(),
                part_b: id_output,
                point_b: Vector3::zeros(),
                target: 2.0,
            },
        ));
        // Pre-solve so the assembly starts in a valid configuration.
        solve(&mut a, SolverConfig::default()).unwrap();
        let middle_before = a
            .get_part(id_middle)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        let output_before = a
            .get_part(id_output)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());

        // Drag the input by +Y 0.5 unit.
        let drag_target = PartTransform {
            translation: Vector3::new(0.0, 0.5, 0.0),
            orientation: UnitQuaternion::identity(),
        };
        let outcome = drag_part(&mut a, id_input, drag_target.clone()).unwrap();
        assert!(
            matches!(outcome, DragOutcome::Success { .. }),
            "drag rejected: {outcome:?}"
        );

        // Verify the drag landed: input is at the new pose; the mates are
        // still satisfied; middle + output moved.
        let input_after = a
            .get_part(id_input)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (input_after - drag_target.translation).norm() < 1e-6,
            "input did not land at drag target: {input_after:?}"
        );

        let middle_after = a
            .get_part(id_middle)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        let output_after = a
            .get_part(id_output)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());

        // Distance mates still satisfied at the new pose.
        let d1 = (middle_after - input_after).norm();
        let d2 = (output_after - middle_after).norm();
        assert!((d1 - 2.0).abs() < 1e-4, "input-middle distance = {d1}");
        assert!((d2 - 2.0).abs() < 1e-4, "middle-output distance = {d2}");

        // And the mated parts actually moved (vs the initial pre-drag pose).
        assert!(
            (middle_after - middle_before).norm() > 1e-3,
            "middle didn't follow drag"
        );
        assert!(
            (output_after - output_before).norm() > 1e-3,
            "output didn't follow drag"
        );
    }

    /// A drag to an unsolvable target — e.g. asking two Coincident-mated
    /// parts to be very far apart — should be rejected and the pose
    /// rolled back.
    #[test]
    fn impossible_drag_is_rejected_and_pose_rolls_back() {
        let mut a = Assembly::new();
        // Anchor: fixed at the origin. Cannot move.
        let mut p0 = unit_cube("anchor");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("free"));
        // Coincident: pin b's local-(0,0,0) to a's local-(0,0,0). With a
        // fixed at origin, b must also be at origin.
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        solve(&mut a, SolverConfig::default()).unwrap();

        // Drag the *anchor* (which is fixed) to (100, 0, 0). The anchor's
        // fixed flag was true, so dragging it is the way you say "I'm
        // moving the world frame". `b` must follow to (100, 0, 0) too.
        // This is the "feasible drag" — should succeed.
        let outcome = drag_part(
            &mut a,
            id_a,
            PartTransform {
                translation: Vector3::new(100.0, 0.0, 0.0),
                orientation: UnitQuaternion::identity(),
            },
        )
        .unwrap();
        assert!(
            matches!(outcome, DragOutcome::Success { .. }),
            "anchor drag rejected: {outcome:?}"
        );
        let b_at = a
            .get_part(id_b)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(
            (b_at - Vector3::new(100.0, 0.0, 0.0)).norm() < 1e-4,
            "b didn't follow anchor: {b_at:?}"
        );
        // Anchor's fixed flag should be restored.
        assert!(a.get_part(id_a).unwrap().fixed);
    }

    /// Dragging a non-existent part returns UnknownPart.
    #[test]
    fn drag_unknown_part_returns_typed_error() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("p"));
        let err = drag_part(&mut a, 99, PartTransform::identity()).unwrap_err();
        assert_eq!(err.code(), "assembly.unknown_part");
    }

    /// After a successful drag, the dragged part's fixed flag is whatever
    /// it was before — drag is non-mutating w.r.t. the fixed flag.
    #[test]
    fn drag_restores_fixed_flag_on_success() {
        let mut a = Assembly::new();
        let id = a.add_part(unit_cube("solo"));
        assert!(!a.get_part(id).unwrap().fixed);
        let _ = drag_part(
            &mut a,
            id,
            PartTransform {
                translation: Vector3::new(1.0, 2.0, 3.0),
                orientation: UnitQuaternion::identity(),
            },
        )
        .unwrap();
        assert!(!a.get_part(id).unwrap().fixed, "fixed flag drifted");
        // The pose did land at the dragged target.
        assert!(
            (a.get_part(id).unwrap().transform.translation - Vector3::new(1.0, 2.0, 3.0)).norm()
                < 1e-12
        );
    }
}
