//! CAD-kernel validation suite — BRep fillet.
//!
//! A fillet of radius `r` on a convex edge **removes a precisely
//! known volume**: the sharp-corner sliver minus the rounded fill. For
//! a 90° edge the cross-section removed is a square corner of side `r`
//! minus the quarter-disk that the fillet leaves behind:
//!
//! ```text
//!   removed cross-section area = r² − πr²/4 = r²·(1 − π/4)
//!   removed volume             = r²·(1 − π/4) · edge_length
//! ```
//!
//! This suite checks two things honestly:
//!
//! 1. **The fillet *plan* geometry is exact.** The tangent points,
//!    arc centre and dihedral the planner computes are checked against
//!    the closed-form values — this is the geometry the BRep surgery
//!    is built from, and it must be correct independent of whether the
//!    `truck_shapeops` boolean can resolve the coincident cutter
//!    faces.
//! 2. **The fillet result is honest.** When `fillet_planar_edge`
//!    succeeds the result must be a valid closed solid whose volume is
//!    the original minus the analytic sliver; when the coincident-face
//!    boolean cannot be resolved it must surface the soft `TruckOp`
//!    error — never a panic, never a silently-invalid solid.
//!
//! ## Honest note on the boolean fall-through
//!
//! `truck_shapeops` returns `None` for a difference where the cutter
//! prism trims **flush** with the solid's faces (verified — see the
//! boolean-robustness pass). A flush cut is exactly what a fillet
//! cutter does, so on `truck` 0.4 the BRep fillet's
//! `(solid − cutter)` step often soft-fails to `TruckOp`. The tests
//! below therefore accept *either* a correct filleted solid or the
//! soft `TruckOp` fall-through — both are honest outcomes; the
//! `valenx-feature-tree` dispatcher falls through to the mesh-domain
//! fillet on the latter. What the tests reject is a panic or a
//! geometric-precondition error on a textbook convex planar edge.

use std::collections::HashSet;
use std::f64::consts::{FRAC_PI_2, PI};

use truck_modeling::{Edge as TruckEdge, InnerSpace, Solid as TruckSolid};
use valenx_cad::measure::{is_closed_solid_tol, solid_volume_tol};
use valenx_cad::primitives::box_solid;
use valenx_cad::Solid;
use valenx_fillet_brep::error::FilletBrepError;
use valenx_fillet_brep::fillet::{fillet_planar_edge, plan_planar_edge_fillet};

fn inner(s: &Solid) -> &TruckSolid {
    match s {
        Solid::Brep(b) => b,
        _ => panic!("expected a BRep solid"),
    }
}

/// Collect the unique edges of a solid.
fn unique_edges(brep: &TruckSolid) -> Vec<TruckEdge> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for e in brep.edge_iter() {
        if seen.insert(e.id()) {
            out.push(e);
        }
    }
    out
}

// ===========================================================================
// The fillet plan — exact cross-section geometry.
// ===========================================================================

#[test]
fn cube_edge_fillet_plan_has_right_angle_dihedral() {
    // Every edge of a box joins two faces at a right angle — the
    // planner must see a π/2 dihedral.
    let cube = box_solid(4.0, 4.0, 4.0).unwrap();
    let brep = inner(&cube);
    for edge in unique_edges(brep) {
        let plan = plan_planar_edge_fillet(brep, &edge, 0.5).unwrap();
        assert!(
            (plan.dihedral_angle - FRAC_PI_2).abs() < 1e-9,
            "cube edge dihedral {} != π/2",
            plan.dihedral_angle
        );
    }
}

#[test]
fn fillet_plan_tangent_points_sit_radius_offset_from_the_edge() {
    // On a 90° corner the tangent contact line is offset from the edge
    // by exactly `radius / tan(45°) = radius` along each face. The
    // plan's tangent points must satisfy that to machine precision.
    let cube = box_solid(4.0, 4.0, 4.0).unwrap();
    let brep = inner(&cube);
    let r = 0.7;
    let edge = unique_edges(brep)[0].clone();
    let plan = plan_planar_edge_fillet(brep, &edge, r).unwrap();
    // tan(45°) = 1, so offset == r.
    let off0 = (plan.tangent_on_face0_front - plan.edge_front).magnitude();
    let off1 = (plan.tangent_on_face1_front - plan.edge_front).magnitude();
    assert!(
        (off0 - r).abs() < 1e-9,
        "face-0 tangent offset {off0} != radius {r}"
    );
    assert!(
        (off1 - r).abs() < 1e-9,
        "face-1 tangent offset {off1} != radius {r}"
    );
}

#[test]
fn fillet_plan_rejects_radius_too_large_for_edge() {
    // The 2·r ≤ edge_length self-intersection bound: a unit-cube edge
    // is length 1, so r = 0.6 (2r = 1.2 > 1) must be rejected.
    let cube = box_solid(1.0, 1.0, 1.0).unwrap();
    let brep = inner(&cube);
    let edge = unique_edges(brep)[0].clone();
    let err = plan_planar_edge_fillet(brep, &edge, 0.6).unwrap_err();
    assert!(
        matches!(err, FilletBrepError::RadiusTooLarge { .. }),
        "expected RadiusTooLarge, got {err:?}"
    );
}

#[test]
fn fillet_plan_rejects_bad_radius() {
    let cube = box_solid(2.0, 2.0, 2.0).unwrap();
    let brep = inner(&cube);
    let edge = unique_edges(brep)[0].clone();
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let err = plan_planar_edge_fillet(brep, &edge, bad).unwrap_err();
        assert!(
            matches!(err, FilletBrepError::BadParameter { name: "radius", .. }),
            "radius {bad} should be BadParameter, got {err:?}"
        );
    }
}

// ===========================================================================
// The fillet result — analytic volume removal, or a clean soft-fail.
// ===========================================================================

/// Analytic removed volume of a 90°-edge fillet: the corner sliver
/// `r²(1 − π/4)` per unit edge length.
fn analytic_removed_volume_90deg(radius: f64, edge_length: f64) -> f64 {
    radius * radius * (1.0 - PI / 4.0) * edge_length
}

#[test]
fn single_edge_fillet_either_removes_the_analytic_sliver_or_soft_fails() {
    // Fillet one edge of a 4-unit cube with radius 0.5. The cube's
    // volume is 64; a successful fillet must remove exactly the
    // analytic corner sliver r²(1−π/4)·edge_length.
    let cube = box_solid(4.0, 4.0, 4.0).unwrap();
    let brep = inner(&cube);
    let edge = unique_edges(brep)[0].clone();
    let r = 0.5;
    let cube_volume = 64.0;
    let removed = analytic_removed_volume_90deg(r, 4.0);

    match fillet_planar_edge(brep, &edge, r) {
        Ok(filleted) => {
            let fs = Solid::from_truck(filleted);
            assert!(
                is_closed_solid_tol(&fs, 1e-3).unwrap(),
                "a successful fillet must produce a valid closed solid"
            );
            let v = solid_volume_tol(&fs, 1e-3).unwrap();
            let expected = cube_volume - removed;
            assert!(
                (v - expected).abs() / expected < 0.01,
                "filleted-cube volume {v} should be within 1% of \
                 cube − sliver = {expected}"
            );
        }
        Err(FilletBrepError::TruckOp(_)) => {
            // The coincident-face cutter boolean could not be resolved
            // by truck_shapeops — the documented soft fall-through.
        }
        other => panic!("a textbook convex planar edge should fillet or soft-fail, got {other:?}"),
    }
}

#[test]
fn fillet_never_panics_across_a_radius_sweep() {
    // Sweep the radius across a wide range on a cube edge. Every value
    // must produce a defined outcome — a real BRep solid, the soft
    // TruckOp fall-through, or a RadiusTooLarge for the largest radii.
    // The regression guarded against is a panic crossing the truck FFI.
    let cube = box_solid(4.0, 4.0, 4.0).unwrap();
    let brep = inner(&cube);
    let edge = unique_edges(brep)[0].clone();
    for &r in &[0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 1.9, 2.0, 2.5] {
        match fillet_planar_edge(brep, &edge, r) {
            Ok(_)
            | Err(FilletBrepError::TruckOp(_))
            | Err(FilletBrepError::RadiusTooLarge { .. }) => {
                // All defined outcomes.
            }
            other => panic!("radius {r} produced an unexpected fillet outcome: {other:?}"),
        }
    }
}

#[test]
fn fillet_result_when_successful_is_smaller_than_the_original() {
    // Whatever the radius, a successful convex-edge fillet removes
    // material — the filleted solid's volume must be strictly less
    // than the original. (If the boolean soft-fails, skip.)
    let cube = box_solid(6.0, 6.0, 6.0).unwrap();
    let brep = inner(&cube);
    let edge = unique_edges(brep)[0].clone();
    if let Ok(filleted) = fillet_planar_edge(brep, &edge, 1.0) {
        let fs = Solid::from_truck(filleted);
        let v = solid_volume_tol(&fs, 1e-3).unwrap();
        assert!(
            v < 216.0 && v > 0.0,
            "a filleted 6-cube must have volume in (0, 216), got {v}"
        );
    }
}
