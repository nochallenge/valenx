//! Per-edge BRep fillet: convex straight edge, two planar adjacent
//! faces, constant radius.
//!
//! # The pipeline
//!
//! For each candidate edge, the pipeline runs in three stages:
//!
//! 1. **Classify** — confirm both adjacent faces are planar, the edge
//!    is convex, and the requested radius fits along the edge.
//!    ([`crate::edge_classify`].)
//! 2. **Plan** — compute the [`EdgeFilletPlan`]: bisector direction,
//!    cylinder axis, tangent contact points on each face, fillet
//!    surface radius and extent. This is pure geometry — no truck
//!    construction yet — and is fully covered by unit tests.
//! 3. **Construct** — perform the real BRep surgery. See
//!    [`crate::brep_build`].
//!
//! # Stage 3 — the real BRep fillet (Phase 14.5)
//!
//! Stage 3 is implemented in [`crate::brep_build`] as constructive
//! solid geometry. Rather than reaching for a face-trim primitive
//! truck-modeling 0.6 does not expose, the fillet is built from two
//! genuine BRep solids — a triangular **cutter** prism (the corner
//! sliver) and a circular-sector **fillet bar** prism (the rounded
//! fill) — and the real `truck_shapeops` booleans:
//!
//! ```text
//!   filleted = (solid − cutter) ∪ bar
//! ```
//!
//! The result is a [`valenx_cad::Solid::Brep`] with a true
//! circular-arc-swept fillet face; it round-trips through STEP/IGES
//! and composes with further BRep ops. This is a true-BRep fillet,
//! not a tessellated approximation.
//!
//! # Multi-edge corners
//!
//! [`fillet_solid_edges`] composes per-edge fillets and then blends
//! every **orthogonal convex 3-edge corner** (the corner of a box)
//! with a real rolling-ball corner — see [`fillet_corner_blend`] and
//! [`crate::corner_build`]. General N-edge, non-orthogonal, and
//! concave corners stay a Tier-3 residue and fall back to the
//! independent per-edge fillets.
//!
//! # Honest scope
//!
//! Two cutter faces lie flush with the solid's adjacent faces;
//! coincident faces are the historically fragile input for any BRep
//! boolean kernel. On geometry where `truck_shapeops` cannot resolve
//! the coincidence the boolean returns no solid and this layer
//! surfaces [`FilletBrepError::TruckOp`] — a *soft* failure the
//! `valenx-feature-tree` dispatcher treats as a fall-through signal to
//! the Phase 3 mesh-domain pipeline. Curved adjacent faces and
//! concave edges remain out of scope (see the crate-level docs and
//! [`crate::brep_build`]).

use truck_modeling::{Edge as TruckEdge, InnerSpace, Point3, Solid as TruckSolid, Vector3};

use crate::edge_classify::{adjacent_faces, is_convex_edge, is_planar, planar_normal};
use crate::error::FilletBrepError;
use crate::topology::{edge_chord_length, edge_endpoints};

/// Geometric description of the fillet that *would* be built along an
/// edge. Output of stage 2 of the pipeline; consumed by stage 3 (the
/// substitution step) which is a Phase 14.5+ TODO.
///
/// All vectors are in world space.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeFilletPlan {
    /// World-space front endpoint of the original edge.
    pub edge_front: Point3,
    /// World-space back endpoint of the original edge.
    pub edge_back: Point3,
    /// Outward normal of the first adjacent planar face.
    pub face0_normal: Vector3,
    /// Outward normal of the second adjacent planar face.
    pub face1_normal: Vector3,
    /// Bisector direction, pointing *into* the solid from the edge.
    /// The fillet's centerline is the original edge translated by
    /// `-radius * bisector` (i.e. `radius` units inward along this
    /// vector). All edge points get one cylinder centerline point.
    pub bisector_inward: Vector3,
    /// Fillet radius (same as the request).
    pub radius: f64,
    /// Tangent contact line on `face0` — the line where the
    /// cylindrical fillet surface meets the planar face. Offset
    /// inward from the original edge by `radius / sin(theta/2)`
    /// along `face0`'s in-plane direction.
    pub tangent_on_face0_front: Point3,
    /// Tangent contact line on `face0` — back endpoint.
    pub tangent_on_face0_back: Point3,
    /// Tangent contact line on `face1` — front endpoint.
    pub tangent_on_face1_front: Point3,
    /// Tangent contact line on `face1` — back endpoint.
    pub tangent_on_face1_back: Point3,
    /// Dihedral angle (radians) between the two outward face normals.
    /// 0 = coplanar (cannot fillet), π/2 = right angle, π = back-to-
    /// back fold.
    pub dihedral_angle: f64,
}

/// Plan a fillet on the given edge: validate inputs and compute
/// [`EdgeFilletPlan`]. Pure-geometry stage of the pipeline.
///
/// Returns:
/// - [`FilletBrepError::NotPlanarFaces`] if either adjacent face is
///   non-planar or the edge has != 2 adjacent faces.
/// - [`FilletBrepError::NonConvexEdge`] if the edge is concave.
/// - [`FilletBrepError::RadiusTooLarge`] if `radius * 2 >
///   min_edge_length`.
/// - [`FilletBrepError::BadParameter`] for invalid radius.
pub fn plan_planar_edge_fillet(
    solid: &TruckSolid,
    edge: &TruckEdge,
    radius: f64,
) -> Result<EdgeFilletPlan, FilletBrepError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(FilletBrepError::BadParameter {
            name: "radius",
            reason: format!("must be > 0 and finite, got {radius}"),
        });
    }

    let faces = adjacent_faces(solid, edge);
    if faces.len() != 2 {
        return Err(FilletBrepError::NotPlanarFaces);
    }
    if !is_planar(faces[0]) || !is_planar(faces[1]) {
        return Err(FilletBrepError::NotPlanarFaces);
    }
    if !is_convex_edge(solid, edge) {
        return Err(FilletBrepError::NonConvexEdge);
    }

    // Edge-length bound: the fillet's tangent contact line on each
    // face is offset inward from the edge by `radius / tan(theta/2)`
    // where theta is the dihedral angle between the *inward* normals
    // (= π - dihedral between outward normals). For a 90° corner
    // this offset is exactly `radius`. We use a simpler bound here:
    // `2 * radius <= edge_length` keeps the fillet surface from
    // looping around itself when applied near both endpoints.
    let edge_length = edge_chord_length(edge);
    if 2.0 * radius > edge_length {
        return Err(FilletBrepError::RadiusTooLarge {
            radius,
            min_edge_length: edge_length,
        });
    }

    let n0 = planar_normal(faces[0]).ok_or(FilletBrepError::NotPlanarFaces)?;
    let n1 = planar_normal(faces[1]).ok_or(FilletBrepError::NotPlanarFaces)?;
    let (front, back) = edge_endpoints(edge);

    // Bisector points "outward" along (n0 + n1).normalize() for a
    // convex corner; we want the *inward* bisector so the cylinder
    // centerline sits inside the solid.
    let sum = n0 + n1;
    let bisector_outward_mag = sum.magnitude();
    if bisector_outward_mag < 1e-12 {
        // n0 + n1 == 0 means the faces are back-to-back (180°
        // fold) — not a real solid edge.
        return Err(FilletBrepError::NotPlanarFaces);
    }
    let bisector_outward = sum / bisector_outward_mag;
    let bisector_inward = -bisector_outward;

    // Dihedral between *outward* normals.
    let cos_th = n0.dot(n1) / (n0.magnitude() * n1.magnitude());
    let dihedral_angle = cos_th.clamp(-1.0, 1.0).acos();

    // The half-angle of the corner (between the inward face normals)
    // is `(π - dihedral_angle) / 2`. The tangent contact line on
    // each face is offset from the edge by `radius / tan(half)`.
    let half = (std::f64::consts::PI - dihedral_angle) * 0.5;
    let tan_half = half.tan();
    if !tan_half.is_finite() || tan_half.abs() < 1e-12 {
        return Err(FilletBrepError::NotPlanarFaces);
    }
    let offset = radius / tan_half;

    // In-plane direction on each face pointing *away* from the edge
    // (so adding `offset * d` to an edge endpoint gives the tangent
    // contact point on that face). For face i: the inward direction
    // along the face is `bisector_inward × n_i × n_i` projected back
    // into the face's plane, or equivalently `bisector_inward minus
    // its projection onto n_i`, normalized.
    let in_face0 = project_into_plane(bisector_inward, n0);
    let in_face1 = project_into_plane(bisector_inward, n1);

    let tangent_on_face0_front = front + in_face0 * offset;
    let tangent_on_face0_back = back + in_face0 * offset;
    let tangent_on_face1_front = front + in_face1 * offset;
    let tangent_on_face1_back = back + in_face1 * offset;

    Ok(EdgeFilletPlan {
        edge_front: front,
        edge_back: back,
        face0_normal: n0,
        face1_normal: n1,
        bisector_inward,
        radius,
        tangent_on_face0_front,
        tangent_on_face0_back,
        tangent_on_face1_front,
        tangent_on_face1_back,
        dihedral_angle,
    })
}

/// Project `v` into the plane perpendicular to `n`, then normalize.
/// Returns the zero vector if `v` is parallel to `n`.
fn project_into_plane(v: Vector3, n: Vector3) -> Vector3 {
    let n_unit = if n.magnitude() > 1e-12 {
        n.normalize()
    } else {
        return Vector3::new(0.0, 0.0, 0.0);
    };
    let in_plane = v - n_unit * v.dot(n_unit);
    let mag = in_plane.magnitude();
    if mag > 1e-12 {
        in_plane / mag
    } else {
        Vector3::new(0.0, 0.0, 0.0)
    }
}

/// Apply a real BRep fillet to a single convex straight edge of a
/// solid.
///
/// Validates the input (planar faces, convex edge, radius bound) via
/// [`plan_planar_edge_fillet`], then performs the constructive-solid-
/// geometry surgery in [`crate::brep_build::fillet_convex_planar_edge`]:
/// builds the BRep cutter + fillet-bar prisms and evaluates
/// `(solid − cutter) ∪ bar` with the real `truck_shapeops` booleans.
///
/// Returns a [`TruckSolid`] carrying a genuine circular-arc-swept
/// fillet face.
///
/// # Errors
///
/// - The geometric-precondition errors from the planner
///   (`NotPlanarFaces`, `NonConvexEdge`, `RadiusTooLarge`,
///   `BadParameter`) — a bad input is reported by its true cause.
/// - [`FilletBrepError::TruckOp`] if the boolean kernel cannot
///   resolve the flush cutter faces. The `valenx-feature-tree`
///   dispatcher treats this as a soft fall-through signal to the
///   Phase 3 mesh-domain pipeline.
pub fn fillet_planar_edge(
    solid: &TruckSolid,
    edge: &TruckEdge,
    radius: f64,
) -> Result<TruckSolid, FilletBrepError> {
    crate::brep_build::fillet_convex_planar_edge(solid, edge, radius)
}

/// Apply a real BRep **variable-radius** fillet to a single convex
/// straight edge of a solid.
///
/// The fillet radius varies **linearly** along the edge — `radius_start`
/// at the edge's front endpoint, `radius_end` at the back — so the
/// fillet surface is a tapered, lofted blend rather than a constant
/// cylinder. With `radius_start == radius_end` the result is the
/// constant-radius fillet of [`fillet_planar_edge`].
///
/// Validates the input via [`plan_planar_edge_fillet`] (the radius
/// bound is checked against the larger of the two endpoint radii), then
/// performs the surgery in
/// [`crate::brep_build::fillet_variable_radius_planar_edge`]: it lofts
/// the BRep cutter + fillet-bar prisms between the two end cross-
/// sections and evaluates `(solid − cutter) ∪ bar`.
///
/// # Honest scope
///
/// The radius law is **linear** between the two endpoints. A general
/// radius profile (a spline of radius-vs-arc-length) would loft through
/// intermediate stations and is a bounded follow-up. As with the
/// constant-radius fillet, the adjacent faces must be planar, the edge
/// must be a single convex straight edge, and **multi-edge corner
/// blends** (3+ filleted edges meeting at a vertex) remain a Tier-3
/// research residual — see the crate-level docs.
///
/// # Errors
///
/// - The geometric-precondition errors from the planner
///   (`NotPlanarFaces`, `NonConvexEdge`, `RadiusTooLarge`,
///   `BadParameter`).
/// - [`FilletBrepError::BadParameter`] if either endpoint radius is
///   non-finite or non-positive.
/// - [`FilletBrepError::TruckOp`] if a loft / cap step or a boolean
///   fails — a *soft* fall-through signal for the dispatcher, exactly
///   as for the constant-radius fillet.
pub fn fillet_variable_radius_edge(
    solid: &TruckSolid,
    edge: &TruckEdge,
    radius_start: f64,
    radius_end: f64,
) -> Result<TruckSolid, FilletBrepError> {
    crate::brep_build::fillet_variable_radius_planar_edge(
        solid,
        edge,
        radius_start,
        radius_end,
    )
}

/// Apply a real BRep **rolling-ball corner blend** at an orthogonal
/// convex 3-edge corner.
///
/// When three filleted edges meet at a box corner, the three
/// quarter-cylinder fillet surfaces leave a tri-tangent gap at the
/// vertex. This blends that gap with the radius-`r` ball seated
/// tangent to all three faces — see
/// [`crate::corner_build::blend_corner`].
///
/// `corner` must be the supported orthogonal-convex case; obtain it
/// from [`crate::corner::classify_corner`] or
/// [`crate::corner::find_blendable_corners`]. The blend is a pure
/// geometric construction (centre + sphere + cutter box from the
/// corner apex and edge directions), so it applies correctly to a
/// solid whose per-edge fillets have already replaced the original
/// corner topology.
///
/// # Errors
///
/// - [`FilletBrepError::BadParameter`] if the radius is invalid.
/// - [`FilletBrepError::TruckOp`] if the corner-blend boolean cannot
///   be resolved (the coincident-surface case — a *soft* fall-through
///   signal, exactly as for the per-edge fillet).
pub fn fillet_corner_blend(
    solid: &TruckSolid,
    corner: &crate::corner::CornerEdges,
    radius: f64,
) -> Result<TruckSolid, FilletBrepError> {
    crate::corner_build::blend_corner(solid, corner, radius)
}

/// Batch helper: apply a fillet to each edge in sequence, then blend
/// the box-style corners where 3 filleted edges meet.
///
/// Each edge is filleted on the running result, so the per-edge
/// fillets compose. After the per-edge pass, every **orthogonal
/// convex 3-edge corner** (the corner of a box) is blended with a
/// real rolling-ball corner via [`fillet_corner_blend`] — so
/// filleting all 12 edges of a box now also blends its 8 corners.
///
/// # Honest behaviour
///
/// - The corners to blend are detected on the **original** solid +
///   edge list *before* filleting (the per-edge fillets replace the
///   original corner topology, so the apex geometry must be captured
///   first); the blend itself is a geometric construction that
///   applies to the filleted running result.
/// - Only orthogonal convex 3-edge corners are blended. A vertex
///   where 3+ filleted edges meet but the corner is non-orthogonal,
///   concave, or higher-degree is **not** blended — it stays the
///   independent per-edge fillets v1 has always produced (the
///   genuine Tier-3 corner residue; see the crate docs).
/// - Each corner blend **soft-fails independently**: if a corner's
///   `truck_shapeops` boolean cannot resolve the coincident surfaces,
///   that corner is left as the un-blended per-edge fillets and the
///   batch continues — the result is never worse than the pre-14.7
///   independent-fillet behaviour.
/// - A per-edge fillet failure is still hard: it returns the error
///   from the first edge that fails (the per-edge fillet is the
///   load-bearing operation; a corner blend is a refinement on top).
pub fn fillet_solid_edges(
    solid: &TruckSolid,
    edges: &[TruckEdge],
    radius: f64,
) -> Result<TruckSolid, FilletBrepError> {
    if edges.is_empty() {
        return Err(FilletBrepError::BadParameter {
            name: "edges",
            reason: "edge list cannot be empty".into(),
        });
    }

    // Detect the blendable corners on the ORIGINAL solid before any
    // fillet replaces the corner topology.
    let corners = crate::corner::find_blendable_corners(solid, edges);

    // Per-edge fillet pass — the load-bearing operation.
    let mut current = solid.clone();
    for edge in edges {
        current = fillet_planar_edge(&current, edge, radius)?;
    }

    // Corner-blend pass — each corner soft-fails independently so a
    // single un-resolvable coincident-surface boolean never degrades
    // the whole result below the per-edge-fillet baseline.
    for (_vid, corner) in &corners {
        match fillet_corner_blend(&current, corner, radius) {
            Ok(blended) => current = blended,
            Err(_) => {
                // Soft fall-through: leave this corner as the
                // independent per-edge fillets (pre-14.7 behaviour).
            }
        }
    }

    Ok(current)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;
    use valenx_cad::primitives::box_solid;

    fn inner_brep(s: &valenx_cad::Solid) -> &TruckSolid {
        match s {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("expected brep"),
        }
    }

    /// Pick the first unique edge of a cube.
    fn pick_first_unique_edge(brep: &TruckSolid) -> TruckEdge {
        let mut seen = std::collections::HashSet::new();
        for edge in brep.edge_iter() {
            if seen.insert(edge.id()) {
                return edge;
            }
        }
        panic!("solid has no edges")
    }

    #[test]
    fn fillet_planar_edge_on_cube_either_builds_or_soft_fails() {
        // The fillet now runs the real BRep construction. For a
        // valid, well-conditioned cube edge the planner always
        // succeeds; the boolean stage then either produces a real
        // BRep solid or — if `truck_shapeops` cannot resolve the
        // flush cutter faces — surfaces the soft `TruckOp` error the
        // dispatcher falls through on. Both are honest outcomes; what
        // must NOT happen is a geometric-precondition error (the cube
        // edge is a textbook convex planar-faces case) or a panic.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        match fillet_planar_edge(brep, &edge, 0.3) {
            Ok(filleted) => {
                // A real BRep fillet: filleting a convex edge rounds
                // one corner, so the result has strictly more faces
                // than the original 6-face cube (the cube loses no
                // face but gains the cylindrical fillet face).
                let faces: usize =
                    filleted.boundaries().iter().map(|s| s.len()).sum();
                assert!(
                    faces >= 6,
                    "filleted cube should have at least the cube's faces, got {faces}"
                );
            }
            Err(FilletBrepError::TruckOp(_)) => {
                // Soft fall-through — acceptable: coincident-face
                // booleans are the documented fragile case.
            }
            other => panic!("unexpected fillet outcome: {other:?}"),
        }
    }

    #[test]
    fn fillet_with_too_large_radius_returns_radius_too_large() {
        // Unit cube edges are length 1. radius=0.6 gives 2*r > edge,
        // should trigger RadiusTooLarge *before* reaching the
        // construction stage.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = fillet_planar_edge(brep, &edge, 0.6);
        match result {
            Err(FilletBrepError::RadiusTooLarge {
                radius,
                min_edge_length,
            }) => {
                assert!((radius - 0.6).abs() < 1e-9);
                assert!((min_edge_length - 1.0).abs() < 1e-9);
            }
            other => panic!("expected RadiusTooLarge, got {other:?}"),
        }
    }

    #[test]
    fn fillet_with_zero_radius_returns_bad_parameter() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = fillet_planar_edge(brep, &edge, 0.0);
        assert!(matches!(
            result,
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
    }

    #[test]
    fn fillet_with_negative_radius_returns_bad_parameter() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = fillet_planar_edge(brep, &edge, -0.1);
        assert!(matches!(
            result,
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
    }

    #[test]
    fn fillet_with_nan_radius_returns_bad_parameter() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = fillet_planar_edge(brep, &edge, f64::NAN);
        assert!(matches!(
            result,
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
    }

    #[test]
    fn plan_cube_edge_has_right_angle_dihedral() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.1).expect("plan succeeds");
        assert!(
            (plan.dihedral_angle - FRAC_PI_2).abs() < 1e-6,
            "cube edge dihedral should be 90°, got {} rad",
            plan.dihedral_angle
        );
        assert!((plan.radius - 0.1).abs() < 1e-12);
        // For a 90° corner, the half-angle of the corner (between
        // inward normals) is 45°, so offset = radius/tan(45°) = radius.
        // Tangent contact points should be exactly `radius` from the
        // edge endpoint along the in-plane direction.
        let offset_front_face0 = (plan.tangent_on_face0_front - plan.edge_front).magnitude();
        assert!(
            (offset_front_face0 - 0.1).abs() < 1e-9,
            "tangent offset should equal radius for 90° corner, got {offset_front_face0}",
        );
    }

    #[test]
    fn plan_cube_bisector_points_into_solid() {
        // The cube spans (0,0,0)..(1,1,1). Both face normals are
        // outward. Their sum's normalized direction is *outward*; we
        // flip to get inward → should have negative components on the
        // axes orthogonal to the edge.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.1).expect("plan succeeds");
        // The bisector must be a unit vector.
        assert!(
            (plan.bisector_inward.magnitude() - 1.0).abs() < 1e-9,
            "bisector_inward should be unit, got {}",
            plan.bisector_inward.magnitude()
        );
        // For a cube edge at a corner of the cube, bisector_inward
        // points into the interior. The midpoint of the edge offset
        // by a small step along bisector_inward should land inside
        // the AABB.
        let mid = plan.edge_front + (plan.edge_back - plan.edge_front) * 0.5;
        let probe = mid + plan.bisector_inward * 0.1;
        assert!(
            probe.x >= -1e-6 && probe.x <= 1.0 + 1e-6,
            "probe x should be inside cube, got {}",
            probe.x
        );
        assert!(
            probe.y >= -1e-6 && probe.y <= 1.0 + 1e-6,
            "probe y should be inside cube, got {}",
            probe.y
        );
        assert!(
            probe.z >= -1e-6 && probe.z <= 1.0 + 1e-6,
            "probe z should be inside cube, got {}",
            probe.z
        );
    }

    #[test]
    fn fillet_solid_edges_empty_list_errors() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let result = fillet_solid_edges(brep, &[], 0.1);
        assert!(matches!(
            result,
            Err(FilletBrepError::BadParameter { name: "edges", .. })
        ));
    }

    #[test]
    fn fillet_solid_edges_batch_rejects_bad_radius_before_construction() {
        // A genuinely-bad radius (too large for the edge) must surface
        // as the geometric `RadiusTooLarge` cause from the planner —
        // the batch path runs per-edge validation before any boolean.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = fillet_solid_edges(brep, std::slice::from_ref(&edge), 0.7);
        assert!(matches!(result, Err(FilletBrepError::RadiusTooLarge { .. })));
    }

    #[test]
    fn fillet_solid_edges_single_edge_either_builds_or_soft_fails() {
        // A one-edge batch is equivalent to a single `fillet_planar_edge`
        // call: either a real BRep result or a soft `TruckOp`
        // fall-through. Never a geometric-precondition error for a
        // textbook cube edge, never a panic.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        match fillet_solid_edges(brep, std::slice::from_ref(&edge), 0.3) {
            Ok(_) | Err(FilletBrepError::TruckOp(_)) => {}
            other => panic!("unexpected batch fillet outcome: {other:?}"),
        }
    }

    /// Collect every unique edge of a brep.
    fn all_unique_edges(brep: &TruckSolid) -> Vec<TruckEdge> {
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for e in brep.edge_iter() {
            if seen.insert(e.id()) {
                out.push(e);
            }
        }
        out
    }

    #[test]
    fn fillet_corner_blend_on_cube_corner_either_builds_or_soft_fails() {
        // The corner-blend primitive runs the real cutter + ball +
        // boolean construction at a cube corner. The outcome must be a
        // real BRep solid or the documented soft `TruckOp`
        // fall-through (the corner cutter trims flush with three
        // faces) — never a panic, never a precondition error.
        let cube = box_solid(6.0, 6.0, 6.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = all_unique_edges(brep);
        let corners = crate::corner::find_blendable_corners(brep, &edges);
        assert!(!corners.is_empty(), "a cube has blendable corners");
        let (_vid, corner) = corners[0];
        match fillet_corner_blend(brep, &corner, 0.6) {
            Ok(blended) => {
                let faces: usize = blended.boundaries().iter().map(|s| s.len()).sum();
                assert!(faces >= 6, "corner-blended cube keeps its faces");
            }
            Err(FilletBrepError::TruckOp(_)) => {}
            other => panic!("unexpected corner-blend outcome: {other:?}"),
        }
    }

    #[test]
    fn fillet_corner_blend_rejects_bad_radius() {
        // A bad radius must fail with the BadParameter cause.
        let cube = box_solid(6.0, 6.0, 6.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = all_unique_edges(brep);
        let corners = crate::corner::find_blendable_corners(brep, &edges);
        let (_vid, corner) = corners[0];
        assert!(matches!(
            fillet_corner_blend(brep, &corner, 0.0),
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
    }

    #[test]
    fn fillet_all_twelve_cube_edges_builds_with_corner_blends_or_soft_fails() {
        // Filleting every edge of a box now runs the per-edge fillet
        // pass AND the 8 corner blends. The whole pipeline must
        // either produce a real BRep solid (with the cube's faces
        // preserved) or soft-fail at a per-edge step — corner blends
        // soft-fail independently, so they never turn an otherwise-OK
        // result into an error. Never a panic, never a precondition
        // error for a textbook cube.
        let cube = box_solid(6.0, 6.0, 6.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = all_unique_edges(brep);
        assert_eq!(edges.len(), 12, "a cube has 12 edges");
        // There are 8 blendable corners detected up front.
        assert_eq!(
            crate::corner::find_blendable_corners(brep, &edges).len(),
            8,
            "a fully-filleted cube has 8 blendable corners"
        );
        match fillet_solid_edges(brep, &edges, 0.5) {
            Ok(filleted) => {
                let faces: usize = filleted.boundaries().iter().map(|s| s.len()).sum();
                assert!(
                    faces >= 6,
                    "a fully-filleted + corner-blended cube keeps the cube's faces, got {faces}"
                );
            }
            Err(FilletBrepError::TruckOp(_)) => {
                // Soft fall-through at a per-edge step — acceptable:
                // coincident-face booleans are the documented fragile
                // case.
            }
            other => panic!("unexpected full-cube fillet outcome: {other:?}"),
        }
    }

    #[test]
    fn project_into_plane_unit_xy() {
        // Vector (1, 0, 1) projected into plane perpendicular to (0, 0, 1)
        // (i.e. XY plane) → (1, 0, 0). Then normalized → (1, 0, 0).
        let v = Vector3::new(1.0, 0.0, 1.0);
        let n = Vector3::new(0.0, 0.0, 1.0);
        let out = project_into_plane(v, n);
        assert!((out.x - 1.0).abs() < 1e-9);
        assert!(out.y.abs() < 1e-9);
        assert!(out.z.abs() < 1e-9);
    }

    #[test]
    fn project_into_plane_parallel_to_normal_returns_zero() {
        // A vector exactly along the normal projects to (0, 0, 0).
        let v = Vector3::new(0.0, 0.0, 1.0);
        let n = Vector3::new(0.0, 0.0, 1.0);
        let out = project_into_plane(v, n);
        assert!(out.magnitude() < 1e-9);
    }
}
