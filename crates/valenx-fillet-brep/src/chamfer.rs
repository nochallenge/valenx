//! Per-edge BRep chamfer for v1: convex straight edge, two planar
//! adjacent faces, flat bevel of constant width.
//!
//! Same pipeline shape as [`crate::fillet`] but the inserted surface
//! is a *planar quad* (the bevel) rather than a cylindrical patch.
//! Same v1 limitation: stage 3 substitution returns
//! [`FilletBrepError::TruckSubstitutionUnavailable`] until upstream
//! truck-modeling exposes face-trim-and-substitute.

use truck_modeling::{Edge as TruckEdge, InnerSpace, Point3, Solid as TruckSolid, Vector3};

use crate::edge_classify::{adjacent_faces, is_convex_edge, is_planar, planar_normal};
use crate::error::FilletBrepError;
use crate::topology::{edge_chord_length, edge_endpoints};

/// Geometric description of the chamfer that *would* be built along
/// an edge. Output of the planning stage.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeChamferPlan {
    /// World-space front endpoint of the original edge.
    pub edge_front: Point3,
    /// World-space back endpoint of the original edge.
    pub edge_back: Point3,
    /// Outward normal of the first adjacent planar face.
    pub face0_normal: Vector3,
    /// Outward normal of the second adjacent planar face.
    pub face1_normal: Vector3,
    /// Chamfer offset distance (same as the request).
    pub distance: f64,
    /// Bevel quad corner: face0 side, edge-front end.
    pub bevel_corner_face0_front: Point3,
    /// Bevel quad corner: face0 side, edge-back end.
    pub bevel_corner_face0_back: Point3,
    /// Bevel quad corner: face1 side, edge-front end.
    pub bevel_corner_face1_front: Point3,
    /// Bevel quad corner: face1 side, edge-back end.
    pub bevel_corner_face1_back: Point3,
}

/// Plan a chamfer on the given edge: validate inputs and compute
/// [`EdgeChamferPlan`]. Pure-geometry stage.
///
/// Errors mirror [`crate::fillet::plan_planar_edge_fillet`].
pub fn plan_planar_edge_chamfer(
    solid: &TruckSolid,
    edge: &TruckEdge,
    distance: f64,
) -> Result<EdgeChamferPlan, FilletBrepError> {
    if !distance.is_finite() || distance <= 0.0 {
        return Err(FilletBrepError::BadParameter {
            name: "distance",
            reason: format!("must be > 0 and finite, got {distance}"),
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

    let edge_length = edge_chord_length(edge);
    if 2.0 * distance > edge_length {
        return Err(FilletBrepError::RadiusTooLarge {
            radius: distance,
            min_edge_length: edge_length,
        });
    }

    let n0 = planar_normal(faces[0]).ok_or(FilletBrepError::NotPlanarFaces)?;
    let n1 = planar_normal(faces[1]).ok_or(FilletBrepError::NotPlanarFaces)?;
    let (front, back) = edge_endpoints(edge);

    // In-plane direction pointing inward along face0 (away from the
    // edge into the face), reached by removing the n1-component from
    // (n0+n1)/2 and flipping. We compute it as `-bisector_outward`
    // projected into each face's plane, normalized.
    let sum = n0 + n1;
    let bis_mag = sum.magnitude();
    if bis_mag < 1e-12 {
        return Err(FilletBrepError::NotPlanarFaces);
    }
    let bisector_inward = -sum / bis_mag;

    let in_face0 = project_into_plane(bisector_inward, n0);
    let in_face1 = project_into_plane(bisector_inward, n1);

    let bevel_corner_face0_front = front + in_face0 * distance;
    let bevel_corner_face0_back = back + in_face0 * distance;
    let bevel_corner_face1_front = front + in_face1 * distance;
    let bevel_corner_face1_back = back + in_face1 * distance;

    Ok(EdgeChamferPlan {
        edge_front: front,
        edge_back: back,
        face0_normal: n0,
        face1_normal: n1,
        distance,
        bevel_corner_face0_front,
        bevel_corner_face0_back,
        bevel_corner_face1_front,
        bevel_corner_face1_back,
    })
}

/// Project `v` into the plane perpendicular to `n`, then normalize.
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

/// Apply a chamfer to a single edge.
///
/// **v1 status:** identical to [`crate::fillet::fillet_planar_edge`]
/// — validates the input then returns
/// [`FilletBrepError::TruckSubstitutionUnavailable`]. The dispatcher
/// in `valenx-feature-tree::ops::chamfer` treats this as a soft
/// error and falls through to the Phase 3 mesh-domain chamfer.
pub fn chamfer_planar_edge(
    solid: &TruckSolid,
    edge: &TruckEdge,
    distance: f64,
) -> Result<TruckSolid, FilletBrepError> {
    let _plan = plan_planar_edge_chamfer(solid, edge, distance)?;
    Err(FilletBrepError::TruckSubstitutionUnavailable)
}

/// Batch helper: apply a chamfer to each edge in sequence.
///
/// Returns the soft error from the first edge that fails. v1 does
/// not handle corner blending when 3+ chamfered edges meet at a
/// vertex (the corner is left as a tri-tangent point).
pub fn chamfer_solid_edges(
    solid: &TruckSolid,
    edges: &[TruckEdge],
    distance: f64,
) -> Result<TruckSolid, FilletBrepError> {
    if edges.is_empty() {
        return Err(FilletBrepError::BadParameter {
            name: "edges",
            reason: "edge list cannot be empty".into(),
        });
    }
    let mut current = solid.clone();
    for edge in edges {
        current = chamfer_planar_edge(&current, edge, distance)?;
    }
    Ok(current)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::primitives::box_solid;

    fn inner_brep(s: &valenx_cad::Solid) -> &TruckSolid {
        match s {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("expected brep"),
        }
    }

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
    fn chamfer_planar_edge_returns_substitution_unavailable() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = chamfer_planar_edge(brep, &edge, 0.1);
        assert!(matches!(
            result,
            Err(FilletBrepError::TruckSubstitutionUnavailable)
        ));
    }

    #[test]
    fn chamfer_too_large_returns_radius_too_large() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = chamfer_planar_edge(brep, &edge, 0.6);
        assert!(matches!(
            result,
            Err(FilletBrepError::RadiusTooLarge { .. })
        ));
    }

    #[test]
    fn chamfer_zero_distance_returns_bad_parameter() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let result = chamfer_planar_edge(brep, &edge, 0.0);
        assert!(matches!(
            result,
            Err(FilletBrepError::BadParameter {
                name: "distance",
                ..
            })
        ));
    }

    #[test]
    fn plan_cube_chamfer_offsets_by_distance() {
        // For a 90° corner, the bevel-corner points should be exactly
        // `distance` from the edge endpoint along the in-plane
        // direction.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_unique_edge(brep);
        let plan = plan_planar_edge_chamfer(brep, &edge, 0.15).expect("plan succeeds");
        let off = (plan.bevel_corner_face0_front - plan.edge_front).magnitude();
        assert!(
            (off - 0.15).abs() < 1e-9,
            "bevel offset should equal distance, got {off}",
        );
        assert!((plan.distance - 0.15).abs() < 1e-12);
    }

    #[test]
    fn chamfer_solid_edges_empty_errors() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let result = chamfer_solid_edges(brep, &[], 0.1);
        assert!(matches!(
            result,
            Err(FilletBrepError::BadParameter { name: "edges", .. })
        ));
    }

    #[test]
    fn chamfer_solid_edges_batch_substitution_unavailable() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let mut edges = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for e in brep.edge_iter() {
            if seen.insert(e.id()) {
                edges.push(e);
                if edges.len() == 3 {
                    break;
                }
            }
        }
        let result = chamfer_solid_edges(brep, &edges, 0.1);
        assert!(matches!(
            result,
            Err(FilletBrepError::TruckSubstitutionUnavailable)
        ));
    }
}
