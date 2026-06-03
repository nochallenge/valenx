//! Edge classification on a truck BRep solid.
//!
//! This module asks the questions Phase 14's fillet/chamfer needs
//! before it tries to operate on an edge:
//!
//! - [`adjacent_faces`]: which faces share this edge? (For a closed
//!   manifold solid this should always be exactly two; we return a
//!   `Vec<&Face>` so degenerate cases — open shells, non-manifold
//!   inputs — fail gracefully later instead of panicking here.)
//! - [`is_planar`]: is this face's surface a plane? (v1 fillets only
//!   handle planar adjacency.)
//! - [`planar_normal`]: extract the normal vector from a planar face.
//! - [`is_convex_edge`]: does this edge form a convex corner? (v1
//!   fillets only handle convex edges — concave-edge fillet would
//!   need a sphere blend.)
//!
//! The classification mirrors the mesh-domain edge classifier in
//! `valenx-fillet::edge_graph` but operates at the BRep level — no
//! tessellation, no vertex welding, no triangle pairs. The price is
//! that truck's API is read-only here; we can ask about faces but
//! can't reach inside and rewrite their boundaries (that's
//! [`crate::fillet`]'s job, which builds new topology rather than
//! mutating the input).

use truck_modeling::{
    Edge as TruckEdge, Face as TruckFace, InnerSpace, ParametricSurface3D, Solid as TruckSolid,
    Surface, Vector3,
};

/// Find every face in `solid` whose boundary includes the given edge
/// (by id — orientation doesn't matter; we want the unoriented edge).
///
/// For a well-formed closed manifold solid this returns exactly 2
/// faces. We return a [`Vec`] rather than a fixed array so non-manifold
/// inputs (open shells, T-junctions) don't panic — the caller can
/// detect `len != 2` and error out cleanly.
pub fn adjacent_faces<'a>(solid: &'a TruckSolid, edge: &TruckEdge) -> Vec<&'a TruckFace> {
    let target_id = edge.id();
    let mut out = Vec::with_capacity(2);
    for face in solid.face_iter() {
        // A face's "absolute boundaries" don't depend on the face's
        // orientation flag — we want the raw wire list so we don't
        // double-count inverted faces.
        let mut found = false;
        for wire in face.absolute_boundaries() {
            for edge_in_wire in wire.edge_iter() {
                if edge_in_wire.id() == target_id {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            out.push(face);
        }
    }
    out
}

/// True if the face's underlying surface is a [`truck_modeling::Plane`].
///
/// Other [`Surface`] variants (B-spline, NURBS, revolution) all return
/// false. v1 only supports plane-plane adjacency for fillets.
pub fn is_planar(face: &TruckFace) -> bool {
    matches!(face.surface(), Surface::Plane(_))
}

/// Extract the outward normal of a planar face.
///
/// The convention is "outward from the solid" — truck's planar face
/// stores its normal aligned with the face's outward orientation, and
/// we honour the face's `orientation()` flag: an inverted face flips
/// the normal so the result still points outward.
///
/// Returns `None` if the face is not planar (call [`is_planar`] first).
pub fn planar_normal(face: &TruckFace) -> Option<Vector3> {
    match face.surface() {
        Surface::Plane(_) => {
            // Plane's normal is constant across (u,v); sample at the
            // origin parameters and flip if the face is inverted.
            let n = face.surface().normal(0.0, 0.0);
            // truck's Face stores the orientation flag *outside* the
            // surface enum; an inverted face presents its boundaries
            // reversed but the surface itself is unchanged. To get the
            // outward normal in the face's *presented* orientation we
            // multiply by the orientation sign.
            //
            // (`Face::orientation()` returns `bool`; true = forward.)
            let sign = if face.orientation() { 1.0 } else { -1.0 };
            Some(n * sign)
        }
        _ => None,
    }
}

/// True if the edge forms a convex corner between its two adjacent
/// planar faces.
///
/// The convexity test follows the same logic as Phase 3's mesh-domain
/// classifier:
///
/// 1. Get the two adjacent face normals `n1`, `n2`. Both must be
///    outward-facing.
/// 2. Take the edge tangent `t = back - front`.
/// 3. Compute the cross product `c = n1 × n2`. For a convex edge, `c`
///    points along `t` (positive dot product); for concave, opposite.
///
/// Returns false if:
/// - The edge has != 2 adjacent faces (open shell / non-manifold).
/// - Either adjacent face is non-planar (we can't take a normal).
/// - The two normals are parallel (the faces are coplanar — there's
///   no edge to fillet).
pub fn is_convex_edge(solid: &TruckSolid, edge: &TruckEdge) -> bool {
    let faces = adjacent_faces(solid, edge);
    if faces.len() != 2 {
        return false;
    }
    let n1 = match planar_normal(faces[0]) {
        Some(n) => n,
        None => return false,
    };
    let n2 = match planar_normal(faces[1]) {
        Some(n) => n,
        None => return false,
    };

    let (front, back) = (edge.front().point(), edge.back().point());
    let t: Vector3 = back - front;
    if t.magnitude2() < 1e-24 {
        // Degenerate edge — can't tell convexity.
        return false;
    }

    let c: Vector3 = n1.cross(n2);
    let dot = c.dot(t);
    // Convex iff cross-product is aligned with tangent. Coplanar
    // faces yield ~0 dot; treat as not convex (caller can't fillet
    // along a flat seam anyway).
    let denom = c.magnitude() * t.magnitude();
    if denom < 1e-12 {
        return false;
    }
    let cos_angle = dot / denom;
    // Use a small positive threshold so genuine "near-coplanar"
    // edges don't sneak into the convex bucket.
    cos_angle > 1e-6
}

/// Dihedral angle (in radians) between the two adjacent planar faces
/// of an edge. Used by the auto-selection helper to filter edges by
/// "sharpness". Returns `None` if the edge has != 2 adjacent faces or
/// either face is non-planar.
///
/// Definition: the angle between the *outward* normals of the two
/// faces. A flat seam (faces coplanar) returns 0; a 90° corner
/// returns π/2; a back-to-back fold returns π.
pub fn edge_dihedral_angle(solid: &TruckSolid, edge: &TruckEdge) -> Option<f64> {
    let faces = adjacent_faces(solid, edge);
    if faces.len() != 2 {
        return None;
    }
    let n1 = planar_normal(faces[0])?;
    let n2 = planar_normal(faces[1])?;
    let cos_t = n1.dot(n2) / (n1.magnitude() * n2.magnitude());
    // acos's domain is [-1, 1]; clamp to avoid NaN from float drift.
    Some(cos_t.clamp(-1.0, 1.0).acos())
}

/// Helper used by tests: build a plane through the origin with the
/// given normal direction (kept module-local; callers shouldn't reach
/// for this from outside).
#[doc(hidden)]
#[cfg(test)]
pub(crate) fn dummy_plane(normal: Vector3) -> truck_modeling::Plane {
    use truck_modeling::{EuclideanSpace, Point3};
    // truck::Plane takes (origin, u_axis, v_axis). For a test plane
    // with `normal`, build any two orthogonal directions that span
    // the plane.
    let n = normal.normalize();
    // Pick an arbitrary axis not parallel to `n`.
    let arbitrary = if n.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = n.cross(arbitrary).normalize();
    let v = n.cross(u).normalize();
    truck_modeling::Plane::new(Point3::origin(), Point3::origin() + u, Point3::origin() + v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::primitives::box_solid;

    /// A unit cube should report 12 distinct edges, each with exactly
    /// 2 adjacent faces. This is the canonical sanity check before we
    /// trust the adjacent-face walk on more complex geometry.
    #[test]
    fn unit_cube_every_edge_has_two_adjacent_faces() {
        let cube = box_solid(1.0, 1.0, 1.0).expect("cube builds");
        let inner = match &cube {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("cube should be brep"),
        };
        let mut seen = std::collections::HashSet::new();
        for edge in inner.edge_iter() {
            if !seen.insert(edge.id()) {
                continue;
            }
            let adj = adjacent_faces(inner, &edge);
            assert_eq!(
                adj.len(),
                2,
                "unit-cube edge should have exactly 2 adjacent faces, got {}",
                adj.len()
            );
        }
    }

    /// A cube has 6 planar faces — every face of the box_solid()
    /// primitive should report `is_planar() == true`.
    #[test]
    fn every_cube_face_is_planar() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let inner = match &cube {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("cube should be brep"),
        };
        for face in inner.face_iter() {
            assert!(is_planar(face), "every cube face should be planar");
            assert!(
                planar_normal(face).is_some(),
                "planar face must have a normal"
            );
        }
    }

    /// Every edge of a cube is convex (faces fold outward).
    #[test]
    fn every_cube_edge_is_convex() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let inner = match &cube {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("cube should be brep"),
        };
        let mut seen = std::collections::HashSet::new();
        let mut tested = 0;
        for edge in inner.edge_iter() {
            if !seen.insert(edge.id()) {
                continue;
            }
            assert!(
                is_convex_edge(inner, &edge),
                "every cube edge should be convex"
            );
            tested += 1;
        }
        assert_eq!(tested, 12, "should have classified 12 unique edges");
    }

    /// Every edge of a cube has a 90° (π/2) dihedral angle between
    /// outward face normals.
    #[test]
    fn cube_edges_have_right_angle_dihedral() {
        use std::f64::consts::FRAC_PI_2;
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let inner = match &cube {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("cube should be brep"),
        };
        let mut seen = std::collections::HashSet::new();
        for edge in inner.edge_iter() {
            if !seen.insert(edge.id()) {
                continue;
            }
            let angle = edge_dihedral_angle(inner, &edge).expect("cube edge has dihedral");
            assert!(
                (angle - FRAC_PI_2).abs() < 1e-6,
                "cube edge dihedral should be 90°, got {angle} rad",
            );
        }
    }

    /// A sphere's faces are non-planar; `is_planar` must return false
    /// for them.
    #[test]
    fn sphere_face_is_not_planar() {
        let s = valenx_cad::primitives::sphere(1.0).unwrap();
        let inner = match &s {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("sphere should be brep"),
        };
        let mut had_non_planar = false;
        for face in inner.face_iter() {
            if !is_planar(face) {
                had_non_planar = true;
                assert!(
                    planar_normal(face).is_none(),
                    "non-planar face should have None for planar_normal"
                );
            }
        }
        assert!(
            had_non_planar,
            "sphere should have at least one non-planar face"
        );
    }

    /// Sanity: `dummy_plane` builds a plane whose normal is the
    /// requested direction (mod sign — cgmath::Plane chooses the
    /// "right" handedness from u × v).
    #[test]
    fn dummy_plane_normal_direction() {
        let n = Vector3::new(0.0, 0.0, 1.0);
        let p = dummy_plane(n);
        // truck::Plane has a normal() method via ParametricSurface3D.
        let actual = ParametricSurface3D::normal(&p, 0.0, 0.0);
        // The normal must be parallel to n (up to sign).
        let parallel = actual.cross(n).magnitude() < 1e-9;
        assert!(parallel, "dummy_plane normal should be parallel to input");
    }
}
