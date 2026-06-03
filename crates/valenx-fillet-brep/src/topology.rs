//! BRep topology helpers: pull the inner truck Solid out of the
//! Valenx wrapper, get edge endpoints, count things.
//!
//! Everything here is borrow-only / read-only — no mutation, no
//! construction of new topology. Constructors live in [`crate::fillet`]
//! and [`crate::chamfer`].

use truck_modeling::{Edge as TruckEdge, Point3, Solid as TruckSolid};

use valenx_cad::Solid;

use crate::error::FilletBrepError;

/// Borrow the inner truck Solid out of a [`valenx_cad::Solid`].
///
/// Returns [`FilletBrepError::MeshBackedSolid`] if the input is a
/// `Solid::Mesh` (no BRep topology to operate on). This is the
/// canonical entry point for every BRep-fillet operation; the
/// dispatcher in `valenx-feature-tree` calls this first and falls
/// through to the mesh-domain pipeline on `Err`.
pub fn extract_brep(solid: &Solid) -> Result<&TruckSolid, FilletBrepError> {
    match solid {
        Solid::Brep(inner) => Ok(inner),
        Solid::Mesh(_) => Err(FilletBrepError::MeshBackedSolid),
    }
}

/// Return the two endpoints of a truck edge — `(front, back)` in the
/// edge's current orientation (so reversing the edge swaps them).
///
/// Pulls the [`Point3`] values out of the front/back vertices.
/// truck's [`truck_modeling::Vertex`] is `Vertex<Point3>` so the
/// extracted points are world-space coordinates.
pub fn edge_endpoints(edge: &TruckEdge) -> (Point3, Point3) {
    (edge.front().point(), edge.back().point())
}

/// Length of an edge measured by its endpoint distance.
///
/// **Caveat:** this is the *chord* length, not the geometric length of
/// the underlying curve. For straight edges (the only kind v1 fillets)
/// the two are identical; for curves a future revision would need to
/// integrate along the parametric curve.
pub fn edge_chord_length(edge: &TruckEdge) -> f64 {
    let (a, b) = edge_endpoints(edge);
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::primitives::box_solid;

    #[test]
    fn extract_brep_returns_inner_for_brep_variant() {
        let cube = box_solid(1.0, 1.0, 1.0).expect("unit cube builds");
        let brep = extract_brep(&cube).expect("brep variant returns Ok");
        // truck::Solid::boundaries() returns the shells.
        assert!(!brep.boundaries().is_empty(), "cube has at least one shell");
    }

    // Note: mesh-backed Solid construction requires valenx-mesh; the
    // corresponding `extract_brep_rejects_mesh_variant` test lives in
    // the feature-tree bridge tests where valenx-mesh is already a
    // transitive dependency. We don't pull valenx-mesh into
    // valenx-fillet-brep just to write one negative test — the
    // dispatch logic is what matters and is covered there.

    #[test]
    fn edge_endpoints_match_constructed_vertices() {
        // Build an edge from two known points and confirm endpoints
        // round-trip. Tests both the truck::builder integration and our
        // (front, back) ordering convention.
        use truck_modeling::{builder, InnerSpace, Point3};
        let a = builder::vertex(Point3::new(1.0, 2.0, 3.0));
        let b = builder::vertex(Point3::new(4.0, 5.0, 6.0));
        let edge = builder::line(&a, &b);
        let (front, back) = edge_endpoints(&edge);
        assert!((front - Point3::new(1.0, 2.0, 3.0)).magnitude() < 1e-9);
        assert!((back - Point3::new(4.0, 5.0, 6.0)).magnitude() < 1e-9);
    }

    #[test]
    fn edge_chord_length_unit_segment() {
        use truck_modeling::{builder, Point3};
        let a = builder::vertex(Point3::new(0.0, 0.0, 0.0));
        let b = builder::vertex(Point3::new(3.0, 4.0, 0.0));
        let edge = builder::line(&a, &b);
        // 3-4-5 triangle.
        assert!((edge_chord_length(&edge) - 5.0).abs() < 1e-9);
    }
}
