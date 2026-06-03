//! BRep ↔ polyhedron conversion.

use valenx_cad::Solid;
use valenx_mesh::Mesh;

use crate::error::MeshPartError;

/// Tessellate a [`Solid`] into a polyhedral [`Mesh`] with an
/// adjustable chord-error budget. Thin shim over
/// [`valenx_cad::solid_to_mesh`] so MeshPart-flavoured callers don't
/// have to depend on `valenx-cad` directly.
///
/// `tolerance` is in the same units as the model (mm by default in
/// Valenx).
pub fn brep_to_polyhedron(solid: &Solid, tolerance: f64) -> Result<Mesh, MeshPartError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(MeshPartError::BadParameter {
            name: "tolerance",
            reason: format!("must be finite > 0, got {tolerance}"),
        });
    }
    valenx_cad::solid_to_mesh(solid, tolerance).map_err(|e| MeshPartError::Cad(e.to_string()))
}

/// Naïve polyhedron → BRep promotion. Each triangle becomes one
/// planar face; the implementation tries to sew the resulting shell
/// into a true [`Solid::Brep`]. **v1 limitation:** truck-modeling
/// 0.6 doesn't expose the planar-face-from-triangle constructor we'd
/// need here, so v1 always returns a mesh-backed [`Solid::Mesh`]
/// wrapped via [`Solid::from_mesh`]. The shape of the API stays
/// stable so Phase 32.5 can slot in the real sewing without touching
/// callers.
///
/// The function never errors — sewing failures return
/// [`MeshPartError::SewingFallback`] **only** when a future
/// implementation actually attempts BRep sewing and fails; v1 always
/// returns `Ok(Solid::Mesh)`.
pub fn polyhedron_to_brep(mesh: &Mesh) -> Result<Solid, MeshPartError> {
    if mesh.nodes.is_empty() {
        return Err(MeshPartError::Empty("mesh"));
    }
    // v1 keep-it-honest: wrap as Solid::Mesh.
    Ok(Solid::from_mesh(mesh.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polyhedron_to_brep_empty_errors() {
        let m = Mesh::new("empty");
        assert!(matches!(
            polyhedron_to_brep(&m),
            Err(MeshPartError::Empty("mesh"))
        ));
    }

    #[test]
    fn polyhedron_to_brep_round_trips_via_mesh() {
        let mut m = Mesh::new("t");
        m.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        let s = polyhedron_to_brep(&m).unwrap();
        // v1 must currently return the mesh-backed variant.
        assert!(matches!(s, Solid::Mesh(_)));
    }
}
