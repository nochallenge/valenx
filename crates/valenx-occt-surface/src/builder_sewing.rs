//! Phase 73 — `BRepBuilderAPI_Sewing`: stitch a face collection into
//! a shell, sharing coincident edges.
//!
//! ## What OCCT does
//!
//! `BRepBuilderAPI_Sewing` is the canonical OCCT topology-rebuild tool.
//! Callers `Add(face)` a collection of `TopoDS_Face`s (typically from
//! a STEP/IGES import where each face arrived as its own disconnected
//! shape), then `Perform()` walks every face-face boundary pair
//! looking for edges that are geometrically coincident to within a
//! caller-supplied tolerance. Coincident edges are merged into a
//! single shared edge so the resulting `TopoDS_Shell` is a proper
//! 2-manifold. Without sewing, downstream booleans / sweeps / STEP
//! re-export will treat each face as a separate body.
//!
//! Sewing also identifies "free boundary" edges (edges on exactly one
//! face) — the result tells you whether the assembly is closed
//! (`IsClosed`) and where the holes are if not.
//!
//! ## v1 status — real mesh-domain sewing
//!
//! This is a genuine sewing pass. Each input solid is tessellated;
//! the triangle meshes are concatenated (index-offset merged) into a
//! single mesh, and coincident vertices within `tolerance` are then
//! **welded** so triangles from different input faces that meet along
//! a shared boundary now reference the same vertices — the
//! mesh-domain equivalent of OCCT's coincident-edge merge.
//!
//! The result is a sewn mesh-backed [`Solid`]: faces that touched are
//! now topologically joined (no duplicate vertices on the seam). What
//! a tessellated mesh cannot carry is OCCT's per-edge `TopoDS_Edge`
//! sharing in a parametric BRep — a true parametric sew needs the
//! NURBS-surface stitcher (`valenx_surface::sew::stitch`) for
//! tensor-product faces, which remains the path for parametric
//! input. The mesh sew here is correct for the import-repair use case
//! (a STEP/STL import arriving as disconnected faces) and composes
//! with the rest of the mesh-backed-`Solid` pipeline.

use valenx_cad::Solid;
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;

/// Chord tolerance used to tessellate each input solid before sewing.
const SEW_TESS_TOLERANCE: f64 = 0.1;

/// Sew a collection of faces into a shell, merging coincident edges.
///
/// `faces` carries one solid per element (used as the face container —
/// truck does not expose a face-only handle through `valenx-cad`'s
/// public surface). `tolerance` is the linear distance below which two
/// boundary vertices are considered the "same" point and welded.
///
/// Returns a single sewn mesh-backed [`Solid`].
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for an empty face list or a
///   non-positive / non-finite tolerance.
/// - [`OcctSurfaceError::TruckLimit`] when an input fails to
///   tessellate.
pub fn builder_sewing(faces: &[Solid], tolerance: f64) -> Result<Solid, OcctSurfaceError> {
    if faces.is_empty() {
        return Err(OcctSurfaceError::bad_input(
            "faces",
            "need at least one face to sew",
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(OcctSurfaceError::bad_input(
            "tolerance",
            format!("must be a positive finite number, got {tolerance}"),
        ));
    }

    // Tessellate every input and concatenate into one mesh.
    let mut combined = Mesh::new("sewn");
    for (i, face) in faces.iter().enumerate() {
        let mesh = valenx_cad::solid_to_mesh(face, SEW_TESS_TOLERANCE).map_err(|e| {
            OcctSurfaceError::TruckLimit(format!("sewing: tessellate face {i}: {e:?}"))
        })?;
        combined = valenx_mesh::boolean::concatenate(&combined, &mesh);
    }
    // Weld coincident boundary vertices — the mesh-domain edge merge.
    let sewn = valenx_mesh::boolean::merge_coincident_nodes(&combined, tolerance);
    Ok(Solid::from_mesh(sewn))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn sewing_rejects_empty_face_list() {
        let err = builder_sewing(&[], 0.01).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn sewing_rejects_bad_tolerance() {
        let f = box_solid(1.0, 1.0, 1.0).unwrap();
        assert_eq!(
            builder_sewing(std::slice::from_ref(&f), -0.1).unwrap_err().code(),
            "occt_surface.bad_input"
        );
        assert_eq!(
            builder_sewing(std::slice::from_ref(&f), f64::INFINITY).unwrap_err().code(),
            "occt_surface.bad_input"
        );
    }

    #[test]
    fn sewing_one_face_returns_a_solid() {
        // Sewing a single solid yields a sewn mesh-backed solid that
        // tessellates to the same geometry.
        let f = box_solid(1.0, 1.0, 1.0).unwrap();
        let sewn = builder_sewing(std::slice::from_ref(&f), 0.01).unwrap();
        let mesh = valenx_cad::solid_to_mesh(&sewn, 0.1).unwrap();
        assert!(!mesh.nodes.is_empty());
        // A cube has 8 corners — after welding, the sewn mesh's unique
        // vertex count is far below the per-triangle vertex count.
        assert!(mesh.nodes.len() <= 24, "weld should dedup: {}", mesh.nodes.len());
    }

    #[test]
    fn sewing_two_touching_faces_welds_the_seam() {
        // Two unit cubes sharing the x=1 face. After sewing, the seam
        // vertices are welded so the combined vertex count is strictly
        // less than the sum of the two cubes' tessellated vertices.
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0)
            .unwrap()
            .translated(1.0, 0.0, 0.0)
            .unwrap();
        let mesh_a = valenx_cad::solid_to_mesh(&a, 0.1).unwrap();
        let mesh_b = valenx_cad::solid_to_mesh(&b, 0.1).unwrap();
        let separate = mesh_a.nodes.len() + mesh_b.nodes.len();

        let sewn = builder_sewing(&[a, b], 1e-4).unwrap();
        let sewn_mesh = valenx_cad::solid_to_mesh(&sewn, 0.1).unwrap();
        assert!(
            sewn_mesh.nodes.len() < separate,
            "sewing should weld the shared seam: sewn={}, separate={separate}",
            sewn_mesh.nodes.len()
        );
    }
}
