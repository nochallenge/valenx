//! Phase 159 — `BRepFeat_Builder::RemoveFace` — remove a face from a
//! solid (creates an open shell).
//!
//! ## What OCCT does
//!
//! Discards `face` from the solid's boundary representation. The
//! result is no longer closed — it's a "shell" with one fewer face,
//! which downstream ops can sew shut, extrude into a solid, or use
//! as a surface for FEM analysis. Used for "split this face away
//! from the body" workflows in CAD UI ("Detach face" in FreeCAD).
//!
//! Importantly: removing a face *keeps the edges and vertices* that
//! the face referenced — they survive as free boundary edges of the
//! remaining shell (see [`crate::shape_analysis_freebounds()`]).
//!
//! ## v1 status
//!
//! Stub — same topology-mutation gap as
//! [`crate::local_op_replace_face()`]. Phase 159.5 ships with Phase
//! 158.5. v1 callers can pre-mesh the solid via
//! [`valenx_cad::solid_to_mesh`], drop the chosen face's triangles
//! by hand, and re-emit as a mesh-backed `Solid::Mesh` — that's the
//! visualisation-grade fallback.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Remove face `face_index` from `solid` (leaves an open shell).
///
/// `face_index` is 0-based into the solid's face iterator order
/// (see [`valenx_cad::Solid::faces`]).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for out-of-range `face_index`.
/// - [`OcctAdvancedError::NotYetImplemented`] otherwise in v1.
pub fn local_op_remove_face(solid: &Solid, face_index: usize) -> Result<Solid, OcctAdvancedError> {
    let face_count = solid.faces();
    if face_count > 0 && face_index >= face_count {
        return Err(OcctAdvancedError::bad_input(
            "face_index",
            format!("{face_index} out of range (solid has {face_count} faces)"),
        ));
    }
    if face_count == 1 {
        return Err(OcctAdvancedError::bad_input(
            "solid",
            "cannot remove the only face of a single-face shell",
        ));
    }
    Err(OcctAdvancedError::not_yet("local_op_remove_face"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn rejects_out_of_range_index() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = local_op_remove_face(&cube, 6).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn stub_with_valid_index() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = local_op_remove_face(&cube, 0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
    }
}
