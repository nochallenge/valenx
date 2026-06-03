//! Phase 158 — `BRepFeat_Builder::ReplaceFace` — replace one face of
//! a solid with a new face geometry.
//!
//! ## What OCCT does
//!
//! `BRepFeat_Builder` is the local-operation API: instead of
//! rebuilding the entire feature tree, it surgically modifies a
//! single face/edge/vertex while preserving the surrounding topology.
//! `ReplaceFace(old_face, new_face)` consumes a face from the input
//! solid and a replacement face (typically a re-fitted surface), then:
//!
//! 1. Identifies all edges shared between `old_face` and its
//!    neighbors.
//! 2. Re-projects those edges onto `new_face`'s surface (computing
//!    new PCurves).
//! 3. Stitches `new_face` into the original solid's topology in
//!    place of `old_face`.
//!
//! Used by FreeCAD's "DraftAngle" tool, SolidWorks' "Move Face", and
//! generally anywhere the user wants to tweak a single face without
//! rolling back the feature history.
//!
//! ## v1 status
//!
//! Stub — the in-place face-replacement requires the same face-
//! mutation infrastructure as the `shape_upgrade_*` modules.
//! truck-modeling exposes face *reading* but not face *replacement*.
//! Phase 158.5 ships once valenx-cad lands a `replace_face` op
//! (likely built on top of an export-rebuild round-trip until truck
//! gains a native operator).

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Replace face `face_index` of `solid` with a new face built from
/// `new_face_vertices` (a closed polyline defining the new face's
/// outer wire).
///
/// `face_index` is 0-based into the solid's face iterator order
/// (see [`valenx_cad::Solid::faces`]).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for out-of-range `face_index`,
///   too few new-face vertices.
/// - [`OcctAdvancedError::NotYetImplemented`] otherwise in v1.
pub fn local_op_replace_face(
    solid: &Solid,
    face_index: usize,
    new_face_vertices: &[[f64; 3]],
) -> Result<Solid, OcctAdvancedError> {
    if new_face_vertices.len() < 3 {
        return Err(OcctAdvancedError::bad_input(
            "new_face_vertices",
            format!(
                "need ≥3 vertices for a face boundary; got {}",
                new_face_vertices.len()
            ),
        ));
    }
    let face_count = solid.faces();
    if face_count > 0 && face_index >= face_count {
        return Err(OcctAdvancedError::bad_input(
            "face_index",
            format!("{face_index} out of range (solid has {face_count} faces)"),
        ));
    }
    Err(OcctAdvancedError::not_yet("local_op_replace_face"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    fn square() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    }

    #[test]
    fn rejects_short_face() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = local_op_replace_face(&cube, 0, &[[0.0; 3], [1.0, 0.0, 0.0]]).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_out_of_range_index() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        // Cube has 6 faces — index 6 is out of range.
        let err = local_op_replace_face(&cube, 6, &square()).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn stub_with_valid_inputs() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = local_op_replace_face(&cube, 0, &square()).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
    }
}
