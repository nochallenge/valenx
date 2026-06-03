//! Phase 93 — `BRepOffsetAPI_DraftAngle` (draft / tapered faces).
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_DraftAngle(shape)` tilts selected faces of a solid
//! by a small angle (the "draft") with respect to a neutral plane —
//! the canonical operation for moulded / cast parts where every
//! vertical wall has to be slightly slanted so the part can be
//! ejected from the mould. Callers `Add(face, direction, angle,
//! neutral_plane)` per face; `Build()` warps the solid.
//!
//! Constraints:
//!
//! - The neutral plane is the plane along which the draft is measured.
//! - The draft direction is the mould-release direction.
//! - Angles are typically in the range `[1°, 5°]`; outside that range
//!   the result is rarely useful for manufacturing.
//!
//! ## v1 status
//!
//! Stub — true topology-preserving draft needs to redirect each face's
//! supporting surface (typically a `Geom_Plane`) and re-stitch the
//! edges; truck does not expose face-redirection primitives. Phase
//! 93.5 will ship the mesh-domain draft using
//! `valenx_blender_mesh_ops::shear`.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Apply a draft angle to the faces of `shape`.
///
/// `face_indices` identify which faces to draft (by their order in
/// `Solid::faces` traversal). `neutral_plane_z` is the world-space
/// Z height of the neutral plane; `direction` is the mould-release
/// direction; `angle_rad` is the draft.
///
/// # Errors
///
/// Always [`OcctSurfaceError::NotYetImplemented`] in v1.
pub fn offset_api_draft_angle(
    _shape: &Solid,
    face_indices: &[usize],
    direction: [f64; 3],
    angle_rad: f64,
    _neutral_plane_z: f64,
) -> Result<Solid, OcctSurfaceError> {
    if face_indices.is_empty() {
        return Err(OcctSurfaceError::bad_input(
            "face_indices",
            "no faces to draft",
        ));
    }
    let dir_len = (direction[0].powi(2) + direction[1].powi(2) + direction[2].powi(2)).sqrt();
    if dir_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "direction",
            "must be a non-zero vector",
        ));
    }
    if !angle_rad.is_finite() || angle_rad.abs() < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "angle_rad",
            format!("must be a non-zero finite angle, got {angle_rad}"),
        ));
    }
    Err(OcctSurfaceError::not_yet("offset_api_draft_angle"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn draft_rejects_empty_face_list() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_draft_angle(&cube, &[], [0.0, 0.0, 1.0], 0.05, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn draft_is_stub_with_valid_inputs() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_draft_angle(&cube, &[0], [0.0, 0.0, 1.0], 0.05, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.not_yet_implemented");
    }
}
