//! Phase 152 — `ShapeUpgrade_FaceDivide` (geometric) — divide a face
//! into multiple smaller faces.
//!
//! ## What OCCT does
//!
//! `ShapeUpgrade_FaceDivide(face)` parameterises the face's surface
//! into an `n_u × n_v` grid of iso-parametric strips and emits each
//! strip as a separate face. The result is topologically equivalent
//! to the input face (same surface, same boundary) but split into
//! quads for downstream meshing or import targets that can't handle
//! arbitrarily-large faces.
//!
//! Distinct from [`crate::shape_upgrade_split_continuity()`] (Phase
//! 150) which only splits at discontinuities — this op always splits,
//! at a caller-chosen grid resolution.
//!
//! ## v1 status
//!
//! Stub — requires the same face mutation infrastructure as the
//! other `shape_upgrade_*` modules. truck-modeling exposes none of
//! it. Phase 152.5 ships once Phase 149.5's face-rebuild API lands.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Apply geometric face division to `solid`.
///
/// `face_index` — 0-based into the solid's face iterator order (see
/// [`valenx_cad::Solid::faces`]).
/// `n_u`, `n_v` — grid resolution (number of strips per parametric
/// dimension; total subfaces = `n_u * n_v`).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for `n_u < 1` or `n_v < 1`.
/// - [`OcctAdvancedError::NotYetImplemented`] otherwise in v1.
pub fn shape_upgrade_face_division(
    _solid: &Solid,
    _face_index: usize,
    n_u: usize,
    n_v: usize,
) -> Result<Solid, OcctAdvancedError> {
    if n_u < 1 || n_v < 1 {
        return Err(OcctAdvancedError::bad_input(
            "n_u/n_v",
            "grid resolution must be ≥1 per dimension",
        ));
    }
    if n_u == 1 && n_v == 1 {
        return Err(OcctAdvancedError::bad_input(
            "n_u/n_v",
            "1×1 grid is a no-op; pick at least 2 in one dimension",
        ));
    }
    Err(OcctAdvancedError::not_yet("shape_upgrade_face_division"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn rejects_zero_grid() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_face_division(&cube, 0, 0, 2).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_unit_grid() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_face_division(&cube, 0, 1, 1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn stub_with_valid_input() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_face_division(&cube, 0, 2, 2).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
    }
}
