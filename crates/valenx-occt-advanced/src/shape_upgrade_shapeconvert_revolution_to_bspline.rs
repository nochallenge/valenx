//! Phase 148 — `ShapeUpgrade_ConvertSurfaceToBezierBasis` for
//! revolution surfaces.
//!
//! ## What OCCT does
//!
//! `ShapeCustom_DirectModification(shape, ConvertToBSplineSurface)`
//! walks every face and replaces analytic surface types
//! (`Geom_SurfaceOfRevolution`, `Geom_CylindricalSurface`,
//! `Geom_ConicalSurface`, `Geom_SphericalSurface`,
//! `Geom_ToroidalSurface`) with an equivalent
//! `Geom_BSplineSurface`. Required for downstream tools that only
//! accept B-spline input (e.g. some FEA mesh generators that don't
//! know how to handle parametric tori).
//!
//! The conversion is exact for cylinders / cones / spheres (rational
//! NURBS can represent quadrics exactly with degree 2 + 3 weights);
//! tori need degree 4 × 4. Revolution surfaces are converted by
//! sampling the meridian curve + the rotation rail and tensoring.
//!
//! ## v1 status
//!
//! Stub — needs face-level surface introspection to detect a
//! `SurfaceOfRevolution` instance, then an exact-conversion algorithm
//! to emit the equivalent NURBS. valenx-surface has the NURBS
//! representation but no analytic-surface-to-NURBS converter. Phase
//! 148.5 ships with the converter (which itself ports the classic
//! "Piegl & Tiller §7.5" algorithm).

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Convert every revolution surface in `solid` to an equivalent
/// B-spline surface.
///
/// # Errors
///
/// Always [`OcctAdvancedError::NotYetImplemented`] in v1.
pub fn shape_upgrade_shapeconvert_revolution_to_bspline(
    _solid: &Solid,
) -> Result<Solid, OcctAdvancedError> {
    Err(OcctAdvancedError::not_yet(
        "shape_upgrade_shapeconvert_revolution_to_bspline",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn stub_with_cube_input() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_shapeconvert_revolution_to_bspline(&cube).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
    }
}
