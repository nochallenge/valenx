//! Phase 84 — `BRepPrimAPI_MakePrism` (linear extrusion of a profile).
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakePrism(profile, vec)` extrudes a `TopoDS_Wire` or
//! `TopoDS_Face` by translation along `gp_Vec`. The result is a
//! `TopoDS_Solid` if the profile was a planar face. Overloads:
//!
//! - `MakePrism(profile, vec, Copy=true)` — clone the profile (default).
//! - `MakePrism(profile, vec, Copy=false)` — consume the profile.
//! - `MakePrism(profile, dir, height)` — direction + scalar length.
//!
//! ## v1 status
//!
//! **Honest implementation** for a closed planar polygon profile via
//! `valenx_cad::prism`. The profile must lie in the XY plane and the
//! extrusion is along +Z by `height`. General "any-direction" prism
//! of an arbitrary BRep face is deferred to Phase 84.5.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Extrude a closed XY polygon along +Z.
///
/// `profile_xy` is the vertex list (not auto-closed; the last point
/// should differ from the first); `height` is the extrusion length.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] when the underlying
/// `valenx-cad` builder rejects the inputs.
pub fn prim_api_prism(profile_xy: &[(f64, f64)], height: f64) -> Result<Solid, OcctSurfaceError> {
    valenx_cad::prism(profile_xy, height)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("prism: {e:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangular_prism_builds() {
        let tri = prim_api_prism(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)], 2.0).unwrap();
        // 2 triangular caps + 3 rectangular sides = 5 faces.
        assert_eq!(tri.faces(), 5);
        assert_eq!(tri.vertices(), 6);
    }

    #[test]
    fn prism_rejects_short_profile() {
        let err = prim_api_prism(&[(0.0, 0.0), (1.0, 0.0)], 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }
}
