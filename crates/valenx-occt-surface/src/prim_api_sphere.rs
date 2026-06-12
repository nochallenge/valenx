//! Phase 81 — `BRepPrimAPI_MakeSphere` (full or restricted-angle).
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeSphere(R)` builds the full sphere; overloads
//! `MakeSphere(R, Angle1, Angle2, Angle3)` produce a sphere segment
//! restricted in latitude (`Angle1`/`Angle2`) and longitude
//! (`Angle3`). The restricted variant is used for hemispheres, orange
//! slices, and the canonical "sphere with one face removed" test
//! shape from the OCCT samples.
//!
//! ## v1 status
//!
//! **Honest implementation** for the full sphere via
//! `valenx_cad::sphere`. The angle-restricted variants would need to
//! revolve a partial arc and close the polar caps — deferred to
//! Phase 81.5 (the math is straightforward but the polar
//! degeneracy needs careful BRep wiring).

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Full sphere centred on the origin.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] for non-positive radius.
pub fn prim_api_sphere(radius: f64) -> Result<Solid, OcctSurfaceError> {
    valenx_cad::sphere(radius).map_err(|e| OcctSurfaceError::TruckLimit(format!("sphere: {e:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sphere_builds() {
        let s = prim_api_sphere(1.0).unwrap();
        assert!(s.faces() > 0);
    }

    #[test]
    fn sphere_rejects_zero_radius() {
        let err = prim_api_sphere(0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }
}
