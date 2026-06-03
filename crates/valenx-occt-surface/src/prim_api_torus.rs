//! Phase 82 — `BRepPrimAPI_MakeTorus`.
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeTorus(R1, R2)` builds a closed torus with
//! major-circle radius `R1` and tube radius `R2` (must satisfy
//! `R2 < R1` to be self-non-intersecting). Overloads exist for
//! oriented `Ax2` placement and for partial tori (sweep angle less
//! than 2π).
//!
//! ## v1 status
//!
//! **Honest implementation** for the full torus via
//! `valenx_cad::torus`. Partial tori are deferred.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Full torus centred on the origin, major axis along +Z.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] when `minor_radius >= major_radius`
/// (self-intersecting torus) or for invalid dimensions.
pub fn prim_api_torus(major_radius: f64, minor_radius: f64) -> Result<Solid, OcctSurfaceError> {
    valenx_cad::torus(major_radius, minor_radius)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("torus: {e:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn torus_builds() {
        let t = prim_api_torus(2.0, 0.5).unwrap();
        assert!(t.faces() > 0);
    }

    #[test]
    fn torus_rejects_self_intersecting() {
        let err = prim_api_torus(1.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }
}
