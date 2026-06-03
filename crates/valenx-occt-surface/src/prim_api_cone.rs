//! Phase 80 — `BRepPrimAPI_MakeCone` (pointed cone or truncated
//! frustum).
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeCone(R1, R2, H)` builds a closed cone whose base
//! disk has radius `R1`, top disk has radius `R2` (zero for a
//! pointed apex), and height `H` along +Z. There's also the
//! `MakeCone(Ax2, R1, R2, H, Angle)` overload for an oriented cone
//! sector — same axis story as cylinder.
//!
//! ## v1 status
//!
//! **Honest implementation** for the axis-aligned case via
//! `valenx_cad::cone`. Set `top_radius = 0.0` for the pointed cone.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Closed cone (or frustum if `top_radius > 0`) along +Z.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] for invalid dimensions.
pub fn prim_api_cone(
    base_radius: f64,
    top_radius: f64,
    height: f64,
) -> Result<Solid, OcctSurfaceError> {
    valenx_cad::cone(base_radius, top_radius, height)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("cone: {e:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pointed_cone_builds() {
        let c = prim_api_cone(1.0, 0.0, 2.0).unwrap();
        assert!(c.faces() > 0);
    }

    #[test]
    fn frustum_builds() {
        let c = prim_api_cone(2.0, 1.0, 3.0).unwrap();
        assert!(c.faces() > 0);
    }

    #[test]
    fn cone_rejects_negative_top_radius() {
        let err = prim_api_cone(1.0, -0.5, 2.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }
}
