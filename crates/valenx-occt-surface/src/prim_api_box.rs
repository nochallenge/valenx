//! Phase 83 — `BRepPrimAPI_MakeBox` (axis-aligned + oriented).
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeBox(dx, dy, dz)` builds a rectangular solid with
//! one corner at the origin and the opposite corner at
//! `(dx, dy, dz)`. Overloads:
//!
//! - `MakeBox(P, dx, dy, dz)` — anchored at `gp_Pnt P`.
//! - `MakeBox(P1, P2)` — diagonal corners.
//! - `MakeBox(Ax2, dx, dy, dz)` — oriented along the user's coord
//!   system.
//!
//! ## v1 status
//!
//! **Honest implementation** for the axis-aligned origin-anchored
//! case via `valenx_cad::box_solid`. The oriented overload is
//! achievable by post-multiplying with `Solid::rotated` /
//! `Solid::translated` but is deferred to Phase 83.5 for a single
//! parameterised constructor.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Axis-aligned box with one corner at the origin.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] for non-positive dimensions.
pub fn prim_api_box(dx: f64, dy: f64, dz: f64) -> Result<Solid, OcctSurfaceError> {
    valenx_cad::box_solid(dx, dy, dz)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("box: {e:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_builds_with_expected_topology() {
        let b = prim_api_box(1.0, 2.0, 3.0).unwrap();
        assert_eq!(b.faces(), 6);
        assert_eq!(b.edges(), 12);
        assert_eq!(b.vertices(), 8);
    }

    #[test]
    fn box_rejects_zero_dimension() {
        let err = prim_api_box(0.0, 1.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }
}
