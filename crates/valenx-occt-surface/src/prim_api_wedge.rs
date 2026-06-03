//! Phase 86 — `BRepPrimAPI_MakeWedge`.
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeWedge(dx, dy, dz, ltx)` builds a wedge — a box
//! whose top face is shorter than its base along the X axis. The
//! parameters are:
//!
//! - `dx, dy, dz` — base footprint and total height (same as
//!   [`crate::prim_api_box()`]).
//! - `ltx` — length of the top edge in X (must satisfy `0 <= ltx <= dx`).
//!
//! The full 6-parameter overload `MakeWedge(dx, dy, dz, xmin, zmin,
//! xmax, zmax)` lets the caller specify the entire trapezoidal
//! cross-section of the top face explicitly.
//!
//! Most commonly used as the negative volume for sloped/chamfered
//! cuts (subtract a wedge from a slab to produce a tapered edge).
//!
//! ## v1 status
//!
//! **Honest implementation** using the 5-parameter simple wedge
//! constructor, built as a prism over a quadrilateral profile (the
//! base rectangle) extruded upward, then cut by a sloped plane —
//! v1 emits this as a prism over a wedge-cross-section polygon
//! using `valenx_cad::prism` after computing the cross section
//! analytically.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Build a wedge with base `dx × dy`, height `dz`, top X-length `ltx`.
///
/// The base face occupies `[0..dx] × [0..dy]` in the XY plane;
/// the top face is shrunk to `[0..ltx] × [0..dy]` at height `dz`,
/// producing a sloped X+ face.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] when `ltx > dx` or any dimension
/// is non-positive;
/// [`OcctSurfaceError::TruckLimit`] when the prism builder fails.
pub fn prim_api_wedge(dx: f64, dy: f64, dz: f64, ltx: f64) -> Result<Solid, OcctSurfaceError> {
    for (name, v) in [("dx", dx), ("dy", dy), ("dz", dz)] {
        if !v.is_finite() || v <= 0.0 {
            return Err(OcctSurfaceError::bad_input(
                "dimension",
                format!("{name} must be positive finite, got {v}"),
            ));
        }
    }
    if !ltx.is_finite() || ltx < 0.0 || ltx > dx {
        return Err(OcctSurfaceError::bad_input(
            "ltx",
            format!("must be in [0, dx={dx}], got {ltx}"),
        ));
    }
    // v1 strategy: take the y=const cross-section view from the side.
    // It's a quadrilateral with vertices:
    //   (0, 0), (dx, 0), (ltx, dz), (0, dz)
    // Extruding this XZ-quad along Y for `dy` would yield the wedge.
    // But `valenx_cad::prism` extrudes an XY profile along +Z, so we
    // build the quad in the XY plane and extrude along Y... except
    // prism only goes +Z. So we extrude the side-view quad (now sitting
    // in the XZ plane mapped to "XY" of the prism API) along +Z by dy
    // and re-interpret axes after the fact.
    //
    // After construction the solid has:
    //   X = original X (base length)
    //   Y = original Z (height) — mismatched with caller expectation
    //   Z = original Y (depth)
    //
    // To match the OCCT convention we rotate the result 90° about the
    // X axis so the height ends up along +Z. That's a clone-and-
    // transform; cheap.
    let quad_xz = [(0.0, 0.0), (dx, 0.0), (ltx, dz), (0.0, dz)];
    let solid = valenx_cad::prism(&quad_xz, dy)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("wedge prism: {e:?}")))?;
    // Rotate 90° about +X to swap Y↔Z.
    let rotated = solid
        .rotated((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), std::f64::consts::FRAC_PI_2)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("wedge rotate: {e}")))?;
    Ok(rotated)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wedge_builds_with_six_faces() {
        let w = prim_api_wedge(2.0, 1.0, 1.0, 1.0).unwrap();
        // A wedge has 6 faces (base + top + 4 sides), same as a box.
        assert_eq!(w.faces(), 6);
        assert_eq!(w.vertices(), 8);
    }

    #[test]
    fn wedge_rejects_oversized_top() {
        let err = prim_api_wedge(1.0, 1.0, 1.0, 2.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn wedge_rejects_zero_dimension() {
        let err = prim_api_wedge(0.0, 1.0, 1.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }
}
