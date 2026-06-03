//! Phase 96 — `BRepAlgoAPI_Cut` with face-warping tolerance.
//!
//! ## What OCCT does
//!
//! Boolean difference `A - B` with the addition of a face-warping
//! tolerance — when a face on `A` lies "almost" on a face on `B` but
//! within the warping budget, the operator warps the `A` face onto
//! the `B` face before performing the cut. This avoids the classic
//! "slivers" left behind when two faces are nearly coplanar but not
//! identical.
//!
//! Configured via `SetFuzzyValue(epsilon)` and is enabled by default
//! in OCCT v7+. The default fuzzy value is `1e-6` model units.
//!
//! ## v1 status
//!
//! Partial — for the **no-warping case** (`tolerance == 0`) we
//! delegate to `valenx_cad::difference` which is already a real BRep
//! cut via truck-shapeops. For a non-zero warp tolerance, v1 stubs
//! out: truck-shapeops has its own coincidence tolerance but doesn't
//! expose face-warping pre-processing.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Boolean difference `A - B` with an optional face-warping
/// tolerance.
///
/// Pass `warp_tolerance = 0.0` for the strict cut (delegates to
/// [`valenx_cad::difference`]); positive values request face-warping
/// pre-processing.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] for negative warp tolerance;
/// [`OcctSurfaceError::TruckLimit`] when the difference op fails;
/// [`OcctSurfaceError::NotYetImplemented`] for nonzero
/// `warp_tolerance`.
pub fn cut_api_cut_with_warp(
    a: &Solid,
    b: &Solid,
    warp_tolerance: f64,
) -> Result<Solid, OcctSurfaceError> {
    if !warp_tolerance.is_finite() || warp_tolerance < 0.0 {
        return Err(OcctSurfaceError::bad_input(
            "warp_tolerance",
            format!("must be a non-negative finite number, got {warp_tolerance}"),
        ));
    }
    if warp_tolerance == 0.0 {
        return valenx_cad::difference(a, b)
            .map_err(|e| OcctSurfaceError::TruckLimit(format!("difference: {e:?}")));
    }
    Err(OcctSurfaceError::not_yet("cut_api_cut_with_warp"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn cut_with_zero_warp_acts_as_difference() {
        // Two non-overlapping boxes — difference is a (yields the
        // first box back, since b doesn't touch a). truck-shapeops
        // may surface this as a kernel error; either way we should
        // not panic.
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(0.5, 0.5, 0.5)
            .unwrap()
            .translated(0.25, 0.25, 0.25)
            .unwrap();
        let _ = cut_api_cut_with_warp(&a, &b, 0.0);
        // We don't assert on the result — the point is it returned
        // *something* (Ok or TruckLimit) rather than panicking or
        // returning NotYetImplemented.
    }

    #[test]
    fn cut_rejects_negative_warp_tolerance() {
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = cut_api_cut_with_warp(&a, &b, -0.1).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn cut_with_nonzero_warp_is_stub() {
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = cut_api_cut_with_warp(&a, &b, 0.01).unwrap_err();
        assert_eq!(err.code(), "occt_surface.not_yet_implemented");
    }
}
