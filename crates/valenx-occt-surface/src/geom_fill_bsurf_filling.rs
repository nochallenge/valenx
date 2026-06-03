//! Phase 75 — `GeomFill_BSplineCurves`: fill a quad of boundary
//! curves with a G1 patch (Coons-like).
//!
//! ## What OCCT does
//!
//! `GeomFill_BSplineCurves` builds a `Geom_BSplineSurface` that
//! interpolates four input `Geom_BSplineCurve` boundary curves, with
//! tangent continuity at the corners. The fill style is one of:
//!
//! - `GeomFill_StretchStyle` — minimum-area patch (default for thin
//!   strips).
//! - `GeomFill_CoonsStyle` — classic bilinearly-blended Coons patch,
//!   exact at the boundary, free interior.
//! - `GeomFill_CurvedStyle` — third-order interior shape minimising a
//!   curvature integral.
//!
//! Used as the workhorse for filling small holes left by trim
//! operations or for stitching adjacent surface patches with G1
//! continuity at the seam.
//!
//! ## v1 status
//!
//! **Honest implementation.** `valenx-surface::coons::fill` already
//! delivers the Coons-style variant of this exact operation: four
//! boundary NURBS curves → one tensor-product NURBS surface that
//! interpolates them. We delegate to it.

use valenx_surface::{coons, NurbsCurve, NurbsSurface};

use crate::error::OcctSurfaceError;

/// Fill four boundary curves with a Coons patch.
///
/// The four boundaries form a topological quad: `c0` and `c1` are
/// "opposite" sides (running along the `v` parameter), `d0` and `d1`
/// are the other pair (running along `u`). All four must be cubic
/// NURBS curves sharing endpoints at the patch corners — matching
/// `valenx_surface::coons::fill` boundary contract.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] wrapping the underlying
/// `SurfaceError` when the four curves don't share corners or have
/// incompatible degree.
pub fn geom_fill_bsurf_filling(
    c0: NurbsCurve,
    c1: NurbsCurve,
    d0: NurbsCurve,
    d1: NurbsCurve,
) -> Result<NurbsSurface, OcctSurfaceError> {
    coons::fill([c0, c1, d0, d1])
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("coons::fill: {e:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// Cubic-Bezier straight line from `a` to `b` with intermediate
    /// CPs at 1/3 and 2/3 — matching valenx-surface's own coons test
    /// helper.
    fn line(a: Vector3<f64>, b: Vector3<f64>) -> NurbsCurve {
        let p1 = a + (b - a) / 3.0;
        let p2 = a + 2.0 * (b - a) / 3.0;
        NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![a, p1, p2, b],
            vec![1.0; 4],
        )
        .unwrap()
    }

    #[test]
    fn fill_unit_square_returns_cubic_surface() {
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p01 = Vector3::new(0.0, 1.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);

        let c0 = line(p00, p01);
        let c1 = line(p10, p11);
        let d0 = line(p00, p10);
        let d1 = line(p01, p11);

        let surf = geom_fill_bsurf_filling(c0, c1, d0, d1)
            .expect("Coons fill on a unit square should succeed");
        assert_eq!(surf.u_degree, 3);
        assert_eq!(surf.v_degree, 3);
    }
}
