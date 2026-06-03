//! Phase 94 — `BRepOffsetAPI_MakeFilling` (N-sided patch fill).
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_MakeFilling()` builds an N-sided patch (any number
//! of boundary curves, not restricted to four) that smoothly fills a
//! hole or stitches together adjacent face boundaries. The algorithm
//! uses a least-squares fit minimising a curvature integral subject
//! to:
//!
//! - `Add(edge, GeomAbs_Shape continuity)` — boundary edges with
//!   their continuity requirement (G0 / G1 / G2 / G3) at the seam.
//! - `Add(point)` — interior points the patch must pass through.
//! - `Add(point, normal)` — interior points the patch must pass
//!   through with a tangent constraint.
//!
//! Unlike [`crate::geom_fill_bsurf_filling()`] (Coons, exactly four
//! sides), this supports arbitrary N including pentagons, triangles,
//! and "filling-the-gap" between four mostly-adjacent face edges. It
//! is the go-to tool for closing up gaps in IGES imports.
//!
//! ## v1 status — real 3- and 4-sided patch fill
//!
//! The single-NURBS-surface return type bounds what is honestly
//! representable: a tensor-product NURBS surface is intrinsically a
//! **4-sided** patch.
//!
//! - **4 boundary curves** — built as a real bilinearly-blended
//!   **Coons patch** via [`valenx_surface::coons::fill`]. The four
//!   curves are taken **in Coons order** `[c0, c1, d0, d1]`: `c0` /
//!   `c1` are the opposite pair along `v = min` / `v = max`, `d0` /
//!   `d1` the opposite pair along `u = min` / `u = max` (this is the
//!   ordering `valenx_surface::coons` documents — not a sequential
//!   edge loop). The patch interpolates all four edges.
//! - **3 boundary curves** — built as a **degenerate Coons patch**:
//!   the `c1` side is collapsed to the shared corner of `d0` and
//!   `d1`, the standard technique for a triangular Coons fill. The
//!   surface interpolates the three real edges; the collapsed side is
//!   a pole. Input order: `[d0, d1, c0]` (the two "rail" edges then
//!   the base edge).
//! - **≥ 5 boundary curves** — a genuine N-sided patch (`N ≥ 5`)
//!   cannot be expressed as one tensor-product NURBS surface without
//!   fabricated structure. Rather than return a fake, this case is a
//!   clear [`OcctSurfaceError::BadInput`] — the caller should split
//!   the N-gon into 4-sided sub-patches or use a mesh-domain fill.
//!   (OCCT's own `MakeFilling` returns a *trimmed* surface for these,
//!   which `valenx`'s untrimmed `NurbsSurface` cannot carry.)
//!
//! `interior_points` are accepted for API parity; the Coons blend is
//! fully determined by its boundary, so they currently act only as a
//! (future) refinement hook.

use valenx_surface::coons::fill;
use valenx_surface::{NurbsCurve, NurbsSurface};

use crate::error::OcctSurfaceError;

/// Fill a 3- or 4-sided hole bounded by `boundary` curves.
///
/// See the module docs for the boundary-curve ordering convention.
/// `interior_points` are optional constraint points (API parity).
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] when the boundary has fewer than 3
/// curves, more than 4 (not representable as one NURBS surface), or
/// the curves do not meet at their corners.
pub fn offset_api_filling(
    boundary: &[NurbsCurve],
    _interior_points: &[[f64; 3]],
) -> Result<NurbsSurface, OcctSurfaceError> {
    match boundary.len() {
        0..=2 => Err(OcctSurfaceError::bad_input(
            "boundary",
            format!("need at least 3 boundary curves, got {}", boundary.len()),
        )),
        3 => fill_triangular(&boundary[0], &boundary[1], &boundary[2]),
        4 => fill([
            boundary[0].clone(),
            boundary[1].clone(),
            boundary[2].clone(),
            boundary[3].clone(),
        ])
        .map_err(|e| OcctSurfaceError::bad_input("boundary", format!("coons fill: {e}"))),
        n => Err(OcctSurfaceError::bad_input(
            "boundary",
            format!(
                "an {n}-sided patch (N≥5) cannot be a single NURBS surface; \
                 split into 4-sided sub-patches or use a mesh-domain fill"
            ),
        )),
    }
}

/// Fill a triangular patch as a degenerate Coons patch. `d0` and `d1`
/// are the two rail edges (sharing the apex vertex); `c0` is the base
/// edge. The `c1` side is collapsed to the apex.
fn fill_triangular(
    d0: &NurbsCurve,
    d1: &NurbsCurve,
    c0: &NurbsCurve,
) -> Result<NurbsSurface, OcctSurfaceError> {
    // The apex is where d0 and d1 both end (their v = max ends).
    let apex_d0 = d0.evaluate(d0.parameter_range().1);
    let apex_d1 = d1.evaluate(d1.parameter_range().1);
    let apex = 0.5 * (apex_d0 + apex_d1);
    // The collapsed c1 side is a zero-length curve at the apex.
    let c1 = NurbsCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![apex, apex],
        vec![1.0, 1.0],
    )
    .map_err(|e| OcctSurfaceError::bad_input("boundary", format!("degenerate edge: {e}")))?;
    // Coons order: [c0, c1, d0, d1].
    fill([c0.clone(), c1, d0.clone(), d1.clone()]).map_err(|e| {
        OcctSurfaceError::bad_input("boundary", format!("triangular coons fill: {e}"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// A straight line as a cubic Bezier — gives an exact degree-3
    /// representation that `coons::fill` reproduces precisely.
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
    fn filling_rejects_too_few_boundaries() {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(1.0, 0.0, 0.0);
        let curves = vec![line(a, b), line(b, a)];
        let err = offset_api_filling(&curves, &[]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn filling_rejects_five_sided() {
        let p = Vector3::new(0.0, 0.0, 0.0);
        let pentagon: Vec<NurbsCurve> = (0..5).map(|_| line(p, p + Vector3::x())).collect();
        let err = offset_api_filling(&pentagon, &[]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
        assert!(err.to_string().contains("N≥5") || err.to_string().contains("5-sided"));
    }

    #[test]
    fn filling_four_sided_unit_square_is_a_flat_patch() {
        // Coons order [c0, c1, d0, d1] for the unit square in z=0.
        let c0 = line(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)); // v=0
        let c1 = line(Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 0.0)); // v=1
        let d0 = line(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)); // u=0
        let d1 = line(Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)); // u=1
        let surf = offset_api_filling(&[c0, c1, d0, d1], &[]).unwrap();
        // Every point on the filled patch lies in the z=0 plane.
        for &(u, v) in &[(0.0, 0.0), (0.5, 0.5), (1.0, 0.3), (0.2, 1.0)] {
            let p = surf.evaluate(u, v);
            assert!(p.z.abs() < 1e-9, "patch off the z=0 plane at ({u},{v}): {}", p.z);
        }
        // The corners interpolate the boundary corners exactly.
        let p00 = surf.evaluate(0.0, 0.0);
        assert!((p00 - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn filling_triangular_patch_interpolates_the_three_edges() {
        // Triangle: rails d0, d1 share the apex (0.5, 1); base c0.
        let d0 = line(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.5, 1.0, 0.0));
        let d1 = line(Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.5, 1.0, 0.0));
        let c0 = line(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        let surf = offset_api_filling(&[d0, d1, c0], &[]).unwrap();
        // The patch lies in z=0.
        let mid = surf.evaluate(0.5, 0.5);
        assert!(mid.z.abs() < 1e-9, "triangular patch off z=0: {}", mid.z);
        // The collapsed (v=1) side is the apex.
        let apex = surf.evaluate(0.5, 1.0);
        assert!(
            (apex - Vector3::new(0.5, 1.0, 0.0)).norm() < 1e-6,
            "v=1 side should collapse to the apex, got {apex:?}"
        );
    }
}
