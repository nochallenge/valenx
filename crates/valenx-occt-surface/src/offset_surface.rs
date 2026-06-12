//! Phase 72 — `BRepFill_OffsetSurface`: offset a surface by a signed
//! distance along its normal.
//!
//! ## What OCCT does
//!
//! `BRepFill_OffsetSurface` (and the closely related
//! `Geom_OffsetSurface`) computes a parametric surface whose value at
//! `(u, v)` is `S(u, v) + d * N(u, v)`, where `N` is the unit normal
//! to `S` at that point. The result is a `Geom_OffsetSurface` BRep
//! face that callers can sew into shells, intersect, or tessellate
//! like any other face. Caveats:
//!
//! - Self-intersection — offsetting a curved surface by a distance
//!   greater than its minimum radius of curvature folds the offset
//!   surface back on itself. OCCT detects this and the calling code
//!   has to either trim the result or reject the operation.
//! - Tolerance — the offset preserves G2 continuity of the base but
//!   `Geom_OffsetSurface` evaluation cost is higher; callers usually
//!   refit with `BRepApprox_Approx` for downstream IGES/STEP export.
//!
//! ## v1 status — real implementation (control-point displacement)
//!
//! truck does not expose an on-the-fly surface-normal hook, so the
//! exact `S + d·N` offset surface is not directly representable. This
//! module ships the **control-point-displacement offset**, the
//! standard approximate offset used by every NURBS kernel that wants
//! to keep the result a NURBS surface:
//!
//! 1. Each control point `P[i][j]` is mapped to its parametric
//!    location via the *Greville abscissae* — the average of the
//!    `degree` knots spanning that control point. The Greville point
//!    is where `P[i][j]` exerts the most influence on the surface.
//! 2. The unit surface normal `N` is evaluated there by central
//!    finite differences of the existing tensor-product evaluator
//!    (`∂S/∂u × ∂S/∂v`, normalised).
//! 3. The control point is displaced: `P'[i][j] = P[i][j] + d·N`.
//!    Weights and the knot vectors are carried over unchanged.
//!
//! ### Accuracy
//!
//! - **Exact** for planar surfaces (the normal is constant, so
//!   displacing every CP by `d·N` translates the whole plane by
//!   exactly `d`).
//! - **Good** for gently curved surfaces where `d` is small relative
//!   to the local radius of curvature.
//! - **Approximate** elsewhere: the offset of a NURBS surface is not
//!   in general a NURBS surface of the same degree, so the result is
//!   the best same-degree NURBS fit. The deviation grows with `d/R`
//!   (offset distance over radius of curvature). For a rigorous
//!   offset, refit with more control points or use a higher-degree
//!   surface — that refit is the documented follow-up.

use nalgebra::Vector3;
use valenx_surface::NurbsSurface;

use crate::error::OcctSurfaceError;

/// Offset a NURBS surface by `distance` along its normal. Positive
/// distance offsets in the +N direction (`∂S/∂u × ∂S/∂v`); negative
/// offsets in -N.
///
/// See the module docs for the algorithm and the accuracy envelope.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] when `distance` is non-finite.
pub fn offset_surface(
    surface: &NurbsSurface,
    distance: f64,
) -> Result<NurbsSurface, OcctSurfaceError> {
    if !distance.is_finite() {
        return Err(OcctSurfaceError::bad_input(
            "distance",
            format!("must be finite, got {distance}"),
        ));
    }
    if distance == 0.0 {
        return Ok(surface.clone());
    }

    let nu = surface.nu();
    let nv = surface.nv();

    // Greville abscissae: the parameter associated with each control
    // point is the average of the `degree` knots after its index.
    let greville_u: Vec<f64> = (0..nu)
        .map(|i| greville(&surface.u_knots, i, surface.u_degree))
        .collect();
    let greville_v: Vec<f64> = (0..nv)
        .map(|j| greville(&surface.v_knots, j, surface.v_degree))
        .collect();

    let (u_lo, u_hi) = surface.u_range();
    let (v_lo, v_hi) = surface.v_range();

    // Displace every control point along the surface normal at its
    // Greville point.
    let mut new_cps: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(nu);
    for (i, gu) in greville_u.iter().enumerate() {
        let mut row: Vec<Vector3<f64>> = Vec::with_capacity(nv);
        for (j, gv) in greville_v.iter().enumerate() {
            let u = gu.clamp(u_lo, u_hi);
            let v = gv.clamp(v_lo, v_hi);
            let n = surface_normal(surface, u, v, u_lo, u_hi, v_lo, v_hi);
            row.push(surface.control_points[i][j] + n * distance);
        }
        new_cps.push(row);
    }

    NurbsSurface::new(
        surface.u_degree,
        surface.v_degree,
        surface.u_knots.clone(),
        surface.v_knots.clone(),
        new_cps,
        surface.weights.clone(),
    )
    .map_err(|e| OcctSurfaceError::bad_input("surface", format!("offset surface invalid: {e}")))
}

/// Greville abscissa for control point `i` of a B-spline with the
/// given knot vector and degree: `(k_{i+1} + … + k_{i+degree}) /
/// degree`.
fn greville(knots: &[f64], i: usize, degree: usize) -> f64 {
    if degree == 0 {
        return knots.get(i).copied().unwrap_or(0.0);
    }
    let mut sum = 0.0;
    for d in 1..=degree {
        sum += knots.get(i + d).copied().unwrap_or(0.0);
    }
    sum / degree as f64
}

/// Unit surface normal at `(u, v)` via central finite differences of
/// the tensor-product evaluator. The step is shrunk near the
/// parameter-domain boundary so the difference stays inside the valid
/// range; one-sided differences are used at the very edge.
#[allow(clippy::too_many_arguments)]
fn surface_normal(
    s: &NurbsSurface,
    u: f64,
    v: f64,
    u_lo: f64,
    u_hi: f64,
    v_lo: f64,
    v_hi: f64,
) -> Vector3<f64> {
    let du = ((u_hi - u_lo).abs() * 1e-4).max(1e-9);
    let dv = ((v_hi - v_lo).abs() * 1e-4).max(1e-9);

    // ∂S/∂u — central where possible, else one-sided.
    let su = if u - du >= u_lo && u + du <= u_hi {
        (s.evaluate(u + du, v) - s.evaluate(u - du, v)) / (2.0 * du)
    } else if u + du <= u_hi {
        (s.evaluate(u + du, v) - s.evaluate(u, v)) / du
    } else {
        (s.evaluate(u, v) - s.evaluate(u - du, v)) / du
    };
    // ∂S/∂v.
    let sv = if v - dv >= v_lo && v + dv <= v_hi {
        (s.evaluate(u, v + dv) - s.evaluate(u, v - dv)) / (2.0 * dv)
    } else if v + dv <= v_hi {
        (s.evaluate(u, v + dv) - s.evaluate(u, v)) / dv
    } else {
        (s.evaluate(u, v) - s.evaluate(u, v - dv)) / dv
    };

    let n = su.cross(&sv);
    let len = n.norm();
    if len < 1e-12 {
        // Degenerate parametrisation (e.g. a pole). Fall back to +Z
        // so the offset still produces a finite surface.
        Vector3::z()
    } else {
        n / len
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// A flat 2x2 NURBS surface in the z=0 plane.
    ///
    /// `control_points[i][j]` is indexed `[i (u)][j (v)]`, so the u
    /// direction must run along +x and the v direction along +y for the
    /// parametric normal `∂S/∂u × ∂S/∂v = x̂ × ŷ` to point at +z — the
    /// same orientation as [`curved_3x3`]. (The earlier layout ran u
    /// along +y and v along +x, which gave a -z parametric normal and
    /// made an offset along the surface normal read as a downward
    /// translation.)
    fn planar_2x2() -> NurbsSurface {
        let cps = vec![
            vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)],
            vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            cps,
            weights,
        )
        .expect("planar 2x2 builds")
    }

    /// A cylindrical-ish curved 3x3 quadratic surface bulging in +z.
    fn curved_3x3() -> NurbsSurface {
        let mut cps = vec![vec![Vector3::zeros(); 3]; 3];
        for (i, row) in cps.iter_mut().enumerate() {
            for (j, p) in row.iter_mut().enumerate() {
                let x = i as f64;
                let y = j as f64;
                // Middle row bulges up.
                let z = if i == 1 { 1.0 } else { 0.0 };
                *p = Vector3::new(x, y, z);
            }
        }
        let weights = vec![vec![1.0; 3]; 3];
        NurbsSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            cps,
            weights,
        )
        .expect("curved 3x3 builds")
    }

    #[test]
    fn offset_surface_rejects_nan_distance() {
        let surf = planar_2x2();
        let err = offset_surface(&surf, f64::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn zero_offset_is_identity() {
        let surf = planar_2x2();
        let out = offset_surface(&surf, 0.0).unwrap();
        for i in 0..surf.nu() {
            for j in 0..surf.nv() {
                assert!((out.control_points[i][j] - surf.control_points[i][j]).norm() < 1e-12);
            }
        }
    }

    #[test]
    fn planar_offset_is_exact_translation() {
        // A flat surface in z=0, offset by +2, must land every point
        // on z=2 exactly (the planar case is exact for CP displacement).
        let surf = planar_2x2();
        let out = offset_surface(&surf, 2.0).unwrap();
        for &(u, v) in &[(0.0, 0.0), (0.5, 0.5), (1.0, 0.3), (0.2, 1.0)] {
            let p = out.evaluate(u, v);
            assert!((p.z - 2.0).abs() < 1e-9, "z={} at ({u},{v})", p.z);
        }
    }

    #[test]
    fn negative_offset_flips_direction() {
        let surf = planar_2x2();
        let out = offset_surface(&surf, -1.5).unwrap();
        let p = out.evaluate(0.5, 0.5);
        assert!((p.z + 1.5).abs() < 1e-9, "z={}", p.z);
    }

    #[test]
    fn curved_offset_moves_surface_outward() {
        // A bulging surface offset by +0.5 must move the mid-point
        // strictly further from the original mid-point — the offset
        // is a real displacement along the local normal.
        let surf = curved_3x3();
        let mid_before = surf.evaluate(0.5, 0.5);
        let out = offset_surface(&surf, 0.5).unwrap();
        let mid_after = out.evaluate(0.5, 0.5);
        let moved = (mid_after - mid_before).norm();
        assert!(
            moved > 0.3,
            "curved offset should displace the surface, moved={moved}"
        );
        // The CP grid keeps its shape.
        assert_eq!(out.nu(), 3);
        assert_eq!(out.nv(), 3);
    }
}
