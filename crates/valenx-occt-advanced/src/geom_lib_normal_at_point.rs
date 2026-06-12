//! Phase 154 — `GeomLib::Normal` for surfaces — return the unit
//! surface normal at parametric `(u, v)`.
//!
//! ## What OCCT does
//!
//! `Geom_Surface::D1(u, v, P, dU, dV)` returns the surface point plus
//! its `u`-partial and `v`-partial derivatives. The outward unit
//! normal is then `n = (dU × dV) / |dU × dV|`. Sign convention: OCCT
//! uses the orientation of the parametric grid (right-hand rule on
//! `(dU, dV)`) — flipping `u` ↔ `v` flips the normal.
//!
//! ## v1 status
//!
//! **Honest v1.** Uses `NurbsSurface::evaluate` to sample three
//! offset points (centre, +ε in u, +ε in v) and computes the normal
//! via finite differences. Analytical derivatives (via control-point
//! shifts) are Phase 154.5 — currently the central-difference variant
//! is fast enough for tessellation-grade accuracy.

use nalgebra::Vector3;
use valenx_surface::NurbsSurface;

use crate::error::OcctAdvancedError;

/// Step size for finite-difference partial derivatives. 1e-5 in
/// parameter space gives 6-7 digits of normal accuracy for cubic
/// NURBS — same precision contract as
/// `NurbsCurve::derivative`.
pub const DERIVATIVE_STEP: f64 = 1e-5;

/// Compute the unit normal at `(u, v)` on `surface`.
///
/// Returns the unit normal vector; surfaces a [`OcctAdvancedError::Defect`]
/// when the cross product magnitude falls below `1e-12` (degenerate
/// parameterisation — surface has zero tangent plane at this point).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] when `(u, v)` falls outside the
///   surface's valid parameter range.
/// - [`OcctAdvancedError::Defect`] for a degenerate tangent plane.
pub fn geom_lib_normal_at_point(
    surface: &NurbsSurface,
    u: f64,
    v: f64,
) -> Result<Vector3<f64>, OcctAdvancedError> {
    let (u_min, u_max) = surface.u_range();
    let (v_min, v_max) = surface.v_range();
    if !u.is_finite() || !v.is_finite() {
        return Err(OcctAdvancedError::bad_input(
            "u/v",
            "parameters must be finite",
        ));
    }
    if u < u_min || u > u_max {
        return Err(OcctAdvancedError::bad_input(
            "u",
            format!("{u} outside [{u_min}, {u_max}]"),
        ));
    }
    if v < v_min || v > v_max {
        return Err(OcctAdvancedError::bad_input(
            "v",
            format!("{v} outside [{v_min}, {v_max}]"),
        ));
    }

    // Central differences in u, v — clamp near the boundary so we
    // don't sample outside the valid range.
    let h_u = DERIVATIVE_STEP * (u_max - u_min);
    let h_v = DERIVATIVE_STEP * (v_max - v_min);
    let u_lo = (u - h_u).max(u_min);
    let u_hi = (u + h_u).min(u_max);
    let v_lo = (v - h_v).max(v_min);
    let v_hi = (v + h_v).min(v_max);

    let dudu = u_hi - u_lo;
    let dvdv = v_hi - v_lo;
    if dudu.abs() < 1e-30 || dvdv.abs() < 1e-30 {
        return Err(OcctAdvancedError::defect(
            format!("(u={u},v={v})"),
            "parameter step degenerated at boundary",
        ));
    }

    let p_u_lo = surface.evaluate(u_lo, v);
    let p_u_hi = surface.evaluate(u_hi, v);
    let p_v_lo = surface.evaluate(u, v_lo);
    let p_v_hi = surface.evaluate(u, v_hi);

    let du = (p_u_hi - p_u_lo) / dudu;
    let dv = (p_v_hi - p_v_lo) / dvdv;
    let n = du.cross(&dv);
    let n_norm = n.norm();
    if n_norm < 1e-12 {
        return Err(OcctAdvancedError::defect(
            format!("(u={u},v={v})"),
            format!("|dU × dV| = {n_norm:.3e} (degenerate tangent plane)"),
        ));
    }
    Ok(n / n_norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_surface::{coons, NurbsCurve};

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

    fn unit_square_xy() -> NurbsSurface {
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p01 = Vector3::new(0.0, 1.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);
        coons::fill([
            line(p00, p01),
            line(p10, p11),
            line(p00, p10),
            line(p01, p11),
        ])
        .unwrap()
    }

    #[test]
    fn rejects_out_of_range_u() {
        let s = unit_square_xy();
        let err = geom_lib_normal_at_point(&s, -0.1, 0.5).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_non_finite_v() {
        let s = unit_square_xy();
        let err = geom_lib_normal_at_point(&s, 0.5, f64::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn flat_xy_plane_normal_is_z() {
        let s = unit_square_xy();
        let n = geom_lib_normal_at_point(&s, 0.5, 0.5).unwrap();
        // Coons fill of a unit-XY square is a flat plane in z=0; its
        // normal should be ±z. Coons orientation gives +z.
        assert!(n.x.abs() < 1e-6, "nx should be ~0, got {}", n.x);
        assert!(n.y.abs() < 1e-6, "ny should be ~0, got {}", n.y);
        assert!(
            (n.z.abs() - 1.0).abs() < 1e-6,
            "|nz| should be ~1, got {}",
            n.z
        );
    }
}
