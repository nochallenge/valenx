//! Surfaces of revolution.
//!
//! [`revolve_z_full`] sweeps a profile curve a full 360° about the Z-axis,
//! producing a NURBS surface whose `u` (around-axis) direction is an **exact**
//! rational-quadratic circle and whose `v` (along-profile) direction reuses the
//! profile's own degree, knots and weights. This is the standard Piegl & Tiller
//! revolved-surface construction (*The NURBS Book*, alg. A8.3) specialised to a
//! full turn: four 90° arcs give nine control points in `u`, degree 2, with the
//! alternating circle weights `1, √2/2, 1, …`.
//!
//! Because a full turn closes the circle, only each profile control point's
//! **distance from the Z-axis** matters — the profile may lie in any half-plane
//! containing the axis (the conventional choice is the `xz` half-plane with
//! `x ≥ 0`). A point on the axis becomes a pole (a degenerate control-point row),
//! which is exactly what a cone apex or a sphere pole needs.
//!
//! The construction is geometrically exact; the only approximation when you
//! *measure* the result (e.g. with [`crate::surface_area`]) is that area's own
//! numerical quadrature. The module tests confirm the canonical analytic areas —
//! cylinder `2πrh`, cone `πr·√(r²+h²)`, sphere `4πr²` — to within that tolerance.

use nalgebra::Vector3;

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;

/// Profile radii at or below this (in model units) are treated as lying on the
/// axis of revolution.
const AXIS_EPS: f64 = 1e-9;

/// Revolve `profile` a full 360° about the Z-axis into a NURBS surface of
/// revolution.
///
/// The resulting surface has `u_degree = 2` (the around-axis circle) and
/// `v_degree = profile.degree`; its `v` knots and per-row weight factors come
/// from the profile, so a rational profile (an arc, a conic) is revolved exactly.
///
/// # Errors
///
/// Returns [`SurfaceError::BadBoundary`] if every profile control point lies on
/// the axis (radius ≤ `AXIS_EPS`), so there is nothing to revolve. Propagates
/// any [`SurfaceError`] from [`NurbsSurface::new`] (which re-validates the
/// generated control net).
pub fn revolve_z_full(profile: &NurbsCurve) -> Result<NurbsSurface, SurfaceError> {
    // A profile that lies entirely on the axis revolves to (at most) a line.
    let has_radius = profile
        .control_points
        .iter()
        .any(|p| (p.x * p.x + p.y * p.y).sqrt() > AXIS_EPS);
    if !has_radius {
        return Err(SurfaceError::BadBoundary(
            "profile lies on the axis of revolution; nothing to revolve".to_string(),
        ));
    }

    // Full circle in u: 4 quadratic arcs, 9 control points, degree 2.
    let u_degree = 2;
    let u_knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];

    let nv = profile.control_points.len();
    let mut control_points: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(9);
    let mut weights: Vec<Vec<f64>> = Vec::with_capacity(9);

    for i in 0..9 {
        // Control points sit every 45°; the odd ones are the arc "shoulders" at
        // the corners of the circumscribing square (radius ×√2, weight √2/2).
        let angle = f64::from(i as u32) * std::f64::consts::FRAC_PI_4;
        let (sin_a, cos_a) = angle.sin_cos();
        let shoulder = i % 2 == 1;
        let radius_mult = if shoulder {
            std::f64::consts::SQRT_2
        } else {
            1.0
        };
        let circle_w = if shoulder {
            std::f64::consts::FRAC_1_SQRT_2
        } else {
            1.0
        };

        let mut row_cp = Vec::with_capacity(nv);
        let mut row_w = Vec::with_capacity(nv);
        for (p, &wp) in profile.control_points.iter().zip(&profile.weights) {
            let radius = (p.x * p.x + p.y * p.y).sqrt();
            row_cp.push(Vector3::new(
                radius * radius_mult * cos_a,
                radius * radius_mult * sin_a,
                p.z,
            ));
            row_w.push(circle_w * wp);
        }
        control_points.push(row_cp);
        weights.push(row_w);
    }

    NurbsSurface::new(
        u_degree,
        profile.degree,
        u_knots,
        profile.knots.clone(),
        control_points,
        weights,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::area::surface_area;
    use std::f64::consts::{PI, SQRT_2};

    /// Straight degree-1 profile from `(x0,0,z0)` to `(x1,0,z1)`.
    fn line_profile(x0: f64, z0: f64, x1: f64, z1: f64) -> NurbsCurve {
        NurbsCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Vector3::new(x0, 0.0, z0), Vector3::new(x1, 0.0, z1)],
            vec![1.0, 1.0],
        )
        .unwrap()
    }

    #[test]
    fn cylinder_area_matches_analytic() {
        let (r, h) = (2.0, 5.0);
        let surface = revolve_z_full(&line_profile(r, 0.0, r, h)).unwrap();
        // Surface is 9 (u) × 2 (v) control points.
        assert_eq!(surface.control_points.len(), 9);
        assert_eq!(surface.control_points[0].len(), 2);
        let area = surface_area(&surface);
        let analytic = 2.0 * PI * r * h; // lateral area of a cylinder
        assert!(
            (area - analytic).abs() / analytic < 0.01,
            "cylinder area {area:.4} vs analytic {analytic:.4}"
        );
    }

    #[test]
    fn cone_area_matches_analytic() {
        // Apex on the axis at height h, base circle of radius r at z = 0.
        let (r, h) = (3.0, 4.0);
        let surface = revolve_z_full(&line_profile(0.0, h, r, 0.0)).unwrap();
        let area = surface_area(&surface);
        let slant = (r * r + h * h).sqrt();
        let analytic = PI * r * slant; // lateral area of a cone (= π·3·5)
        assert!(
            (area - analytic).abs() / analytic < 0.01,
            "cone area {area:.4} vs analytic {analytic:.4}"
        );
    }

    #[test]
    fn sphere_area_matches_analytic() {
        // East semicircle of radius R in the xz half-plane, as two 90° rational
        // arcs (5 control points): north pole -> equator -> south pole.
        let r = 2.0;
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let profile = NurbsCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
            vec![
                Vector3::new(0.0, 0.0, r),  // north pole (on axis)
                Vector3::new(r, 0.0, r),    // shoulder
                Vector3::new(r, 0.0, 0.0),  // equator
                Vector3::new(r, 0.0, -r),   // shoulder
                Vector3::new(0.0, 0.0, -r), // south pole (on axis)
            ],
            vec![1.0, s, 1.0, s, 1.0],
        )
        .unwrap();
        let surface = revolve_z_full(&profile).unwrap();
        let area = surface_area(&surface);
        let analytic = 4.0 * PI * r * r; // area of a sphere
        assert!(
            (area - analytic).abs() / analytic < 0.02,
            "sphere area {area:.4} vs analytic {analytic:.4}"
        );
    }

    #[test]
    fn weights_carry_the_circle_pattern() {
        let surface = revolve_z_full(&line_profile(1.0, 0.0, 1.0, 1.0)).unwrap();
        // Corner rows weight 1, shoulder rows √2/2 (profile weights are 1 here).
        assert!((surface.weights[0][0] - 1.0).abs() < 1e-12);
        assert!((surface.weights[1][0] - (SQRT_2 / 2.0)).abs() < 1e-12);
        assert!((surface.weights[2][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn profile_on_axis_is_an_error() {
        // A profile with zero radius everywhere has nothing to revolve.
        let err = revolve_z_full(&line_profile(0.0, 0.0, 0.0, 3.0)).unwrap_err();
        assert!(matches!(err, SurfaceError::BadBoundary(_)));
    }
}
