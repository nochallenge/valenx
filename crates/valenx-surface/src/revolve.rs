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

/// Revolve `profile` about the Z-axis through `angle` radians (in `(0, 2π]`)
/// into a NURBS surface of revolution.
///
/// The around-axis (`u`) direction is an **exact** rational arc built from
/// `⌈angle / 90°⌉` quadratic segments of equal span (Piegl & Tiller, *The NURBS
/// Book*, alg. A8.3): `2·narcs + 1` control points, degree 2, with the arc
/// shoulders at radius `r / cos(Δθ/2)` and weight `cos(Δθ/2)`. The along-profile
/// (`v`) direction reuses the profile's degree, knots and weights, so a rational
/// profile (an arc, a conic) is revolved exactly. Only each profile point's
/// distance from the Z-axis is swept, so the profile may lie in any half-plane
/// containing the axis.
///
/// # Errors
///
/// Returns [`SurfaceError::BadBoundary`] if `angle` is not in `(0, 2π]`, or if
/// every profile control point lies on the axis (radius ≤ `AXIS_EPS`) so there is
/// nothing to revolve. Propagates any [`SurfaceError`] from [`NurbsSurface::new`].
pub fn revolve_z(profile: &NurbsCurve, angle: f64) -> Result<NurbsSurface, SurfaceError> {
    if !(angle > 0.0 && angle <= std::f64::consts::TAU + AXIS_EPS) {
        return Err(SurfaceError::BadBoundary(format!(
            "revolution angle {angle} must be in (0, 2π]"
        )));
    }
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

    // Split the sweep into narcs ≤ 90° quadratic arcs.
    let narcs = (angle / std::f64::consts::FRAC_PI_2).ceil().max(1.0) as usize;
    let d_theta = angle / narcs as f64;
    let shoulder_w = (d_theta / 2.0).cos(); // arc shoulder weight
    let n_u = 2 * narcs + 1; // number of u control points

    // Clamped quadratic knot vector with each interior arc boundary doubled.
    let mut u_knots = vec![0.0, 0.0, 0.0];
    for k in 1..narcs {
        let v = k as f64 / narcs as f64;
        u_knots.push(v);
        u_knots.push(v);
    }
    u_knots.push(1.0);
    u_knots.push(1.0);
    u_knots.push(1.0);

    let nv = profile.control_points.len();
    let mut control_points: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(n_u);
    let mut weights: Vec<Vec<f64>> = Vec::with_capacity(n_u);

    for i in 0..n_u {
        // Even i: a corner on the circle at angle (i/2)·Δθ, radius r, weight 1.
        // Odd  i: an arc shoulder at the half-step angle, radius r/cos(Δθ/2),
        //         weight cos(Δθ/2).
        let shoulder = i % 2 == 1;
        let theta = i as f64 * 0.5 * d_theta;
        let (sin_t, cos_t) = theta.sin_cos();
        let radius_mult = if shoulder { 1.0 / shoulder_w } else { 1.0 };
        let circle_w = if shoulder { shoulder_w } else { 1.0 };

        let mut row_cp = Vec::with_capacity(nv);
        let mut row_w = Vec::with_capacity(nv);
        for (p, &wp) in profile.control_points.iter().zip(&profile.weights) {
            let radius = (p.x * p.x + p.y * p.y).sqrt();
            row_cp.push(Vector3::new(
                radius * radius_mult * cos_t,
                radius * radius_mult * sin_t,
                p.z,
            ));
            row_w.push(circle_w * wp);
        }
        control_points.push(row_cp);
        weights.push(row_w);
    }

    NurbsSurface::new(
        2,
        profile.degree,
        u_knots,
        profile.knots.clone(),
        control_points,
        weights,
    )
}

/// Revolve `profile` a full 360° about the Z-axis — a convenience wrapper for
/// [`revolve_z`] with `angle = 2π` (the around-axis circle is the exact
/// nine-control-point, four-arc NURBS circle).
///
/// # Errors
///
/// As [`revolve_z`].
pub fn revolve_z_full(profile: &NurbsCurve) -> Result<NurbsSurface, SurfaceError> {
    revolve_z(profile, std::f64::consts::TAU)
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

    #[test]
    fn partial_cylinder_area_scales_with_angle() {
        // Revolving a vertical line through θ gives a partial cylinder whose
        // lateral area is the arc length r·θ times the height: θ·r·h.
        let (r, h) = (2.0, 5.0);
        for &theta in &[std::f64::consts::FRAC_PI_2, PI, std::f64::consts::PI * 1.5] {
            let surface = revolve_z(&line_profile(r, 0.0, r, h), theta).unwrap();
            let area = surface_area(&surface);
            let analytic = theta * r * h;
            assert!(
                (area - analytic).abs() / analytic < 0.01,
                "θ={theta:.4}: area {area:.4} vs analytic {analytic:.4}"
            );
        }
    }

    #[test]
    fn narcs_sets_the_control_point_count() {
        // 2·⌈θ / 90°⌉ + 1 control points around the axis.
        let line = line_profile(1.0, 0.0, 1.0, 1.0);
        for (theta, expected) in [(1.0, 3), (2.0, 5), (4.0, 7), (std::f64::consts::TAU, 9)] {
            let s = revolve_z(&line, theta).unwrap();
            assert_eq!(s.control_points.len(), expected, "θ={theta:.4}");
        }
    }

    #[test]
    fn full_wrapper_matches_two_pi_revolve() {
        let line = line_profile(2.0, 0.0, 2.0, 5.0);
        let a_full = surface_area(&revolve_z_full(&line).unwrap());
        let a_2pi = surface_area(&revolve_z(&line, std::f64::consts::TAU).unwrap());
        assert!(
            (a_full - a_2pi).abs() < 1e-9,
            "full wrapper {a_full} vs revolve_z(2π) {a_2pi}"
        );
    }

    #[test]
    fn out_of_range_angle_is_an_error() {
        let line = line_profile(1.0, 0.0, 1.0, 1.0);
        assert!(revolve_z(&line, 0.0).is_err());
        assert!(revolve_z(&line, -1.0).is_err());
        assert!(revolve_z(&line, std::f64::consts::TAU + 0.5).is_err());
    }
}
