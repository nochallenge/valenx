//! Surface area and area-weighted centroid of a `NurbsSurface`
//! (Phase 19, geometric-properties extension).
//!
//! The surface area of a parametric surface `S(u, v)` is the double
//! integral of the magnitude of the cross product of its first partial
//! derivatives,
//!
//! `A = ∫∫ |S_u × S_v| du dv`,
//!
//! taken over the valid parameter rectangle `u_range × v_range`. The
//! integrand `|S_u × S_v|` is the local area-scaling (Jacobian) of the
//! parameterisation: it is exactly the area of the differential
//! parallelogram spanned by the two tangent vectors. The area-weighted
//! centroid uses the same differential area element as the weight,
//!
//! `c = (1/A) ∫∫ S(u, v) |S_u × S_v| du dv`.
//!
//! Both integrals are approximated here by a fine **composite midpoint
//! rule** on a uniform `samples × samples` grid of cells covering the
//! parameter rectangle. The two partial derivatives at each cell
//! centre are obtained by **central finite differences** of
//! `NurbsSurface::evaluate`, with the difference step clamped so the
//! stencil never leaves the valid knot range. This keeps the routine
//! self-contained — it relies only on the public `evaluate`,
//! `u_range`, and `v_range` API — at the cost of the small truncation
//! error inherent to a finite-difference tangent.
//!
//! ## Convergence
//!
//! The midpoint rule converges as `O(1/samples²)` for a smooth
//! integrand. For a flat (planar) surface the parameter Jacobian is
//! constant and the rule is exact regardless of `samples`. For a
//! curved surface the integrand varies and a moderate grid (the
//! default `64`) brings the relative error to well under one percent
//! on the smooth test cases in this module. Surfaces with sharp
//! curvature, rational weight spikes, or interior knot lines may need
//! a larger `samples` value for the same accuracy.
//!
//! ## Scope caveat
//!
//! This is a research / preliminary-design grade quadrature, intended
//! for mass-property estimates and sanity checks. It is **not** an
//! adaptive, error-bounded integrator and makes **no** claim of parity
//! with the validated mass-property engines in CATIA, Ansys, or Adams:
//! it does not refine adaptively, does not specially handle degenerate
//! (zero-area) parameter regions or seams, and reports no error bound.
//! Treat the returned values as approximations whose accuracy is set
//! by the chosen sample density.

use nalgebra::Vector3;

use crate::nurbs_surface::NurbsSurface;

/// Default number of integration cells per parametric direction.
///
/// The quadrature evaluates the surface on a `DEFAULT_SAMPLES ×
/// DEFAULT_SAMPLES` grid of midpoints. `64` is a balance between
/// accuracy and cost that resolves the smooth test cases to better
/// than one percent.
pub const DEFAULT_SAMPLES: usize = 64;

/// Result of an area / centroid computation over a `NurbsSurface`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaResult {
    /// Approximate surface area, `∫∫ |S_u × S_v| du dv`, in the same
    /// length-squared units as the surface control points.
    pub area: f64,
    /// Approximate area-weighted centroid (the geometric centre of
    /// mass of an infinitely thin shell of uniform surface density).
    /// If the computed `area` is numerically zero the centroid falls
    /// back to the surface point at the centre of the parameter
    /// rectangle.
    pub centroid: Vector3<f64>,
}

/// Estimate the surface area of `surface` by composite-midpoint
/// quadrature of `|S_u × S_v|` over `u_range × v_range`, using
/// `DEFAULT_SAMPLES` cells per direction.
///
/// This is a convenience wrapper over `area_and_centroid` that
/// discards the centroid. See the module documentation for the method
/// and its accuracy caveats.
pub fn surface_area(surface: &NurbsSurface) -> f64 {
    area_and_centroid_with_samples(surface, DEFAULT_SAMPLES).area
}

/// Estimate both the surface area and the area-weighted centroid of
/// `surface` using `DEFAULT_SAMPLES` cells per parametric direction.
///
/// See the module documentation for the quadrature method and its
/// accuracy caveats.
pub fn area_and_centroid(surface: &NurbsSurface) -> AreaResult {
    area_and_centroid_with_samples(surface, DEFAULT_SAMPLES)
}

/// Like `area_and_centroid`, but with a caller-chosen number of
/// integration cells per parametric direction.
///
/// `samples` is the count of midpoint cells along **each** of the u
/// and v directions, so the integrand is evaluated at `samples ×
/// samples` interior points. A value of `0` is treated as `1`. Larger
/// values reduce the quadrature error (roughly as `1/samples²` for a
/// smooth surface) at proportionally higher cost.
pub fn area_and_centroid_with_samples(surface: &NurbsSurface, samples: usize) -> AreaResult {
    let n = samples.max(1);
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();

    let du = (u1 - u0) / n as f64;
    let dv = (v1 - v0) / n as f64;
    let cell = du * dv;

    // Central-difference steps: a small fraction of the cell size,
    // bounded below so a tiny (or zero-length) parameter span cannot
    // collapse the stencil to a zero step.
    let hu = finite_diff_step(du, u1 - u0);
    let hv = finite_diff_step(dv, v1 - v0);

    let mut area = 0.0_f64;
    let mut weighted = Vector3::zeros();

    for iu in 0..n {
        // Midpoint of the iu-th cell in u.
        let u = u0 + (iu as f64 + 0.5) * du;
        for iv in 0..n {
            let v = v0 + (iv as f64 + 0.5) * dv;

            let s_u = partial_u(surface, u, v, hu, u0, u1);
            let s_v = partial_v(surface, u, v, hv, v0, v1);
            // |S_u × S_v| is the local area scaling; times the cell
            // area du·dv gives this cell's contribution.
            let jac = s_u.cross(&s_v).norm();
            let da = jac * cell;

            area += da;
            weighted += surface.evaluate(u, v) * da;
        }
    }

    let centroid = if area > f64::EPSILON {
        weighted / area
    } else {
        // Degenerate (near-zero-area) surface: fall back to the
        // geometric centre of the parameter rectangle so the centroid
        // is still a sensible point on the surface.
        surface.evaluate(0.5 * (u0 + u1), 0.5 * (v0 + v1))
    };

    AreaResult { area, centroid }
}

/// Choose a central-difference step: a small fraction of the local
/// cell size, but never larger than a fraction of the whole span and
/// never collapsing to zero for a degenerate span.
fn finite_diff_step(cell_size: f64, span: f64) -> f64 {
    let by_cell = cell_size * 0.5;
    let by_span = span.abs() * 1.0e-3;
    let candidate = by_cell.min(by_span.max(0.0));
    if candidate > 0.0 {
        candidate
    } else {
        // Span is zero / degenerate — use an absolute fallback so the
        // stencil still has a finite width.
        1.0e-6
    }
}

/// Central finite difference of `S` in u at `(u, v)`, clamped so the
/// stencil stays inside `[u0, u1]`.
fn partial_u(surface: &NurbsSurface, u: f64, v: f64, h: f64, u0: f64, u1: f64) -> Vector3<f64> {
    let lo = (u - h).max(u0);
    let hi = (u + h).min(u1);
    let denom = hi - lo;
    if denom <= 0.0 {
        return Vector3::zeros();
    }
    (surface.evaluate(hi, v) - surface.evaluate(lo, v)) / denom
}

/// Central finite difference of `S` in v at `(u, v)`, clamped so the
/// stencil stays inside `[v0, v1]`.
fn partial_v(surface: &NurbsSurface, u: f64, v: f64, h: f64, v0: f64, v1: f64) -> Vector3<f64> {
    let lo = (v - h).max(v0);
    let hi = (v + h).min(v1);
    let denom = hi - lo;
    if denom <= 0.0 {
        return Vector3::zeros();
    }
    (surface.evaluate(u, hi) - surface.evaluate(u, lo)) / denom
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Planar bicubic surface whose image is exactly the unit square
    /// `[0, 1] × [0, 1]` in the z = 0 plane, with the interior control
    /// points placed at the 1/3, 2/3 lattice so the parameterisation
    /// is the identity `S(u, v) = (u, v, 0)`. Area must be exactly 1.
    fn unit_plane() -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let x = i as f64 / 3.0;
                (0..4)
                    .map(|j| Vector3::new(x, j as f64 / 3.0, 0.0))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    /// Quarter-cylinder of radius `r` and height `h` whose axis is the
    /// y-axis: a rational quadratic quarter-arc in the xz plane (the
    /// standard 3-CP control polygon with middle weight `√2/2`) lofted
    /// linearly along y from 0 to `h`. Its lateral surface is the set
    /// `{ (r cosθ, y, r sinθ) : θ ∈ [0, π/2], y ∈ [0, h] }`.
    fn quarter_cylinder(r: f64, h: f64) -> NurbsSurface {
        let s2 = 2.0_f64.sqrt() / 2.0;
        // u direction (degree 1, linear) = along the y axis.
        let row_y0 = vec![
            Vector3::new(r, 0.0, 0.0),
            Vector3::new(r, 0.0, r),
            Vector3::new(0.0, 0.0, r),
        ];
        let row_yh = vec![
            Vector3::new(r, h, 0.0),
            Vector3::new(r, h, r),
            Vector3::new(0.0, h, r),
        ];
        // v direction (degree 2, rational) = the quarter arc.
        NurbsSurface::new(
            1,
            2,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![row_y0, row_yh],
            vec![vec![1.0, s2, 1.0], vec![1.0, s2, 1.0]],
        )
        .unwrap()
    }

    #[test]
    fn unit_plane_area_is_one() {
        // The parameter Jacobian is constant for a plane, so the
        // midpoint rule is exact: area = 1 to machine precision.
        let a = surface_area(&unit_plane());
        assert!((a - 1.0).abs() < 1.0e-9, "unit plane area = {a}, want 1.0");
    }

    #[test]
    fn unit_plane_centroid_is_square_centre() {
        let res = area_and_centroid(&unit_plane());
        let want = Vector3::new(0.5, 0.5, 0.0);
        assert!(
            (res.centroid - want).norm() < 1.0e-9,
            "unit plane centroid = {:?}, want {want:?}",
            res.centroid
        );
    }

    #[test]
    fn quarter_cylinder_lateral_area_matches_closed_form() {
        // Closed form: lateral area of a quarter cylinder
        //   = arc length × height = (π/2)·r · h.
        let r = 2.0;
        let h = 3.0;
        let a = surface_area(&quarter_cylinder(r, h));
        let expected = (PI / 2.0) * r * h;
        let rel = (a - expected).abs() / expected;
        assert!(
            rel < 1.0e-2,
            "quarter-cylinder area = {a}, want {expected} (rel err {rel})"
        );
    }

    #[test]
    fn quarter_cylinder_centroid_matches_closed_form() {
        // For the lateral surface S(θ, y) = (r cosθ, y, r sinθ) with
        // θ ∈ [0, π/2], y ∈ [0, h] and area element r dθ dy, the
        // area-weighted centroid is
        //   x̄ = z̄ = 2 r / π,   ȳ = h / 2.
        let r = 2.0;
        let h = 3.0;
        let res = area_and_centroid(&quarter_cylinder(r, h));
        let want = Vector3::new(2.0 * r / PI, h / 2.0, 2.0 * r / PI);
        // ~1% tolerance, scaled by the radius for the lateral coords.
        let tol = 1.0e-2 * r;
        assert!(
            (res.centroid - want).norm() < tol,
            "quarter-cylinder centroid = {:?}, want {want:?}",
            res.centroid
        );
    }

    #[test]
    fn area_grows_with_radius_and_height() {
        // Sanity monotonicity: doubling the height doubles the area;
        // doubling the radius doubles the area (lateral area is linear
        // in both for a quarter cylinder).
        let base = surface_area(&quarter_cylinder(1.0, 1.0));
        let taller = surface_area(&quarter_cylinder(1.0, 2.0));
        let wider = surface_area(&quarter_cylinder(2.0, 1.0));
        assert!((taller - 2.0 * base).abs() / (2.0 * base) < 1.0e-2);
        assert!((wider - 2.0 * base).abs() / (2.0 * base) < 1.0e-2);
    }

    #[test]
    fn explicit_sample_count_refines_toward_closed_form() {
        // A coarse grid is already close; a finer grid is no worse.
        // (Monotone convergence is not guaranteed for every surface,
        // but for this smooth case the finer grid should not regress.)
        let s = quarter_cylinder(1.5, 2.5);
        let expected = (PI / 2.0) * 1.5 * 2.5;
        let coarse = area_and_centroid_with_samples(&s, 8).area;
        let fine = area_and_centroid_with_samples(&s, 128).area;
        let err_coarse = (coarse - expected).abs();
        let err_fine = (fine - expected).abs();
        assert!(
            err_fine <= err_coarse + 1.0e-9,
            "refinement regressed: coarse err {err_coarse}, fine err {err_fine}"
        );
        assert!(
            err_fine / expected < 1.0e-3,
            "fine rel err {}",
            err_fine / expected
        );
    }

    #[test]
    fn zero_samples_is_treated_as_one() {
        // Must not divide by zero or panic; a single midpoint cell on
        // the exact plane still integrates to the true area of 1.
        let res = area_and_centroid_with_samples(&unit_plane(), 0);
        assert!((res.area - 1.0).abs() < 1.0e-9, "area = {}", res.area);
    }
}
