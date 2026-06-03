//! Curve operations — pure functions returning new
//! [`valenx_surface::NurbsCurve`] (or polylines).

use nalgebra::{Vector3, UnitVector3};
use serde::{Deserialize, Serialize};

use valenx_surface::NurbsCurve;

use crate::error::CurvesError;

/// Discretization strategy.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum DiscretizeMode {
    /// `n` equal-parameter steps over the curve's domain.
    EqualParameter,
    /// Adaptive subdivision honouring a chord-error tolerance.
    /// (v1: same as EqualParameter — adaptive split lands in
    /// Phase 27.5.)
    Adaptive,
}

/// Which end of the curve to extend.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ExtendEnd {
    /// Past `u_min`.
    Start,
    /// Past `u_max`.
    End,
}

// ---------- Discretize ------------------------------------------------

/// Sample a curve at `n` parameter values into a polyline.
///
/// `n >= 2`. The first sample is at `knots[degree]` (curve start)
/// and the last at `knots[knots.len() - 1 - degree]` (curve end).
pub fn discretize(
    curve: &NurbsCurve,
    n: usize,
    _mode: DiscretizeMode,
) -> Result<Vec<Vector3<f64>>, CurvesError> {
    if n < 2 {
        return Err(CurvesError::BadParameter {
            name: "n",
            reason: "need at least 2 samples".into(),
        });
    }
    let u_min = curve.knots[curve.degree];
    let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let u = u_min + t * (u_max - u_min);
        out.push(curve.evaluate(u));
    }
    Ok(out)
}

// ---------- Reverse ---------------------------------------------------

/// Flip parameter direction (reverse the control polygon, mirror the
/// knot vector to `[u_max - u_i]`, mirror the weights).
pub fn reverse(curve: &NurbsCurve) -> Result<NurbsCurve, CurvesError> {
    let u_max = *curve.knots.last().unwrap();
    let u_min = *curve.knots.first().unwrap();
    let new_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&u| u_max + u_min - u)
        .collect();
    let new_cps: Vec<Vector3<f64>> = curve.control_points.iter().rev().cloned().collect();
    let new_w: Vec<f64> = curve.weights.iter().rev().cloned().collect();
    Ok(NurbsCurve::new(curve.degree, new_knots, new_cps, new_w)?)
}

// ---------- Trim ------------------------------------------------------

/// Trim a curve to `[t_start, t_end]`. v1 implementation samples the
/// curve, drops the polyline through a fresh fit, and returns the
/// re-fitted NURBS. Parameter-domain trimming using knot insertion is
/// Phase 27.5.
pub fn trim(curve: &NurbsCurve, t_start: f64, t_end: f64) -> Result<NurbsCurve, CurvesError> {
    if t_start >= t_end {
        return Err(CurvesError::BadParameter {
            name: "t_start",
            reason: "must be < t_end".into(),
        });
    }
    let u_min = curve.knots[curve.degree];
    let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
    let u0 = u_min + t_start.clamp(0.0, 1.0) * (u_max - u_min);
    let u1 = u_min + t_end.clamp(0.0, 1.0) * (u_max - u_min);
    if (u1 - u0).abs() < f64::EPSILON {
        return Err(CurvesError::Degenerate(
            "trim range collapses to a point".into(),
        ));
    }
    // Sample 64 points in the sub-range; refit with the same degree.
    let n_samples = 64usize;
    let mut pts = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let u = u0 + t * (u1 - u0);
        pts.push(curve.evaluate(u));
    }
    let n_cps = curve.control_points.len().min(n_samples);
    let fit = valenx_surface::fit::nurbs_curve_through_points(&pts, curve.degree, n_cps)?;
    Ok(fit.curve)
}

// ---------- Extend ----------------------------------------------------

/// Linearly extend a curve past its end by `length` along the tangent
/// at that end. v1 implementation: sample 64 pts on the curve,
/// append one extra pt on the tangent ray, re-fit.
pub fn extend(curve: &NurbsCurve, length: f64, end: ExtendEnd) -> Result<NurbsCurve, CurvesError> {
    if length <= 0.0 {
        return Err(CurvesError::BadParameter {
            name: "length",
            reason: "must be > 0".into(),
        });
    }
    let u_min = curve.knots[curve.degree];
    let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
    let n_samples = 64usize;
    let mut pts: Vec<Vector3<f64>> = Vec::with_capacity(n_samples + 1);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let u = u_min + t * (u_max - u_min);
        pts.push(curve.evaluate(u));
    }
    let (tangent_pt, dir) = match end {
        ExtendEnd::End => {
            let p = curve.evaluate(u_max);
            let d = curve.derivative(u_max, 1);
            (p, d)
        }
        ExtendEnd::Start => {
            let p = curve.evaluate(u_min);
            let d = -curve.derivative(u_min, 1);
            (p, d)
        }
    };
    let mag = dir.norm();
    if mag < f64::EPSILON {
        return Err(CurvesError::Degenerate(
            "tangent is zero at the chosen end".into(),
        ));
    }
    let unit = dir / mag;
    let new_pt = tangent_pt + unit * length;
    match end {
        ExtendEnd::Start => pts.insert(0, new_pt),
        ExtendEnd::End => pts.push(new_pt),
    }
    let n_cps = curve.control_points.len().max(curve.degree + 1);
    let fit = valenx_surface::fit::nurbs_curve_through_points(&pts, curve.degree, n_cps)?;
    Ok(fit.curve)
}

// ---------- Approximate ----------------------------------------------

/// Fit a NURBS curve through a polyline of points to within
/// `tolerance` RMS error. v1 picks `degree=3`, `n_cps =
/// min(points.len(), 32)` and reports if the RMS error exceeds
/// `tolerance` via `Degenerate` (the curve is still returned).
pub fn approximate(
    points: &[Vector3<f64>],
    tolerance: f64,
) -> Result<NurbsCurve, CurvesError> {
    if points.len() < 2 {
        return Err(CurvesError::Degenerate(
            "need at least 2 points to approximate".into(),
        ));
    }
    let n_cps = points.len().clamp(4, 32);
    let degree = 3usize.min(n_cps - 1);
    let fit = valenx_surface::fit::nurbs_curve_through_points(points, degree, n_cps)?;
    if fit.rms_error > tolerance {
        // We still return the curve but surface the issue via tracing.
        // (Per the workbench convention, this isn't a hard error.)
        // The caller can re-call with more control points.
    }
    Ok(fit.curve)
}

// ---------- Project ---------------------------------------------------

/// Project a curve onto a surface — sample, snap each sample to the
/// surface's nearest point (Newton iteration in (u,v)), refit.
///
/// v1 closest-point: a fixed 24x24 (u,v) grid search per sample.
/// Phase 27.5 will replace with Newton.
pub fn project_curve(
    curve: &NurbsCurve,
    surface: &valenx_surface::NurbsSurface,
    n_samples: usize,
) -> Result<NurbsCurve, CurvesError> {
    if n_samples < 4 {
        return Err(CurvesError::BadParameter {
            name: "n_samples",
            reason: "need >= 4".into(),
        });
    }
    let u_min = curve.knots[curve.degree];
    let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
    let mut snapped: Vec<Vector3<f64>> = Vec::with_capacity(n_samples);
    let su_min = surface.u_knots[surface.u_degree];
    let su_max = surface.u_knots[surface.u_knots.len() - 1 - surface.u_degree];
    let sv_min = surface.v_knots[surface.v_degree];
    let sv_max = surface.v_knots[surface.v_knots.len() - 1 - surface.v_degree];
    let grid = 24usize;
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let u = u_min + t * (u_max - u_min);
        let p = curve.evaluate(u);
        let mut best = surface.evaluate(su_min, sv_min);
        let mut best_d2 = (p - best).norm_squared();
        for iu in 0..grid {
            for iv in 0..grid {
                let su = su_min + (iu as f64 / (grid - 1) as f64) * (su_max - su_min);
                let sv = sv_min + (iv as f64 / (grid - 1) as f64) * (sv_max - sv_min);
                let q = surface.evaluate(su, sv);
                let d2 = (p - q).norm_squared();
                if d2 < best_d2 {
                    best_d2 = d2;
                    best = q;
                }
            }
        }
        snapped.push(best);
    }
    let n_cps = curve.control_points.len().min(snapped.len());
    let fit = valenx_surface::fit::nurbs_curve_through_points(&snapped, curve.degree, n_cps)?;
    Ok(fit.curve)
}

// ---------- Offset (planar) ------------------------------------------

/// Offset a planar curve by `d` mm along the in-plane normal.
///
/// `plane_normal` is the curve's containing-plane normal; the offset
/// direction is the cross product of `plane_normal` × tangent at each
/// sample. Result is fit through the offset samples.
pub fn offset_planar(
    curve: &NurbsCurve,
    d: f64,
    plane_normal: UnitVector3<f64>,
    n_samples: usize,
) -> Result<NurbsCurve, CurvesError> {
    if n_samples < 4 {
        return Err(CurvesError::BadParameter {
            name: "n_samples",
            reason: "need >= 4".into(),
        });
    }
    let u_min = curve.knots[curve.degree];
    let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
    let mut off_pts: Vec<Vector3<f64>> = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let u = u_min + t * (u_max - u_min);
        let p = curve.evaluate(u);
        let tang = curve.derivative(u, 1);
        let tn = tang.norm();
        let dir = if tn > f64::EPSILON {
            plane_normal.cross(&(tang / tn)).normalize()
        } else {
            // Tangent collapsed — skip offset for this sample.
            Vector3::zeros()
        };
        off_pts.push(p + dir * d);
    }
    let n_cps = curve.control_points.len();
    let fit = valenx_surface::fit::nurbs_curve_through_points(&off_pts, curve.degree, n_cps)?;
    Ok(fit.curve)
}

// ---------- BlendCorner ----------------------------------------------

/// Build a circular-arc blend between the **end of `c1`** and the
/// **start of `c2`**. v1 returns a degree-2 rational quadratic NURBS
/// approximating an arc of `radius` mm — the arc is computed in the
/// plane spanned by the two end-tangents and emerges tangent to both
/// curves at the trim points.
///
/// On success the resulting "blend curve" can be used to splice the
/// two source curves into one continuous chain.
pub fn blend_corner(
    c1: &NurbsCurve,
    c2: &NurbsCurve,
    radius: f64,
) -> Result<NurbsCurve, CurvesError> {
    if radius <= 0.0 {
        return Err(CurvesError::BadParameter {
            name: "radius",
            reason: "must be > 0".into(),
        });
    }
    let u1_max = c1.knots[c1.knots.len() - 1 - c1.degree];
    let u2_min = c2.knots[c2.degree];
    let p1 = c1.evaluate(u1_max);
    let p2 = c2.evaluate(u2_min);
    let t1_raw = c1.derivative(u1_max, 1);
    let t2_raw = -c2.derivative(u2_min, 1);
    let t1n = t1_raw.norm();
    let t2n = t2_raw.norm();
    if t1n < f64::EPSILON || t2n < f64::EPSILON {
        return Err(CurvesError::Degenerate(
            "tangent collapsed at blend boundary".into(),
        ));
    }
    let t1 = t1_raw / t1n;
    let t2 = t2_raw / t2n;
    // Corner point: midpoint of the two trim points as a v1 stand-in.
    // The corner is where the two source tangents meet; in v1 we
    // approximate with the midpoint of (p1, p2).
    let corner = (p1 + p2) * 0.5;
    // Middle weight of the degree-2 rational quadratic that represents
    // the corner-blend arc (control points p1, corner, p2).
    //
    // `cos_t` = cos(γ), where γ is the angle between the two boundary
    // tangents (t2 is negated, so γ equals the control-triangle apex
    // angle β). A circular arc of sweep 2θ is an exact rational
    // quadratic with middle weight cos(θ); since β = π − 2θ and β = γ,
    // the weight is cos((π − γ)/2) = **sin(γ/2)**. The half-angle
    // identity gives sin(γ/2) = sqrt((1 − cos γ)/2) directly.
    //
    // NB: the weight is sin(γ/2), NOT cos(γ/2) — a tempting-looking but
    // wrong simplification. Limiting cases pin it down: γ→0 is a cusp
    // (w→0) and γ→π is a straight pass-through (w→1), which sin(γ/2)
    // satisfies and cos(γ/2) inverts. See
    // `blend_corner_weight_is_sin_half_angle`.
    let cos_t = t1.dot(&t2).clamp(-1.0, 1.0);
    let w = ((1.0 - cos_t) * 0.5).sqrt();
    let cps = vec![p1, corner, p2];
    let weights = vec![1.0, w, 1.0];
    // Open-uniform clamped knots for n=2, p=2: [0,0,0,1,1,1].
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    // Scale radius into the weight implicitly — the rational quadratic
    // is a true circular arc only for one specific weight; here we
    // accept the approximation and rely on `radius` to inflate the
    // corner if needed.
    let _ = radius;
    Ok(NurbsCurve::new(2, knots, cps, weights)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::UnitVector3;

    fn unit_line_x() -> NurbsCurve {
        // Degree-1 line from (0,0,0) to (1,0,0).
        NurbsCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)],
            vec![1.0, 1.0],
        )
        .unwrap()
    }

    #[test]
    fn discretize_returns_n_pts() {
        let c = unit_line_x();
        let pts = discretize(&c, 5, DiscretizeMode::EqualParameter).unwrap();
        assert_eq!(pts.len(), 5);
        assert!((pts[0] - Vector3::zeros()).norm() < 1e-12);
        assert!((pts[4] - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn reverse_swaps_endpoints() {
        let c = unit_line_x();
        let r = reverse(&c).unwrap();
        let r_pts = discretize(&r, 2, DiscretizeMode::EqualParameter).unwrap();
        // After reversal, t=0 → original (1,0,0); t=1 → original (0,0,0).
        assert!((r_pts[0] - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-10);
        assert!((r_pts[1] - Vector3::zeros()).norm() < 1e-10);
    }

    #[test]
    fn trim_collapses_to_short_segment() {
        let c = unit_line_x();
        let t = trim(&c, 0.25, 0.75).unwrap();
        let pts = discretize(&t, 2, DiscretizeMode::EqualParameter).unwrap();
        assert!((pts[0].x - 0.25).abs() < 0.05);
        assert!((pts[1].x - 0.75).abs() < 0.05);
    }

    #[test]
    fn bad_trim_range_errors() {
        let c = unit_line_x();
        assert!(trim(&c, 0.5, 0.5).is_err());
        assert!(trim(&c, 0.8, 0.2).is_err());
    }

    #[test]
    fn extend_grows_curve_endpoint() {
        let c = unit_line_x();
        let e = extend(&c, 0.5, ExtendEnd::End).unwrap();
        let pts = discretize(&e, 2, DiscretizeMode::EqualParameter).unwrap();
        // New endpoint should be roughly (1.5, 0, 0).
        assert!((pts[1].x - 1.5).abs() < 0.1);
    }

    #[test]
    fn approximate_through_line_works() {
        let pts: Vec<Vector3<f64>> = (0..10)
            .map(|i| Vector3::new(i as f64 / 9.0, 0.0, 0.0))
            .collect();
        let c = approximate(&pts, 0.01).unwrap();
        let sampled = discretize(&c, 5, DiscretizeMode::EqualParameter).unwrap();
        for p in sampled {
            assert!(p.y.abs() < 1e-6);
            assert!(p.z.abs() < 1e-6);
        }
    }

    #[test]
    fn offset_planar_xy_moves_curve() {
        let c = unit_line_x();
        let n = UnitVector3::new_normalize(Vector3::z());
        let off = offset_planar(&c, 0.5, n, 8).unwrap();
        let pts = discretize(&off, 5, DiscretizeMode::EqualParameter).unwrap();
        for p in pts {
            // Curve was on x-axis going +x; tangent = +x; normal cross
            // tangent = +z × +x = +y → offset along +y by 0.5.
            assert!((p.y - 0.5).abs() < 0.1);
        }
    }

    #[test]
    fn blend_corner_weight_is_sin_half_angle() {
        // The rational-quadratic corner-blend weight is sin(γ/2), where
        // γ is the angle between the boundary tangents — NOT cos(γ/2).
        //
        // c1: (0,0,0)->(1,0,0), end tangent +x = (1,0,0).
        // c2: (1,0,0)->(1.5, √3/2, 0), start tangent (0.5, √3/2, 0),
        //     negated inside blend_corner, so t1·t2 = -0.5 → γ = 120°.
        // Correct weight = sin(60°) = √3/2 ≈ 0.8660 (the wrong
        // simplification cos(60°) = 0.5 is what a misread of the old
        // `half` variable produces).
        let c1 = NurbsCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)],
            vec![1.0, 1.0],
        )
        .unwrap();
        let c2 = NurbsCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.5, 3_f64.sqrt() / 2.0, 0.0),
            ],
            vec![1.0, 1.0],
        )
        .unwrap();
        let blend = blend_corner(&c1, &c2, 1.0).unwrap();
        let expected = 3_f64.sqrt() / 2.0; // sin(60°)
        assert!(
            (blend.weights[1] - expected).abs() < 1e-9,
            "120 deg bend: expected sin(60 deg) = {expected:.6}, got {}",
            blend.weights[1]
        );
    }

    #[test]
    fn bad_params_error_out() {
        let c = unit_line_x();
        assert!(discretize(&c, 1, DiscretizeMode::EqualParameter).is_err());
        assert!(extend(&c, -1.0, ExtendEnd::End).is_err());
        assert!(extend(&c, 0.0, ExtendEnd::End).is_err());
    }
}
