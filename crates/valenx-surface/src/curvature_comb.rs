//! **Curvature comb** for Class-A curve fairing — the porcupine/curvature plot
//! a surfacing tool draws to judge whether a curve flows smoothly.
//!
//! A curvature comb samples the curve and, at each sample `P(uᵢ)`, draws a
//! "tooth" of length proportional to the curvature `κ(uᵢ)`, pointing along the
//! curve's principal normal `N` (toward the centre of curvature). The tip
//! points are `P(uᵢ) + scale·κ(uᵢ)·N(uᵢ)`. The envelope of the tooth tips is
//! the comb: a smooth, monotone comb means a fair curve; kinks, spikes or
//! flips in the comb reveal curvature discontinuities the curve itself hides.
//!
//! The curvature magnitude comes from the validated
//! [`NurbsCurve::curvature`](crate::nurbs_curve::NurbsCurve::curvature)
//! (`κ = |C′ × C″| / |C′|³`); the tooth *direction* is the principal normal,
//! the unit acceleration component perpendicular to the tangent, taken from a
//! centred finite-difference stencil consistent with that curvature.
//!
//! Validated against analytic curves: a circular arc of radius `r` gives a
//! uniform comb of height `scale/r` (constant `κ = 1/r`); a straight line gives
//! a flat comb (`κ ≈ 0`).
//!
//! Honest scope: a pointwise curvature-comb diagnostic from the curve's
//! derivatives — research-grade. It is the data behind the plot, not a rendered
//! viewport overlay, and a step toward, not an equal of, CATIA-class curve
//! diagnostics.

use nalgebra::Vector3;

use crate::nurbs_curve::NurbsCurve;

/// A sampled curvature comb: the spine points on the curve, the comb-tooth tip
/// points, and the curvature at each sample.
#[derive(Clone, Debug)]
pub struct CurvatureComb {
    /// Parameters `uᵢ` at which the curve was sampled.
    pub params: Vec<f64>,
    /// Spine points `P(uᵢ)` on the curve.
    pub samples: Vec<Vector3<f64>>,
    /// Comb-tooth tip points `P(uᵢ) + scale·κ(uᵢ)·N(uᵢ)`.
    pub teeth: Vec<Vector3<f64>>,
    /// Curvature `κ(uᵢ)` at each sample.
    pub curvatures: Vec<f64>,
}

impl CurvatureComb {
    /// The largest sampled curvature `max κ(uᵢ)` — the tightest bend.
    pub fn max_curvature(&self) -> f64 {
        self.curvatures.iter().copied().fold(0.0, f64::max)
    }

    /// The parameter `u` at which the largest sampled curvature occurs (the
    /// worst spot to inspect for fairing). Returns `0.0` for an empty comb.
    pub fn max_curvature_param(&self) -> f64 {
        self.params
            .iter()
            .zip(&self.curvatures)
            .fold(
                (0.0, f64::MIN),
                |(bu, bk), (&u, &k)| {
                    if k > bk {
                        (u, k)
                    } else {
                        (bu, bk)
                    }
                },
            )
            .0
    }
}

/// Build a [`CurvatureComb`] by sampling `curve` at `samples` evenly-spaced
/// parameters across its valid range. Each tooth has length `scale·κ` and
/// points along the principal normal (toward the centre of curvature).
///
/// `samples` is clamped to at least 2; `scale` sets the visual exaggeration of
/// the teeth (it does not affect [`CurvatureComb::curvatures`] or
/// [`CurvatureComb::max_curvature`], which are the true curvature).
pub fn curvature_comb(curve: &NurbsCurve, samples: usize, scale: f64) -> CurvatureComb {
    let n = samples.max(2);
    let (u_min, u_max) = curve.parameter_range();
    let mut params = Vec::with_capacity(n);
    let mut pts = Vec::with_capacity(n);
    let mut teeth = Vec::with_capacity(n);
    let mut curvatures = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let u = u_min + t * (u_max - u_min);
        let p = curve.evaluate(u);
        let k = curve.curvature(u);
        let normal = principal_normal(curve, u);
        let tip = p + normal * (scale * k);
        params.push(u);
        pts.push(p);
        teeth.push(tip);
        curvatures.push(k);
    }
    CurvatureComb {
        params,
        samples: pts,
        teeth,
        curvatures,
    }
}

/// Unit **principal normal** `N(u)` — the direction of the acceleration
/// component perpendicular to the velocity, which points toward the centre of
/// curvature. Uses the same centred 3-point stencil as
/// [`NurbsCurve::curvature`](crate::nurbs_curve::NurbsCurve::curvature) so the
/// tooth direction is consistent with its length. Returns the zero vector at a
/// straight or singular point (where the principal normal is undefined).
fn principal_normal(curve: &NurbsCurve, u: f64) -> Vector3<f64> {
    let (u_min, u_max) = curve.parameter_range();
    let span = u_max - u_min;
    let h = (span * 1e-3).max(1e-6);
    if span < 4.0 * h {
        return Vector3::zeros();
    }
    let uc = u.max(u_min + h).min(u_max - h);
    let c_minus = curve.evaluate(uc - h);
    let c_0 = curve.evaluate(uc);
    let c_plus = curve.evaluate(uc + h);
    let d1 = (c_plus - c_minus) / (2.0 * h);
    let d2 = (c_plus - 2.0 * c_0 + c_minus) / (h * h);
    let speed = d1.norm();
    if speed < 1e-12 {
        return Vector3::zeros();
    }
    let tangent = d1 / speed;
    // Acceleration component perpendicular to the tangent → toward the centre
    // of curvature.
    let perp = d2 - tangent * d2.dot(&tangent);
    let pn = perp.norm();
    if pn < 1e-12 {
        Vector3::zeros()
    } else {
        perp / pn
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_1_SQRT_2;

    /// A rational-quadratic quarter circle of radius `r` in the XY plane,
    /// centred at the origin (exact circular arc, `κ ≡ 1/r`).
    fn quarter_circle(r: f64) -> NurbsCurve {
        NurbsCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Vector3::new(r, 0.0, 0.0),
                Vector3::new(r, r, 0.0),
                Vector3::new(0.0, r, 0.0),
            ],
            vec![1.0, FRAC_1_SQRT_2, 1.0],
        )
        .expect("valid quarter circle")
    }

    /// A straight degree-1 segment from `a` to `b`.
    fn segment(a: Vector3<f64>, b: Vector3<f64>) -> NurbsCurve {
        NurbsCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![a, b], vec![1.0, 1.0])
            .expect("valid segment")
    }

    #[test]
    fn circle_has_uniform_comb_of_height_scale_over_r() {
        let r = 2.5;
        let scale = 0.4;
        let comb = curvature_comb(&quarter_circle(r), 12, scale);
        let kappa = 1.0 / r;
        for (i, &k) in comb.curvatures.iter().enumerate() {
            assert!(
                (k - kappa).abs() < 0.02 * kappa,
                "sample {i}: κ {k} vs analytic 1/r {kappa}"
            );
            // Comb height = |tooth − spine| = scale·κ, uniform around the arc.
            let height = (comb.teeth[i] - comb.samples[i]).norm();
            let expected = scale * kappa;
            assert!(
                (height - expected).abs() < 0.02 * expected,
                "sample {i}: comb height {height} vs scale/r {expected}"
            );
        }
        assert!((comb.max_curvature() - kappa).abs() < 0.02 * kappa);
    }

    #[test]
    fn circle_teeth_point_toward_the_centre() {
        // The principal normal points to the centre of curvature — for a
        // circle centred at the origin, every tooth tip is closer to the
        // origin than its spine point.
        let r = 3.0;
        let comb = curvature_comb(&quarter_circle(r), 9, 0.5);
        for i in 0..comb.samples.len() {
            let spine_d = comb.samples[i].norm();
            let tip_d = comb.teeth[i].norm();
            assert!(
                tip_d < spine_d,
                "sample {i}: tip dist {tip_d} should be < spine dist {spine_d}"
            );
        }
    }

    #[test]
    fn straight_line_has_a_flat_comb() {
        let comb = curvature_comb(
            &segment(Vector3::new(-1.0, 0.5, 0.0), Vector3::new(3.0, 0.5, 0.0)),
            10,
            1.0,
        );
        for (i, &k) in comb.curvatures.iter().enumerate() {
            assert!(k.abs() < 1e-6, "sample {i}: line κ {k} should be ~0");
            // Flat comb: teeth coincide with the spine.
            let height = (comb.teeth[i] - comb.samples[i]).norm();
            assert!(height < 1e-6, "sample {i}: line comb height {height} ≠ 0");
        }
        assert!(comb.max_curvature() < 1e-6);
    }

    #[test]
    fn max_curvature_param_lands_in_range() {
        let comb = curvature_comb(&quarter_circle(1.0), 8, 1.0);
        let (u_min, u_max) = quarter_circle(1.0).parameter_range();
        let up = comb.max_curvature_param();
        assert!(up >= u_min && up <= u_max, "max-κ param {up} out of range");
    }
}
