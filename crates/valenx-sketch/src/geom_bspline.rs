//! 2-D B-Spline primitive for the sketcher.
//!
//! Phase 12A — same Cox-de Boor recursion as `valenx-surface::NurbsCurve`
//! but on 2-D control points (no rational weights in v1; all weights = 1).
//! Control points live as variable-indexed pairs in [`crate::sketch::Sketch::vars`]
//! so the solver can drive them with constraints just like other primitives.

use serde::{Deserialize, Serialize};

use crate::geom::Point2;

/// 2-D B-Spline curve.
///
/// Control points are stored as [`Point2`]s into the sketch parameter
/// vector. The knot vector and degree are stored by value (they are
/// not solver-variable for v1 — degree/knot edits go through a
/// rebuild rather than a constraint).
///
/// `PartialEq` is derived structurally — the float-valued `knots` and
/// `weights` use IEEE 754 semantics so a NaN never compares equal to
/// itself (curves with NaN params fail to dedupe in undo snapshots).
/// The sketcher never lets a NaN reach these fields in practice.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BSpline2 {
    /// Polynomial degree `p` (1 = polyline, 3 = cubic, etc.).
    pub degree: usize,
    /// Knot vector — non-decreasing, length `n_cp + degree + 1`.
    pub knots: Vec<f64>,
    /// Control points — variable-indexed 2-D points.
    pub control_points: Vec<Point2>,
    /// Per-control-point weights. Same length as `control_points`.
    /// All-ones reduces to a non-rational B-spline.
    pub weights: Vec<f64>,
}

impl BSpline2 {
    /// Number of control points.
    pub fn n_control_points(&self) -> usize {
        self.control_points.len()
    }

    /// Valid parameter range: `[knots[degree], knots[n]]`.
    pub fn parameter_range(&self) -> (f64, f64) {
        let n = self.control_points.len();
        (self.knots[self.degree], self.knots[n])
    }

    /// Evaluate the B-spline at parameter `u`. Returns `[x, y]`.
    pub fn evaluate(&self, vars: &[f64], u: f64) -> [f64; 2] {
        let span = find_knot_span(u, &self.knots, self.degree, self.control_points.len());
        let basis = basis_functions(span, u, self.degree, &self.knots);
        let mut nx = 0.0_f64;
        let mut ny = 0.0_f64;
        let mut den = 0.0_f64;
        for (i, b) in basis.iter().enumerate() {
            let cp_idx = span - self.degree + i;
            let w = self.weights[cp_idx];
            let wb = w * b;
            let (px, py) = self.control_points[cp_idx].read(vars);
            nx += px * wb;
            ny += py * wb;
            den += wb;
        }
        if den.abs() < 1e-30 {
            [nx, ny]
        } else {
            [nx / den, ny / den]
        }
    }

    /// First derivative at `u` via central finite differences.
    pub fn derivative(&self, vars: &[f64], u: f64) -> [f64; 2] {
        let (u_min, u_max) = self.parameter_range();
        let h = 1e-4_f64.max((u_max - u_min) * 1e-5);
        let u_lo = (u - h).max(u_min);
        let u_hi = (u + h).min(u_max);
        let denom = u_hi - u_lo;
        if denom.abs() < 1e-30 {
            return [0.0, 0.0];
        }
        let lo = self.evaluate(vars, u_lo);
        let hi = self.evaluate(vars, u_hi);
        [(hi[0] - lo[0]) / denom, (hi[1] - lo[1]) / denom]
    }

    /// Start and end points (parameter min and max).
    pub fn endpoints(&self, vars: &[f64]) -> ([f64; 2], [f64; 2]) {
        let (u_min, u_max) = self.parameter_range();
        (self.evaluate(vars, u_min), self.evaluate(vars, u_max))
    }

    /// Closest parameter `u` on the curve to a 2-D point.
    ///
    /// ## Algorithm (12.5 — multi-seed safeguarded Newton)
    ///
    /// The closest-point problem is `min_u |P(u) − target|²`, whose
    /// stationarity condition is `f(u) = (P(u) − target) · P'(u) = 0`.
    /// `f` can have several roots (a wiggly curve has several local
    /// extrema of distance), so a single Newton run from one seed can
    /// land on a far local minimum. This routine:
    ///
    /// 1. **Coarse scan** at `4·n_cp` (≥ 32) uniform samples, recording
    ///    *every* local minimum of the squared distance — each is a
    ///    basin that may contain the global closest point.
    /// 2. From each local-minimum seed runs a **safeguarded Newton**
    ///    that brackets the root of `f` and falls back to bisection
    ///    whenever a Newton step would leave the bracket or fail to
    ///    decrease `|f|`. The bracket guarantees convergence even when
    ///    `f'` is small or wrong-signed — the "adaptive bracketing"
    ///    the Phase 12 v1 lacked.
    /// 3. Returns the seed whose refined parameter gives the smallest
    ///    distance globally.
    pub fn closest_param(&self, vars: &[f64], target: [f64; 2]) -> f64 {
        let (u_min, u_max) = self.parameter_range();
        if u_max - u_min < 1e-30 {
            return u_min;
        }
        let dist2 = |u: f64| {
            let p = self.evaluate(vars, u);
            let dx = p[0] - target[0];
            let dy = p[1] - target[1];
            dx * dx + dy * dy
        };

        // --- Step 1: coarse scan, collecting every local minimum. ---
        let n_samples = (4 * self.control_points.len()).max(32);
        let sample_u = |i: usize| u_min + (u_max - u_min) * (i as f64 / n_samples as f64);
        let d2: Vec<f64> = (0..=n_samples).map(|i| dist2(sample_u(i))).collect();

        let mut seeds: Vec<f64> = Vec::new();
        for i in 0..=n_samples {
            let here = d2[i];
            let left = if i == 0 { f64::INFINITY } else { d2[i - 1] };
            let right = if i == n_samples {
                f64::INFINITY
            } else {
                d2[i + 1]
            };
            // A local minimum (endpoints count when they beat their
            // single neighbour).
            if here <= left && here <= right {
                seeds.push(sample_u(i));
            }
        }
        if seeds.is_empty() {
            // Degenerate flat scan — seed the global-min sample.
            let best = (0..=n_samples)
                .min_by(|&a, &b| d2[a].partial_cmp(&d2[b]).unwrap())
                .unwrap_or(0);
            seeds.push(sample_u(best));
        }

        // --- Step 2+3: refine each seed, keep the global best. ---
        let mut best_u = seeds[0];
        let mut best_d2 = dist2(best_u);
        for &seed in &seeds {
            let refined = self.refine_closest(vars, target, seed, u_min, u_max);
            let rd2 = dist2(refined);
            if rd2 < best_d2 {
                best_d2 = rd2;
                best_u = refined;
            }
        }
        best_u
    }

    /// Safeguarded Newton-bisection refinement of one closest-point
    /// seed. Brackets the root of `f(u) = (P(u) − target) · P'(u)`
    /// around `seed` and converges with a guaranteed-progress hybrid.
    fn refine_closest(
        &self,
        vars: &[f64],
        target: [f64; 2],
        seed: f64,
        u_min: f64,
        u_max: f64,
    ) -> f64 {
        // f(u) = (P − target) · P'  (zero at a distance extremum).
        let f = |u: f64| {
            let p = self.evaluate(vars, u);
            let d = self.derivative(vars, u);
            (p[0] - target[0]) * d[0] + (p[1] - target[1]) * d[1]
        };
        // Establish a bracket [lo, hi] around `seed` where f changes
        // sign — expand outward from the seed by a fraction of the span.
        let span = u_max - u_min;
        let mut lo = seed;
        let mut hi = seed;
        let mut f_lo = f(seed);
        let mut f_hi = f_lo;
        let mut step = span / 64.0;
        let mut bracketed = false;
        for _ in 0..8 {
            let new_lo = (lo - step).max(u_min);
            let new_hi = (hi + step).min(u_max);
            f_lo = f(new_lo);
            f_hi = f(new_hi);
            lo = new_lo;
            hi = new_hi;
            if f_lo * f_hi <= 0.0 {
                bracketed = true;
                break;
            }
            step *= 2.0;
        }
        if !bracketed {
            // No sign change found near the seed — `f` keeps one sign,
            // so the extremum is at whichever bracket end has the
            // smaller |f| (the curve is monotone in distance here).
            return if f_lo.abs() <= f_hi.abs() { lo } else { hi };
        }
        // Orient the bracket so f(lo) < 0 <= f(hi).
        if f_lo > 0.0 {
            std::mem::swap(&mut lo, &mut hi);
            std::mem::swap(&mut f_lo, &mut f_hi);
        }
        let mut u = 0.5 * (lo + hi);
        for _ in 0..40 {
            let fu = f(u);
            if fu.abs() < 1e-13 {
                break;
            }
            // Shrink the bracket using the sign of f(u).
            if fu < 0.0 {
                lo = u;
            } else {
                hi = u;
            }
            // Newton step on f, with f' by central difference.
            let h = (span * 1e-5).max(1e-7);
            let fp = (f((u + h).min(u_max)) - f((u - h).max(u_min)))
                / (2.0 * h).max(1e-30);
            let newton = if fp.abs() > 1e-14 {
                u - fu / fp
            } else {
                f64::NAN
            };
            // Accept Newton only if it stays strictly inside the
            // bracket; otherwise bisect. This is the safeguard.
            let bracket_lo = lo.min(hi);
            let bracket_hi = lo.max(hi);
            let next = if newton.is_finite()
                && newton > bracket_lo
                && newton < bracket_hi
            {
                newton
            } else {
                0.5 * (lo + hi)
            };
            if (next - u).abs() < 1e-12 {
                u = next;
                break;
            }
            u = next;
        }
        u.clamp(u_min, u_max)
    }
}

/// Find the knot span `k` such that `knots[k] <= u < knots[k+1]`.
pub fn find_knot_span(u: f64, knots: &[f64], degree: usize, n_cp: usize) -> usize {
    let n = n_cp;
    if u >= knots[n] {
        return n - 1;
    }
    if u <= knots[degree] {
        return degree;
    }
    for k in degree..n {
        if knots[k] <= u && u < knots[k + 1] {
            return k;
        }
    }
    n - 1
}

/// Cox-de Boor basis functions for span `k`. Returns `degree + 1` values.
pub fn basis_functions(span: usize, u: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
    let p = degree;
    let mut n = vec![0.0_f64; p + 1];
    let mut left = vec![0.0_f64; p + 1];
    let mut right = vec![0.0_f64; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() < 1e-30 {
                0.0
            } else {
                n[r] / denom
            };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::Sketch;

    fn bezier_knots() -> Vec<f64> {
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    }

    #[test]
    fn parameter_range_for_clamped_cubic() {
        let mut s = Sketch::new();
        let p0 = s.add_point(0.0, 0.0);
        let p1 = s.add_point(1.0, 1.0);
        let p2 = s.add_point(2.0, 1.0);
        let p3 = s.add_point(3.0, 0.0);
        let cps = vec![
            s.point_at(p0).unwrap(),
            s.point_at(p1).unwrap(),
            s.point_at(p2).unwrap(),
            s.point_at(p3).unwrap(),
        ];
        let curve = BSpline2 {
            degree: 3,
            knots: bezier_knots(),
            control_points: cps,
            weights: vec![1.0; 4],
        };
        assert_eq!(curve.parameter_range(), (0.0, 1.0));
    }

    #[test]
    fn clamped_endpoints_match_first_and_last_control_point() {
        let mut s = Sketch::new();
        let p0 = s.add_point(0.0, 0.0);
        let p1 = s.add_point(1.0, 1.0);
        let p2 = s.add_point(2.0, 1.0);
        let p3 = s.add_point(3.0, 0.0);
        let cps = vec![
            s.point_at(p0).unwrap(),
            s.point_at(p1).unwrap(),
            s.point_at(p2).unwrap(),
            s.point_at(p3).unwrap(),
        ];
        let curve = BSpline2 {
            degree: 3,
            knots: bezier_knots(),
            control_points: cps,
            weights: vec![1.0; 4],
        };
        let (start, end) = curve.endpoints(&s.vars);
        assert!((start[0] - 0.0).abs() < 1e-10);
        assert!((start[1] - 0.0).abs() < 1e-10);
        assert!((end[0] - 3.0).abs() < 1e-10);
        assert!((end[1] - 0.0).abs() < 1e-10);
    }

    /// Task 2: degree-3 Bezier midpoint matches de Casteljau.
    /// For a cubic Bezier with control points P0..P3 at u=0.5,
    /// de Casteljau gives: (P0 + 3P1 + 3P2 + P3) / 8.
    #[test]
    fn cubic_bezier_midpoint_matches_de_casteljau() {
        let mut s = Sketch::new();
        let p0 = s.add_point(0.0, 0.0);
        let p1 = s.add_point(1.0, 2.0);
        let p2 = s.add_point(3.0, 2.0);
        let p3 = s.add_point(4.0, 0.0);
        let cps = vec![
            s.point_at(p0).unwrap(),
            s.point_at(p1).unwrap(),
            s.point_at(p2).unwrap(),
            s.point_at(p3).unwrap(),
        ];
        let curve = BSpline2 {
            degree: 3,
            knots: bezier_knots(),
            control_points: cps,
            weights: vec![1.0; 4],
        };
        let mid = curve.evaluate(&s.vars, 0.5);
        let expected_x = (0.0 + 3.0 * 1.0 + 3.0 * 3.0 + 4.0) / 8.0; // 2.0
        let expected_y = (0.0 + 3.0 * 2.0 + 3.0 * 2.0 + 0.0) / 8.0; // 1.5
        assert!(
            (mid[0] - expected_x).abs() < 1e-10,
            "x: got {} want {}",
            mid[0],
            expected_x
        );
        assert!(
            (mid[1] - expected_y).abs() < 1e-10,
            "y: got {} want {}",
            mid[1],
            expected_y
        );
    }

    /// Build a degree-3 Bezier from four explicit control points.
    fn cubic_bezier(s: &mut Sketch, pts: [[f64; 2]; 4]) -> BSpline2 {
        let cps: Vec<Point2> = pts
            .iter()
            .map(|p| {
                let id = s.add_point(p[0], p[1]);
                s.point_at(id).unwrap()
            })
            .collect();
        BSpline2 {
            degree: 3,
            knots: bezier_knots(),
            control_points: cps,
            weights: vec![1.0; 4],
        }
    }

    #[test]
    fn closest_param_on_a_point_already_on_the_curve() {
        // 12.5: the closest parameter to a point lying exactly on the
        // curve must reproduce that point.
        let mut s = Sketch::new();
        let curve = cubic_bezier(
            &mut s,
            [[0.0, 0.0], [1.0, 2.0], [3.0, 2.0], [4.0, 0.0]],
        );
        let on_curve = curve.evaluate(&s.vars, 0.37);
        let u = curve.closest_param(&s.vars, on_curve);
        let recovered = curve.evaluate(&s.vars, u);
        let err = ((recovered[0] - on_curve[0]).powi(2)
            + (recovered[1] - on_curve[1]).powi(2))
        .sqrt();
        assert!(err < 1e-7, "closest point off by {err}");
    }

    #[test]
    fn closest_param_picks_global_minimum_on_a_wiggly_curve() {
        // 12.5: an S-shaped curve has two distance basins for a target
        // placed near the far lobe. A single-seed Newton from the
        // global-nearest coarse sample could still land in the wrong
        // basin if the seed is coarse; multi-seed safeguarded Newton
        // must find the true global closest point.
        let mut s = Sketch::new();
        // An S-curve: control points zig-zag so the curve has two
        // bulges.
        let curve = cubic_bezier(
            &mut s,
            [[0.0, 0.0], [0.0, 4.0], [4.0, -4.0], [4.0, 0.0]],
        );
        // Target near the *end* of the curve (u close to 1).
        let near_end = curve.evaluate(&s.vars, 0.92);
        // Nudge the target slightly off the curve.
        let target = [near_end[0] + 0.05, near_end[1] + 0.05];
        let u = curve.closest_param(&s.vars, target);
        // The recovered closest point must be genuinely the nearest:
        // compare against a fine brute-force scan.
        let proj = curve.evaluate(&s.vars, u);
        let got_d2 = (proj[0] - target[0]).powi(2) + (proj[1] - target[1]).powi(2);
        let mut brute_d2 = f64::INFINITY;
        for i in 0..=2000 {
            let uu = i as f64 / 2000.0;
            let p = curve.evaluate(&s.vars, uu);
            let d2 = (p[0] - target[0]).powi(2) + (p[1] - target[1]).powi(2);
            brute_d2 = brute_d2.min(d2);
        }
        assert!(
            got_d2 <= brute_d2 + 1e-6,
            "closest_param missed the global minimum: got d²={got_d2}, brute d²={brute_d2}"
        );
    }

    #[test]
    fn closest_param_clamps_to_endpoint_for_a_far_off_target() {
        // A target far past the curve's end should project to (near)
        // the endpoint parameter.
        let mut s = Sketch::new();
        let curve = cubic_bezier(
            &mut s,
            [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
        );
        // The curve is the x-axis segment [0,3]; target way past x=3.
        let u = curve.closest_param(&s.vars, [100.0, 0.0]);
        assert!((u - 1.0).abs() < 1e-6, "expected u≈1, got {u}");
    }
}
