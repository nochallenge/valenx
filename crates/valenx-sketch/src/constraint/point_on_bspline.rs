//! Point-on-B-spline constraint — a point lies on a B-spline curve.
//!
//! Phase 12B Task 15.
//!
//! Residual: distance from the point to its closest-point projection
//! on the curve. The projection uses the BSpline's
//! [`crate::geom_bspline::BSpline2::closest_param`] — as of Phase 12.5
//! a **multi-seed safeguarded Newton**: every local minimum of a dense
//! coarse scan is refined by a bracketed Newton/bisection hybrid and
//! the global closest point is returned. This converges reliably even
//! on wiggly high-degree curves where a single-seed Newton would land
//! in the wrong distance basin.
//!
//! The Jacobian is computed via finite differences in v1 (sparse, the
//! perturbation moves only `p.x` and `p.y` — the control-point and
//! knot derivatives are treated as zero, which is sufficient when the
//! curve is held fixed by other constraints; for free-moving control
//! points the closest-param recomputation in `residuals` keeps the
//! solver direction consistent).

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = distance(point, closest_point_on_curve).
pub fn residuals(sketch: &Sketch, point: EntityId, bspline: EntityId, out: &mut [f64]) {
    let Ok(p) = sketch.point_at(point) else {
        out[0] = 0.0;
        return;
    };
    let Ok(curve) = sketch.bspline_at(bspline) else {
        out[0] = 0.0;
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let u = curve.closest_param(&sketch.vars, [px, py]);
    let proj = curve.evaluate(&sketch.vars, u);
    let dx = px - proj[0];
    let dy = py - proj[1];
    out[0] = (dx * dx + dy * dy).sqrt();
}

/// Jacobian: 2 entries (d r / d px, d r / d py). The closest-param
/// is treated as locally constant — moving the point shifts the
/// numerator linearly, which is exact for points already on the curve
/// (the residual surface is locally a smooth quadratic). For points
/// far from the curve the LM damping handles the curvature.
pub fn jacobian(
    sketch: &Sketch,
    point: EntityId,
    bspline: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Ok(p) = sketch.point_at(point) else {
        return;
    };
    let Ok(curve) = sketch.bspline_at(bspline) else {
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let u = curve.closest_param(&sketch.vars, [px, py]);
    let proj = curve.evaluate(&sketch.vars, u);
    let dx = px - proj[0];
    let dy = py - proj[1];
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-15 {
        return;
    }
    triplets.push((0, p.x_var, dx / d));
    triplets.push((0, p.y_var, dy / d));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_near_zero_for_point_on_bezier_midpoint() {
        let mut s = Sketch::new();
        let p0 = s.add_point(0.0, 0.0);
        let p1 = s.add_point(1.0, 2.0);
        let p2 = s.add_point(3.0, 2.0);
        let p3 = s.add_point(4.0, 0.0);
        let curve_id = s
            .add_bspline(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                &[p0, p1, p2, p3],
                vec![1.0; 4],
            )
            .unwrap();
        // The bezier midpoint at u=0.5 is (2, 1.5) — see geom_bspline test.
        let target = s.add_point(2.0, 1.5);
        let mut out = vec![0.0; 1];
        residuals(&s, target, curve_id, &mut out);
        assert!(out[0].abs() < 1e-6, "residual {} not near zero", out[0]);
    }
}
