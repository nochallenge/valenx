//! Point-on-ellipse constraint — a point lies on an ellipse's
//! perimeter.
//!
//! Phase 12B Task 16.
//!
//! Residual: implicit ellipse equation evaluated at the point. For
//! an axis-aligned ellipse centred at the origin with semi-axes (a,
//! b): `(x/a)^2 + (y/b)^2 - 1 = 0`. For a general orientation we
//! rotate the point into the ellipse's local frame first.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual: `(p_local.x / a)^2 + (p_local.y / b)^2 - 1` where
/// `p_local` is the point in the ellipse's local frame.
pub fn residuals(sketch: &Sketch, point: EntityId, ellipse: EntityId, out: &mut [f64]) {
    let (Ok(p), Ok(e)) = (sketch.point_at(point), sketch.ellipse_at(ellipse)) else {
        out[0] = 0.0;
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let (cx, cy) = e.center_xy(&sketch.vars);
    let (mx, my) = e.major_axis(&sketch.vars);
    let a = (mx * mx + my * my).sqrt();
    let b = e.minor_radius(&sketch.vars);
    if a < 1e-15 || b < 1e-15 {
        out[0] = 0.0;
        return;
    }
    let dx = px - cx;
    let dy = py - cy;
    // Local x = dot(d, major_unit), local y = dot(d, perp_major_unit).
    let ux = mx / a;
    let uy = my / a;
    let vx = -uy;
    let vy = ux;
    let local_x = dx * ux + dy * uy;
    let local_y = dx * vx + dy * vy;
    out[0] = (local_x / a).powi(2) + (local_y / b).powi(2) - 1.0;
}

/// Jacobian: finite-difference in v1 — gives 6 entries (point.x/y plus
/// ellipse center.x/y plus major.x/y plus minor_radius). Numerical
/// approximation keeps the constraint cheap; analytic form is a v1.5
/// task once the constraint sees heavy use.
pub fn jacobian(
    sketch: &Sketch,
    point: EntityId,
    ellipse: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(p), Ok(e)) = (sketch.point_at(point), sketch.ellipse_at(ellipse)) else {
        return;
    };
    let vars = vec![
        p.x_var,
        p.y_var,
        e.center.x_var,
        e.center.y_var,
        e.major_x_var,
        e.major_y_var,
        e.minor_radius_var,
    ];
    let mut base = vec![0.0; 1];
    residuals(sketch, point, ellipse, &mut base);
    let r0 = base[0];
    let h = 1e-7;
    let mut perturbed = sketch.clone();
    for var in vars {
        let saved = perturbed.vars[var];
        perturbed.vars[var] = saved + h;
        let mut r = vec![0.0; 1];
        residuals(&perturbed, point, ellipse, &mut r);
        let d = (r[0] - r0) / h;
        if d.abs() > 1e-15 {
            triplets.push((0, var, d));
        }
        perturbed.vars[var] = saved;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_point_on_axis_aligned_ellipse() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let e = s.add_ellipse(c, (2.0, 0.0), 1.0).unwrap();
        let p = s.add_point(2.0, 0.0);
        let mut out = vec![0.0; 1];
        residuals(&s, p, e, &mut out);
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_zero_at_minor_axis_point() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let e = s.add_ellipse(c, (2.0, 0.0), 1.0).unwrap();
        let p = s.add_point(0.0, 1.0);
        let mut out = vec![0.0; 1];
        residuals(&s, p, e, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
