//! Point-on-circle constraint — a point lies on a circle's perimeter.
//!
//! Phase 12B Task 13.
//!
//! Residual: distance(point, center) - radius.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = sqrt((px - cx)^2 + (py - cy)^2) - radius.
pub fn residuals(sketch: &Sketch, point: EntityId, circle: EntityId, out: &mut [f64]) {
    let (Ok(p), Ok(c)) = (sketch.point_at(point), sketch.circle_at(circle)) else {
        out[0] = 0.0;
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let (cx, cy) = c.center.read(&sketch.vars);
    let r = c.radius(&sketch.vars);
    let dx = px - cx;
    let dy = py - cy;
    out[0] = (dx * dx + dy * dy).sqrt() - r;
}

/// Jacobian: 5 entries — d r/d {px, py, cx, cy, radius}.
pub fn jacobian(
    sketch: &Sketch,
    point: EntityId,
    circle: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(p), Ok(c)) = (sketch.point_at(point), sketch.circle_at(circle)) else {
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let (cx, cy) = c.center.read(&sketch.vars);
    let dx = px - cx;
    let dy = py - cy;
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-15 {
        return;
    }
    triplets.push((0, p.x_var, dx / d));
    triplets.push((0, p.y_var, dy / d));
    triplets.push((0, c.center.x_var, -dx / d));
    triplets.push((0, c.center.y_var, -dy / d));
    triplets.push((0, c.radius_var, -1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_point_on_circle() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let circle = s.add_circle(c, 5.0).unwrap();
        let p = s.add_point(3.0, 4.0);
        let mut out = vec![0.0; 1];
        residuals(&s, p, circle, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
