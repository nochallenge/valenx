//! Point-on-arc constraint — a point lies on an arc's circular path.
//!
//! Phase 12B Task 14.
//!
//! Residual: same as point-on-circle (distance(point, center) -
//! radius). Angular bounds are not enforced in v1 — the solver
//! drives the point onto the underlying circle; the user is expected
//! to add a separate angular constraint if needed.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = sqrt((px - cx)^2 + (py - cy)^2) - radius.
pub fn residuals(sketch: &Sketch, point: EntityId, arc: EntityId, out: &mut [f64]) {
    let (Ok(p), Ok(a)) = (sketch.point_at(point), sketch.arc_at(arc)) else {
        out[0] = 0.0;
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let (cx, cy) = a.center.read(&sketch.vars);
    let r = a.radius(&sketch.vars);
    let dx = px - cx;
    let dy = py - cy;
    out[0] = (dx * dx + dy * dy).sqrt() - r;
}

/// Jacobian: 5 entries.
pub fn jacobian(
    sketch: &Sketch,
    point: EntityId,
    arc: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(p), Ok(a)) = (sketch.point_at(point), sketch.arc_at(arc)) else {
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let (cx, cy) = a.center.read(&sketch.vars);
    let dx = px - cx;
    let dy = py - cy;
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-15 {
        return;
    }
    triplets.push((0, p.x_var, dx / d));
    triplets.push((0, p.y_var, dy / d));
    triplets.push((0, a.center.x_var, -dx / d));
    triplets.push((0, a.center.y_var, -dy / d));
    triplets.push((0, a.radius_var, -1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_point_on_arc_path() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let arc = s.add_arc(c, 2.0, 0.0, std::f64::consts::PI).unwrap();
        let p = s.add_point(2.0, 0.0);
        let mut out = vec![0.0; 1];
        residuals(&s, p, arc, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
