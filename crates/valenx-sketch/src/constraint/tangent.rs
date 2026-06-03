//! Tangent constraint — line tangent to a circle (line-circle case).
//!
//! v1 supports only line-circle. For circle-circle the constraint is a
//! v1 gap: it writes zero residuals and zero Jacobian (no-op).
//!
//! Math:
//!   signed_dist = ((ex-sx) * (cy-sy) - (ey-sy) * (cx-sx)) / |line|
//!   r = |signed_dist| - radius
//!
//! Jacobian uses sign(signed_dist) and the quotient rule on the
//! signed-distance fraction.

use crate::geom::{Entity, EntityId};
use crate::sketch::Sketch;

/// Residual r = distance(centre, line) - radius for line-circle.
/// For circle-circle (unsupported in v1), returns zero.
pub fn residuals(sketch: &Sketch, line_or_circle_a: EntityId, circle_b: EntityId, out: &mut [f64]) {
    // Resolve line + circle (only line-circle is supported in v1).
    let (line, circle) = match resolve_line_circle(sketch, line_or_circle_a, circle_b) {
        Some(pair) => pair,
        None => {
            out[0] = 0.0;
            return;
        }
    };
    let ((sx, sy), (ex, ey)) = line.endpoints(&sketch.vars);
    let (cx, cy) = circle.center.read(&sketch.vars);
    let radius = circle.radius(&sketch.vars);
    let dx = ex - sx;
    let dy = ey - sy;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-15 {
        out[0] = 0.0;
        return;
    }
    let signed = (dx * (cy - sy) - dy * (cx - sx)) / len;
    out[0] = signed.abs() - radius;
}

/// Jacobian for line-circle tangent. For circle-circle (unsupported), no-op.
pub fn jacobian(
    sketch: &Sketch,
    line_or_circle_a: EntityId,
    circle_b: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (line, circle) = match resolve_line_circle(sketch, line_or_circle_a, circle_b) {
        Some(pair) => pair,
        None => return,
    };
    let ((sx, sy), (ex, ey)) = line.endpoints(&sketch.vars);
    let (cx, cy) = circle.center.read(&sketch.vars);
    let dx = ex - sx;
    let dy = ey - sy;
    let l_sq = dx * dx + dy * dy;
    let l = l_sq.sqrt();
    if l < 1e-15 {
        return;
    }
    let n = dx * (cy - sy) - dy * (cx - sx);
    let signed = n / l;
    let sign = if signed >= 0.0 { 1.0 } else { -1.0 };

    // ds/dvar = (dn/dvar * L - n * dL/dvar) / L²
    // dr/dvar = sign(signed) * ds/dvar (for endpoint+center vars)
    // dr/dradius_var = -1
    let push = |triplets: &mut Vec<(usize, usize, f64)>, var: usize, dn: f64, dl: f64| {
        let dsdv = (dn * l - n * dl) / l_sq;
        triplets.push((0, var, sign * dsdv));
    };

    // start.x: dn = dy - (cy - sy);  dl = -dx / l
    push(triplets, line.start.x_var, dy - (cy - sy), -dx / l);
    // start.y: dn = (cx - sx) - dx;  dl = -dy / l
    push(triplets, line.start.y_var, (cx - sx) - dx, -dy / l);
    // end.x: dn = cy - sy; dl = dx / l
    push(triplets, line.end.x_var, cy - sy, dx / l);
    // end.y: dn = -(cx - sx); dl = dy / l
    push(triplets, line.end.y_var, -(cx - sx), dy / l);
    // circle.center.x: dn = -dy; dl = 0
    push(triplets, circle.center.x_var, -dy, 0.0);
    // circle.center.y: dn = dx; dl = 0
    push(triplets, circle.center.y_var, dx, 0.0);
    // radius: dr/dradius_var = -1
    triplets.push((0, circle.radius_var, -1.0));
}

/// Resolve the (line, circle) pair from the two entity ids. Either may
/// be the line. Returns None if neither is a line, or if neither is a
/// circle (e.g. circle-circle, which is a v1 gap).
fn resolve_line_circle(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
) -> Option<(crate::geom::Line2, crate::geom::Circle2)> {
    let ea = sketch.entities.get(a.0.wrapping_sub(1))?;
    let eb = sketch.entities.get(b.0.wrapping_sub(1))?;
    match (ea, eb) {
        (Entity::Line(l), Entity::Circle(c)) => Some((*l, *c)),
        (Entity::Circle(c), Entity::Line(l)) => Some((*l, *c)),
        _ => None, // circle-circle and other combos: v1 gap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_line_tangent_to_circle() {
        // Line y = 0 (along x-axis), circle centred at (0, 1) radius 1.
        // Distance from (0,1) to y=0 line is 1 = radius. Tangent.
        let mut s = Sketch::new();
        let a = s.add_point(-1.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        let c = s.add_point(0.0, 1.0);
        let circle = s.add_circle(c, 1.0).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, line, circle, &mut out);
        assert!(out[0].abs() < 1e-12, "got {}", out[0]);
    }

    #[test]
    fn residual_nonzero_when_not_tangent() {
        // Line y=0; circle centre (0,2) radius 1. Distance = 2, residual = 1.
        let mut s = Sketch::new();
        let a = s.add_point(-1.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        let c = s.add_point(0.0, 2.0);
        let circle = s.add_circle(c, 1.0).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, line, circle, &mut out);
        assert!((out[0] - 1.0).abs() < 1e-12, "got {}", out[0]);
    }

    #[test]
    fn jacobian_has_seven_entries() {
        // 4 line endpoint coords + 2 centre coords + 1 radius = 7.
        let mut s = Sketch::new();
        let a = s.add_point(-1.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        let c = s.add_point(0.0, 2.0);
        let circle = s.add_circle(c, 1.0).unwrap();
        let mut t = Vec::new();
        jacobian(&s, line, circle, &mut t);
        assert_eq!(t.len(), 7);
    }

    #[test]
    fn jacobian_matches_finite_difference() {
        // Same setup as nonzero test. Perturb cy.
        let mut s = Sketch::new();
        let a = s.add_point(-1.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        let c = s.add_point(0.0, 2.0);
        let circle = s.add_circle(c, 1.0).unwrap();

        let mut out0 = vec![0.0; 1];
        residuals(&s, line, circle, &mut out0);
        let r0 = out0[0];

        // cy is var index 5 (a=0,1; b=2,3; c=4,5; radius=6).
        let h = 1e-7;
        let mut s2 = s.clone();
        s2.vars[5] += h;
        let mut out1 = vec![0.0; 1];
        residuals(&s2, line, circle, &mut out1);
        let numerical = (out1[0] - r0) / h;

        let mut t = Vec::new();
        jacobian(&s, line, circle, &mut t);
        let analytical = t
            .iter()
            .find(|(_, v, _)| *v == 5)
            .map(|(_, _, d)| *d)
            .unwrap();
        assert!(
            (numerical - analytical).abs() < 1e-4,
            "num={numerical} ana={analytical}"
        );
    }

    #[test]
    fn circle_circle_is_noop_v1_gap() {
        let mut s = Sketch::new();
        let c1 = s.add_point(0.0, 0.0);
        let circle1 = s.add_circle(c1, 1.0).unwrap();
        let c2 = s.add_point(5.0, 0.0);
        let circle2 = s.add_circle(c2, 1.0).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, circle1, circle2, &mut out);
        assert_eq!(out[0], 0.0);
        let mut t = Vec::new();
        jacobian(&s, circle1, circle2, &mut t);
        assert!(t.is_empty());
    }
}
