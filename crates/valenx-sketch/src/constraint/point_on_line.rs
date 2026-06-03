//! Point-on-line constraint — a point lies on the infinite line
//! through a line segment.
//!
//! Phase 12B Task 12.
//!
//! Residual: signed perpendicular distance from `p` to the line
//! through `(s.start, s.end)`. Equivalent to the 2-D cross product
//! `(p - start) × (end - start)`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = (px - sx)*(ey - sy) - (py - sy)*(ex - sx).
pub fn residuals(sketch: &Sketch, point: EntityId, line: EntityId, out: &mut [f64]) {
    let (Ok(p), Ok(l)) = (sketch.point_at(point), sketch.line_at(line)) else {
        out[0] = 0.0;
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let ((sx, sy), (ex, ey)) = l.endpoints(&sketch.vars);
    out[0] = (px - sx) * (ey - sy) - (py - sy) * (ex - sx);
}

/// Jacobian: 6 entries from derivatives of the cross product.
pub fn jacobian(
    sketch: &Sketch,
    point: EntityId,
    line: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(p), Ok(l)) = (sketch.point_at(point), sketch.line_at(line)) else {
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let ((sx, sy), (ex, ey)) = l.endpoints(&sketch.vars);
    let dx = ex - sx;
    let dy = ey - sy;
    // d r / d px = dy, d r / d py = -dx
    triplets.push((0, p.x_var, dy));
    triplets.push((0, p.y_var, -dx));
    // d r / d sx = -dy + (py - sy)
    // r = (px - sx)*(ey - sy) - (py - sy)*(ex - sx)
    // d/dsx = -(ey - sy) - (py - sy)*(-1) = -dy + (py - sy)
    triplets.push((0, l.start.x_var, -dy + (py - sy)));
    // d/dsy = (px - sx)*(-1) - (-1)*(ex - sx) = -(px - sx) + dx
    triplets.push((0, l.start.y_var, dx - (px - sx)));
    // d/dex = -(py - sy)
    triplets.push((0, l.end.x_var, -(py - sy)));
    // d/dey = (px - sx)
    triplets.push((0, l.end.y_var, px - sx));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_point_on_line() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(4.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        let p = s.add_point(2.0, 0.0);
        let mut out = vec![0.0; 1];
        residuals(&s, p, line, &mut out);
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_positive_when_point_above_line() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(4.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        let p = s.add_point(2.0, 1.0);
        let mut out = vec![0.0; 1];
        residuals(&s, p, line, &mut out);
        // (2 - 0)*(0 - 0) - (1 - 0)*(4 - 0) = -4
        assert!((out[0] + 4.0).abs() < 1e-12);
    }
}
