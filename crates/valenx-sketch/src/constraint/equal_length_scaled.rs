//! EqualLengthScaled — `|a| = |b| * factor`.
//!
//! Phase 12B Task 25.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = |a| - factor * |b|.
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, factor: f64, out: &mut [f64]) {
    let (Ok(la), Ok(lb)) = (sketch.line_at(a), sketch.line_at(b)) else {
        out[0] = 0.0;
        return;
    };
    out[0] = la.length(&sketch.vars) - factor * lb.length(&sketch.vars);
}

/// Jacobian: 8 entries (4 per line, ±dx/d for each endpoint).
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    factor: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(la), Ok(lb)) = (sketch.line_at(a), sketch.line_at(b)) else {
        return;
    };
    let ((sx, sy), (ex, ey)) = la.endpoints(&sketch.vars);
    let adx = ex - sx;
    let ady = ey - sy;
    let alen = (adx * adx + ady * ady).sqrt();
    if alen > 1e-15 {
        triplets.push((0, la.start.x_var, -adx / alen));
        triplets.push((0, la.start.y_var, -ady / alen));
        triplets.push((0, la.end.x_var, adx / alen));
        triplets.push((0, la.end.y_var, ady / alen));
    }
    let ((sx, sy), (ex, ey)) = lb.endpoints(&sketch.vars);
    let bdx = ex - sx;
    let bdy = ey - sy;
    let blen = (bdx * bdx + bdy * bdy).sqrt();
    if blen > 1e-15 {
        triplets.push((0, lb.start.x_var, factor * bdx / blen));
        triplets.push((0, lb.start.y_var, factor * bdy / blen));
        triplets.push((0, lb.end.x_var, -factor * bdx / blen));
        triplets.push((0, lb.end.y_var, -factor * bdy / blen));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_a_equals_b_times_factor() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(2.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(1.0, 0.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, 2.0, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
