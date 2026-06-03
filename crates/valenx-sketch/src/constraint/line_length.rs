//! LineLength — a line has a specified Euclidean length.
//!
//! Phase 12B Task 20. Alias of [`super::distance`] applied to the
//! line's endpoints, exposed as a first-class constraint for UI
//! ergonomics (the user selects the line, not its two endpoints).
//!
//! Residual: `sqrt((bx - ax)^2 + (by - ay)^2) - target`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = |line| - target.
pub fn residuals(sketch: &Sketch, line: EntityId, target: f64, out: &mut [f64]) {
    let Ok(l) = sketch.line_at(line) else {
        out[0] = 0.0;
        return;
    };
    let len = l.length(&sketch.vars);
    out[0] = len - target;
}

/// Jacobian: 4 entries (same form as distance).
pub fn jacobian(
    sketch: &Sketch,
    line: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Ok(l) = sketch.line_at(line) else { return };
    let ((sx, sy), (ex, ey)) = l.endpoints(&sketch.vars);
    let dx = ex - sx;
    let dy = ey - sy;
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-15 {
        return;
    }
    triplets.push((0, l.start.x_var, -dx / d));
    triplets.push((0, l.start.y_var, -dy / d));
    triplets.push((0, l.end.x_var, dx / d));
    triplets.push((0, l.end.y_var, dy / d));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_length_matches() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        let line = s.add_line(a, b).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, line, 5.0, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
