//! Vertical constraint — a line's endpoints have the same x.
//!
//! Residual: r = end.x - start.x.
//! Jacobian: ∂r/∂start.x = -1, ∂r/∂end.x = +1.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = end.x - start.x.
pub fn residuals(sketch: &Sketch, line: EntityId, out: &mut [f64]) {
    let l = match sketch.line_at(line) {
        Ok(l) => l,
        Err(_) => {
            out[0] = 0.0;
            return;
        }
    };
    let ((sx, _), (ex, _)) = l.endpoints(&sketch.vars);
    out[0] = ex - sx;
}

/// Jacobian: ∂r/∂start.x = -1, ∂r/∂end.x = +1.
pub fn jacobian(sketch: &Sketch, line: EntityId, triplets: &mut Vec<(usize, usize, f64)>) {
    let Ok(l) = sketch.line_at(line) else { return };
    triplets.push((0, l.start.x_var, -1.0));
    triplets.push((0, l.end.x_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_line_already_vertical() {
        let mut s = Sketch::new();
        let a = s.add_point(3.0, 0.0);
        let b = s.add_point(3.0, 5.0);
        let line = s.add_line(a, b).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, line, &mut out);
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_equals_dx_when_not_vertical() {
        let mut s = Sketch::new();
        let a = s.add_point(1.0, 0.0);
        let b = s.add_point(4.0, 5.0);
        let line = s.add_line(a, b).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, line, &mut out);
        assert!((out[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn jacobian_has_two_entries_with_correct_signs() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 1.0);
        let line = s.add_line(a, b).unwrap();
        let mut t = Vec::new();
        jacobian(&s, line, &mut t);
        assert_eq!(t.len(), 2);
        // a.x_var = 0, b.x_var = 2
        assert!(t.contains(&(0, 0, -1.0)));
        assert!(t.contains(&(0, 2, 1.0)));
    }
}
