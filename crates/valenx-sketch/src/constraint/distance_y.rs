//! DistanceY — signed vertical distance between two points.
//!
//! Phase 12B Task 19.
//!
//! Residual: `b.y - a.y - target`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual = b.y - a.y - target.
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, target: f64, out: &mut [f64]) {
    let (Ok(pa), Ok(pb)) = (sketch.point_at(a), sketch.point_at(b)) else {
        out[0] = 0.0;
        return;
    };
    let (_, ay) = pa.read(&sketch.vars);
    let (_, by) = pb.read(&sketch.vars);
    out[0] = by - ay - target;
}

/// Jacobian: 2 entries.
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(pa), Ok(pb)) = (sketch.point_at(a), sketch.point_at(b)) else {
        return;
    };
    triplets.push((0, pa.y_var, -1.0));
    triplets.push((0, pb.y_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_dy_matches_target() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 1.0);
        let b = s.add_point(7.0, 4.0);
        let mut out = vec![0.0; 1];
        residuals(&s, a, b, 3.0, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
