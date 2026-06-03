//! DistanceX — signed horizontal distance between two points.
//!
//! Phase 12B Task 18.
//!
//! Residual: `b.x - a.x - target`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual = b.x - a.x - target.
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, target: f64, out: &mut [f64]) {
    let (Ok(pa), Ok(pb)) = (sketch.point_at(a), sketch.point_at(b)) else {
        out[0] = 0.0;
        return;
    };
    let (ax, _) = pa.read(&sketch.vars);
    let (bx, _) = pb.read(&sketch.vars);
    out[0] = bx - ax - target;
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
    triplets.push((0, pa.x_var, -1.0));
    triplets.push((0, pb.x_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_dx_matches_target() {
        let mut s = Sketch::new();
        let a = s.add_point(1.0, 0.0);
        let b = s.add_point(4.0, 7.0);
        let mut out = vec![0.0; 1];
        residuals(&s, a, b, 3.0, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
