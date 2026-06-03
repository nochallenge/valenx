//! Distance constraint — Euclidean distance between two points equals target.
//!
//! Residual: r = sqrt((bx-ax)² + (by-ay)²) - target.
//!
//! Jacobian (4 entries):
//!   dr/d(a.x) = -dx/d, dr/d(a.y) = -dy/d
//!   dr/d(b.x) = +dx/d, dr/d(b.y) = +dy/d
//! where d = current distance.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = |b - a| - target.
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, target: f64, out: &mut [f64]) {
    let (pa, pb) = match (sketch.point_at(a), sketch.point_at(b)) {
        (Ok(pa), Ok(pb)) => (pa, pb),
        _ => {
            out[0] = 0.0;
            return;
        }
    };
    let (ax, ay) = pa.read(&sketch.vars);
    let (bx, by) = pb.read(&sketch.vars);
    let dx = bx - ax;
    let dy = by - ay;
    out[0] = (dx * dx + dy * dy).sqrt() - target;
}

/// Jacobian — 4 entries from chain rule on sqrt.
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
    let (ax, ay) = pa.read(&sketch.vars);
    let (bx, by) = pb.read(&sketch.vars);
    let dx = bx - ax;
    let dy = by - ay;
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-15 {
        return;
    }
    triplets.push((0, pa.x_var, -dx / d));
    triplets.push((0, pa.y_var, -dy / d));
    triplets.push((0, pb.x_var, dx / d));
    triplets.push((0, pb.y_var, dy / d));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_distance_matches_target() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0); // distance 5
        let mut out = vec![0.0; 1];
        residuals(&s, a, b, 5.0, &mut out);
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_equals_actual_minus_target() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0); // distance 5
        let mut out = vec![0.0; 1];
        residuals(&s, a, b, 2.0, &mut out);
        assert!((out[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn jacobian_has_four_entries() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        let mut t = Vec::new();
        jacobian(&s, a, b, 5.0, &mut t);
        assert_eq!(t.len(), 4);
    }

    #[test]
    fn jacobian_matches_finite_difference() {
        let mut s = Sketch::new();
        let a = s.add_point(1.0, 2.0);
        let b = s.add_point(4.0, 6.0);
        let target = 3.0;

        let mut out0 = vec![0.0; 1];
        residuals(&s, a, b, target, &mut out0);
        let r0 = out0[0];

        // Perturb b.x (var index = 2) by h.
        let h = 1e-7;
        let mut s2 = s.clone();
        s2.vars[2] += h;
        let mut out1 = vec![0.0; 1];
        residuals(&s2, a, b, target, &mut out1);
        let numerical = (out1[0] - r0) / h;

        let mut t = Vec::new();
        jacobian(&s, a, b, target, &mut t);
        let analytical = t
            .iter()
            .find(|(_, v, _)| *v == 2)
            .map(|(_, _, d)| *d)
            .unwrap();
        assert!((numerical - analytical).abs() < 1e-4);
    }
}
