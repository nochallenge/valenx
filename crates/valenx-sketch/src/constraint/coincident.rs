//! Coincident constraint: two points share xy.
//!
//! Residuals: (b.x - a.x, b.y - a.y). Solver drives both to zero.
//! Jacobian: ∂r₀/∂a.x = -1, ∂r₀/∂b.x = +1, similar for y.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residuals = (b.x - a.x, b.y - a.y).
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, out: &mut [f64]) {
    let pa = match sketch.point_at(a) {
        Ok(p) => p,
        Err(_) => {
            out[0] = 0.0;
            out[1] = 0.0;
            return;
        }
    };
    let pb = match sketch.point_at(b) {
        Ok(p) => p,
        Err(_) => {
            out[0] = 0.0;
            out[1] = 0.0;
            return;
        }
    };
    let (ax, ay) = pa.read(&sketch.vars);
    let (bx, by) = pb.read(&sketch.vars);
    out[0] = bx - ax;
    out[1] = by - ay;
}

/// Jacobian: per-coordinate ±1 entries on the involved variables.
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(pa), Ok(pb)) = (sketch.point_at(a), sketch.point_at(b)) else {
        return;
    };
    triplets.push((0, pa.x_var, -1.0));
    triplets.push((0, pb.x_var, 1.0));
    triplets.push((1, pa.y_var, -1.0));
    triplets.push((1, pb.y_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residuals_zero_when_points_already_coincident() {
        let mut s = Sketch::new();
        let a = s.add_point(1.0, 2.0);
        let b = s.add_point(1.0, 2.0);
        let mut out = vec![0.0; 2];
        residuals(&s, a, b, &mut out);
        assert_eq!(out, vec![0.0, 0.0]);
    }

    #[test]
    fn residuals_equal_delta_when_separated() {
        let mut s = Sketch::new();
        let a = s.add_point(1.0, 2.0);
        let b = s.add_point(4.0, 7.0);
        let mut out = vec![0.0; 2];
        residuals(&s, a, b, &mut out);
        assert_eq!(out, vec![3.0, 5.0]);
    }

    #[test]
    fn jacobian_has_four_entries_with_correct_signs() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 1.0);
        let mut t = Vec::new();
        jacobian(&s, a, b, &mut t);
        // 4 entries: ∂r0/∂a.x = -1, ∂r0/∂b.x = +1, ∂r1/∂a.y = -1, ∂r1/∂b.y = +1
        assert_eq!(t.len(), 4);
        assert!(t.contains(&(0, 0, -1.0))); // a.x var = 0
        assert!(t.contains(&(0, 2, 1.0))); // b.x var = 2
        assert!(t.contains(&(1, 1, -1.0))); // a.y var = 1
        assert!(t.contains(&(1, 3, 1.0))); // b.y var = 3
    }
}
