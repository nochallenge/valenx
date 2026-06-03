//! Symmetric constraint — two points are symmetric about a third
//! point (the midpoint), i.e. `c = (a + b) / 2`.
//!
//! Phase 12B Task 10.
//!
//! Residuals: `(a.x + b.x - 2*c.x, a.y + b.y - 2*c.y)`. Drives the
//! midpoint of (a, b) to coincide with `c`. The v1 form takes a
//! symmetry *point*; the FreeCAD-style line symmetry can be expressed
//! by additionally constraining `c` to lie on the line via
//! [`super::point_on_line`].

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residuals = (a.x + b.x - 2*c.x, a.y + b.y - 2*c.y).
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, c: EntityId, out: &mut [f64]) {
    let (pa, pb, pc) = match (sketch.point_at(a), sketch.point_at(b), sketch.point_at(c)) {
        (Ok(pa), Ok(pb), Ok(pc)) => (pa, pb, pc),
        _ => {
            out[0] = 0.0;
            out[1] = 0.0;
            return;
        }
    };
    let (ax, ay) = pa.read(&sketch.vars);
    let (bx, by) = pb.read(&sketch.vars);
    let (cx, cy) = pc.read(&sketch.vars);
    out[0] = ax + bx - 2.0 * cx;
    out[1] = ay + by - 2.0 * cy;
}

/// Jacobian: 6 entries per row.
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    c: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(pa), Ok(pb), Ok(pc)) = (sketch.point_at(a), sketch.point_at(b), sketch.point_at(c))
    else {
        return;
    };
    triplets.push((0, pa.x_var, 1.0));
    triplets.push((0, pb.x_var, 1.0));
    triplets.push((0, pc.x_var, -2.0));
    triplets.push((1, pa.y_var, 1.0));
    triplets.push((1, pb.y_var, 1.0));
    triplets.push((1, pc.y_var, -2.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_c_is_midpoint() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(4.0, 6.0);
        let c = s.add_point(2.0, 3.0);
        let mut out = vec![0.0; 2];
        residuals(&s, a, b, c, &mut out);
        assert!(out[0].abs() < 1e-12 && out[1].abs() < 1e-12);
    }
}
