//! Perpendicular constraint — two lines meet at 90°.
//!
//! Residual: r = axd * bxd + ayd * byd (dot product of directions).
//!
//! Jacobian (8 entries):
//!   ∂r/∂a.start.x = -bxd, ∂r/∂a.end.x = +bxd
//!   ∂r/∂a.start.y = -byd, ∂r/∂a.end.y = +byd
//!   ∂r/∂b.start.x = -axd, ∂r/∂b.end.x = +axd
//!   ∂r/∂b.start.y = -ayd, ∂r/∂b.end.y = +ayd

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = dot(a.dir, b.dir).
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, out: &mut [f64]) {
    let (la, lb) = match (sketch.line_at(a), sketch.line_at(b)) {
        (Ok(la), Ok(lb)) => (la, lb),
        _ => {
            out[0] = 0.0;
            return;
        }
    };
    let (axd, ayd) = la.direction(&sketch.vars);
    let (bxd, byd) = lb.direction(&sketch.vars);
    out[0] = axd * bxd + ayd * byd;
}

/// Jacobian — 8 entries from chain rule on the dot product.
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(la), Ok(lb)) = (sketch.line_at(a), sketch.line_at(b)) else {
        return;
    };
    let (axd, ayd) = la.direction(&sketch.vars);
    let (bxd, byd) = lb.direction(&sketch.vars);

    triplets.push((0, la.start.x_var, -bxd));
    triplets.push((0, la.end.x_var, bxd));
    triplets.push((0, la.start.y_var, -byd));
    triplets.push((0, la.end.y_var, byd));
    triplets.push((0, lb.start.x_var, -axd));
    triplets.push((0, lb.end.x_var, axd));
    triplets.push((0, lb.start.y_var, -ayd));
    triplets.push((0, lb.end.y_var, ayd));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_lines_already_perpendicular() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(1.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(0.0, 1.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, &mut out);
        // a.dir = (1,0), b.dir = (0,1) -> dot = 0
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_nonzero_when_lines_not_perpendicular() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(2.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(3.0, 0.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, &mut out);
        // a.dir = (2,0), b.dir = (3,0) -> dot = 6
        assert!((out[0] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn jacobian_has_eight_entries() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(1.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(0.0, 1.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut t = Vec::new();
        jacobian(&s, la, lb, &mut t);
        assert_eq!(t.len(), 8);
    }

    #[test]
    fn jacobian_matches_finite_difference() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(2.0, 0.5);
        let b0 = s.add_point(1.0, 1.0);
        let b1 = s.add_point(3.0, 4.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();

        let mut out0 = vec![0.0; 1];
        residuals(&s, la, lb, &mut out0);
        let r0 = out0[0];

        // Perturb a1.x (var index = 2) by h.
        let h = 1e-7;
        let mut s2 = s.clone();
        s2.vars[2] += h;
        let mut out1 = vec![0.0; 1];
        residuals(&s2, la, lb, &mut out1);
        let numerical = (out1[0] - r0) / h;

        let mut t = Vec::new();
        jacobian(&s, la, lb, &mut t);
        let analytical = t
            .iter()
            .find(|(_, v, _)| *v == 2)
            .map(|(_, _, d)| *d)
            .unwrap();
        assert!((numerical - analytical).abs() < 1e-4);
    }
}
