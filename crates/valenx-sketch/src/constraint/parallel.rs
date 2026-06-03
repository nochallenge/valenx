//! Parallel constraint — two lines have parallel direction vectors.
//!
//! Residual: r = axd * byd - ayd * bxd (2-D cross product of directions).
//! Where axd = a.end.x - a.start.x, ayd = a.end.y - a.start.y, etc.
//!
//! Jacobian (8 entries, chain rule):
//!   ∂r/∂a.start.x = -byd, ∂r/∂a.end.x = +byd
//!   ∂r/∂a.start.y = +bxd, ∂r/∂a.end.y = -bxd
//!   ∂r/∂b.start.x = +ayd, ∂r/∂b.end.x = -ayd
//!   ∂r/∂b.start.y = -axd, ∂r/∂b.end.y = +axd

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = cross(a.dir, b.dir).
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
    out[0] = axd * byd - ayd * bxd;
}

/// Jacobian — 8 entries from chain rule.
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

    // d r / d a.start.x = -byd, d r / d a.end.x = +byd
    triplets.push((0, la.start.x_var, -byd));
    triplets.push((0, la.end.x_var, byd));
    // d r / d a.start.y = +bxd, d r / d a.end.y = -bxd
    triplets.push((0, la.start.y_var, bxd));
    triplets.push((0, la.end.y_var, -bxd));
    // d r / d b.start.x = +ayd, d r / d b.end.x = -ayd
    triplets.push((0, lb.start.x_var, ayd));
    triplets.push((0, lb.end.x_var, -ayd));
    // d r / d b.start.y = -axd, d r / d b.end.y = +axd
    triplets.push((0, lb.start.y_var, -axd));
    triplets.push((0, lb.end.y_var, axd));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_lines_already_parallel() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(2.0, 1.0);
        let b0 = s.add_point(5.0, 5.0);
        let b1 = s.add_point(9.0, 7.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, &mut out);
        // a.dir = (2,1), b.dir = (4,2) -> cross = 2*2 - 1*4 = 0
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_nonzero_when_lines_not_parallel() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(1.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(0.0, 1.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, &mut out);
        // a.dir = (1,0), b.dir = (0,1) -> cross = 1*1 - 0*0 = 1
        assert!((out[0] - 1.0).abs() < 1e-12);
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
        // a = horizontal (1,0), b = (1,1) -> dir = (-1,1) at b1, sketched as two points
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

        // Perturb b1.x (var index = 7) by h.
        let h = 1e-7;
        let mut s2 = s.clone();
        s2.vars[7] += h;
        let mut out1 = vec![0.0; 1];
        residuals(&s2, la, lb, &mut out1);
        let numerical = (out1[0] - r0) / h;

        let mut t = Vec::new();
        jacobian(&s, la, lb, &mut t);
        let analytical = t
            .iter()
            .find(|(_, v, _)| *v == 7)
            .map(|(_, _, d)| *d)
            .unwrap();
        assert!((numerical - analytical).abs() < 1e-4);
    }
}
