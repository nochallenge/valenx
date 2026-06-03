//! EqualLength constraint — two line segments have equal length.
//!
//! Residual: r = |a| - |b| where |line| = sqrt(dx² + dy²).
//!
//! Jacobian (8 entries, chain rule):
//!   dr/d(a.start.x) = -adx/la, dr/d(a.end.x) = +adx/la
//!   dr/d(a.start.y) = -ady/la, dr/d(a.end.y) = +ady/la
//!   dr/d(b.start.x) = +bdx/lb, dr/d(b.end.x) = -bdx/lb
//!   dr/d(b.start.y) = +bdy/lb, dr/d(b.end.y) = -bdy/lb

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = a.length - b.length.
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, out: &mut [f64]) {
    let (la, lb) = match (sketch.line_at(a), sketch.line_at(b)) {
        (Ok(la), Ok(lb)) => (la, lb),
        _ => {
            out[0] = 0.0;
            return;
        }
    };
    out[0] = la.length(&sketch.vars) - lb.length(&sketch.vars);
}

/// Jacobian — 8 entries from chain rule on sqrt(dx² + dy²).
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(la), Ok(lb)) = (sketch.line_at(a), sketch.line_at(b)) else {
        return;
    };
    let (adx, ady) = la.direction(&sketch.vars);
    let (bdx, bdy) = lb.direction(&sketch.vars);
    let la_len = (adx * adx + ady * ady).sqrt();
    let lb_len = (bdx * bdx + bdy * bdy).sqrt();
    if la_len < 1e-15 || lb_len < 1e-15 {
        return;
    }

    triplets.push((0, la.start.x_var, -adx / la_len));
    triplets.push((0, la.end.x_var, adx / la_len));
    triplets.push((0, la.start.y_var, -ady / la_len));
    triplets.push((0, la.end.y_var, ady / la_len));
    triplets.push((0, lb.start.x_var, bdx / lb_len));
    triplets.push((0, lb.end.x_var, -bdx / lb_len));
    triplets.push((0, lb.start.y_var, bdy / lb_len));
    triplets.push((0, lb.end.y_var, -bdy / lb_len));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_lengths_already_equal() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(3.0, 4.0); // length 5
        let b0 = s.add_point(10.0, 10.0);
        let b1 = s.add_point(13.0, 14.0); // length 5
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, &mut out);
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_equals_length_difference() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(3.0, 4.0); // length 5
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(1.0, 0.0); // length 1
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, &mut out);
        assert!((out[0] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn jacobian_has_eight_entries() {
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(3.0, 4.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(1.0, 0.0);
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
        let a1 = s.add_point(3.0, 4.0);
        let b0 = s.add_point(1.0, 1.0);
        let b1 = s.add_point(4.0, 2.0);
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
