//! Angle constraint — two lines meet at a fixed angle (radians).
//!
//! Dot/cross form, avoiding atan2 branch cuts:
//!   r = dot(a.dir, b.dir) - |a.dir| * |b.dir| * cos(target)
//!
//! Jacobian: 8 entries.
//!   dr = dD - cos(target) * (lb * d(la) + la * d(lb))
//! where D = dot, la = |a.dir|, lb = |b.dir|.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = dot(a.dir, b.dir) - |a.dir| * |b.dir| * cos(target).
pub fn residuals(sketch: &Sketch, a: EntityId, b: EntityId, target: f64, out: &mut [f64]) {
    let (la, lb) = match (sketch.line_at(a), sketch.line_at(b)) {
        (Ok(la), Ok(lb)) => (la, lb),
        _ => {
            out[0] = 0.0;
            return;
        }
    };
    let (adx, ady) = la.direction(&sketch.vars);
    let (bdx, bdy) = lb.direction(&sketch.vars);
    let dot = adx * bdx + ady * bdy;
    let la_len = (adx * adx + ady * ady).sqrt();
    let lb_len = (bdx * bdx + bdy * bdy).sqrt();
    out[0] = dot - la_len * lb_len * target.cos();
}

/// Jacobian — 8 entries from chain rule.
pub fn jacobian(
    sketch: &Sketch,
    a: EntityId,
    b: EntityId,
    target: f64,
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
    let c = target.cos();

    // For each variable, the contribution is dD/dvar - c * (lb * d(la)/dvar + la * d(lb)/dvar).
    // Note dla/dvar = 0 for vars on line b, and vice versa.

    // a.start.x: dD = -bdx; d(la) = -adx/la; d(lb) = 0
    triplets.push((0, la.start.x_var, -bdx - c * lb_len * (-adx / la_len)));
    // a.end.x: dD = bdx; d(la) = adx/la
    triplets.push((0, la.end.x_var, bdx - c * lb_len * (adx / la_len)));
    // a.start.y: dD = -bdy; d(la) = -ady/la
    triplets.push((0, la.start.y_var, -bdy - c * lb_len * (-ady / la_len)));
    // a.end.y: dD = bdy; d(la) = ady/la
    triplets.push((0, la.end.y_var, bdy - c * lb_len * (ady / la_len)));
    // b.start.x: dD = -adx; d(lb) = -bdx/lb
    triplets.push((0, lb.start.x_var, -adx - c * la_len * (-bdx / lb_len)));
    // b.end.x: dD = adx; d(lb) = bdx/lb
    triplets.push((0, lb.end.x_var, adx - c * la_len * (bdx / lb_len)));
    // b.start.y: dD = -ady; d(lb) = -bdy/lb
    triplets.push((0, lb.start.y_var, -ady - c * la_len * (-bdy / lb_len)));
    // b.end.y: dD = ady; d(lb) = bdy/lb
    triplets.push((0, lb.end.y_var, ady - c * la_len * (bdy / lb_len)));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_angle_matches_target() {
        // a along +x, b along +y, target = pi/2.
        // dot = 0; |a|=|b|=1; cos(pi/2) = 0. r = 0 - 1*1*0 = 0.
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(1.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(0.0, 1.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, std::f64::consts::FRAC_PI_2, &mut out);
        assert!(out[0].abs() < 1e-12, "got {}", out[0]);
    }

    #[test]
    fn residual_zero_when_parallel_and_target_zero() {
        // a along +x, b along +x. target = 0. dot = 1; |a|=|b|=1; cos(0)=1.
        // r = 1 - 1*1*1 = 0.
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(1.0, 0.0);
        let b0 = s.add_point(5.0, 5.0);
        let b1 = s.add_point(6.0, 5.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, 0.0, &mut out);
        assert!(out[0].abs() < 1e-12, "got {}", out[0]);
    }

    #[test]
    fn residual_nonzero_when_angle_does_not_match() {
        // a along +x, b along +x. target = pi/2. r = 1 - 1*1*0 = 1.
        let mut s = Sketch::new();
        let a0 = s.add_point(0.0, 0.0);
        let a1 = s.add_point(1.0, 0.0);
        let b0 = s.add_point(0.0, 0.0);
        let b1 = s.add_point(1.0, 0.0);
        let la = s.add_line(a0, a1).unwrap();
        let lb = s.add_line(b0, b1).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, la, lb, std::f64::consts::FRAC_PI_2, &mut out);
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
        jacobian(&s, la, lb, std::f64::consts::FRAC_PI_2, &mut t);
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
        let target = 0.7;

        let mut out0 = vec![0.0; 1];
        residuals(&s, la, lb, target, &mut out0);
        let r0 = out0[0];

        // Perturb b1.x (var index = 7) by h.
        let h = 1e-7;
        let mut s2 = s.clone();
        s2.vars[7] += h;
        let mut out1 = vec![0.0; 1];
        residuals(&s2, la, lb, target, &mut out1);
        let numerical = (out1[0] - r0) / h;

        let mut t = Vec::new();
        jacobian(&s, la, lb, target, &mut t);
        let analytical = t
            .iter()
            .find(|(_, v, _)| *v == 7)
            .map(|(_, _, d)| *d)
            .unwrap();
        assert!(
            (numerical - analytical).abs() < 1e-4,
            "num={numerical} ana={analytical}"
        );
    }
}
