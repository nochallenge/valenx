//! EllipseRadiusB — semi-minor magnitude equals a specified value.
//!
//! Phase 12B Task 24.
//!
//! Residual: `minor_radius - target`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual = minor_radius - target.
pub fn residuals(sketch: &Sketch, ellipse: EntityId, target: f64, out: &mut [f64]) {
    let Ok(e) = sketch.ellipse_at(ellipse) else {
        out[0] = 0.0;
        return;
    };
    out[0] = e.minor_radius(&sketch.vars) - target;
}

/// Jacobian: 1 entry.
pub fn jacobian(
    sketch: &Sketch,
    ellipse: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Ok(e) = sketch.ellipse_at(ellipse) else {
        return;
    };
    triplets.push((0, e.minor_radius_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_minor_radius_matches() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let e = s.add_ellipse(c, (3.0, 0.0), 1.5).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, e, 1.5, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
