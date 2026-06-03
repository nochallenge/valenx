//! EllipseRadiusA — semi-major axis length equals a specified value.
//!
//! Phase 12B Task 23.
//!
//! Residual: `sqrt(mx^2 + my^2) - target`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = ||major_axis|| - target.
pub fn residuals(sketch: &Sketch, ellipse: EntityId, target: f64, out: &mut [f64]) {
    let Ok(e) = sketch.ellipse_at(ellipse) else {
        out[0] = 0.0;
        return;
    };
    out[0] = e.major_radius(&sketch.vars) - target;
}

/// Jacobian: 2 entries (∂/∂mx, ∂/∂my via chain rule on sqrt).
pub fn jacobian(
    sketch: &Sketch,
    ellipse: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Ok(e) = sketch.ellipse_at(ellipse) else {
        return;
    };
    let (mx, my) = e.major_axis(&sketch.vars);
    let a = (mx * mx + my * my).sqrt();
    if a < 1e-15 {
        return;
    }
    triplets.push((0, e.major_x_var, mx / a));
    triplets.push((0, e.major_y_var, my / a));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_major_radius_matches() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let e = s.add_ellipse(c, (3.0, 0.0), 1.0).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, e, 3.0, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
