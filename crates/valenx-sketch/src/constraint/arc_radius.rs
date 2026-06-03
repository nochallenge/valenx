//! ArcRadius — an arc has a specified radius.
//!
//! Phase 12B Task 21.
//!
//! Residual: `arc.radius - target`. Same shape as
//! [`super::radius`] when applied to a circle; this variant exists so
//! the UI can carry a stable "arc dimension" semantic.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual = radius - target.
pub fn residuals(sketch: &Sketch, arc: EntityId, target: f64, out: &mut [f64]) {
    let Ok(a) = sketch.arc_at(arc) else {
        out[0] = 0.0;
        return;
    };
    out[0] = a.radius(&sketch.vars) - target;
}

/// Jacobian: 1 entry on the radius variable.
pub fn jacobian(
    sketch: &Sketch,
    arc: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Ok(a) = sketch.arc_at(arc) else { return };
    triplets.push((0, a.radius_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_radius_matches() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let arc = s.add_arc(c, 2.5, 0.0, std::f64::consts::PI).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, arc, 2.5, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
