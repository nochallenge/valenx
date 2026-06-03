//! ArcAngle — an arc sweep equals a specified angle.
//!
//! Phase 12B Task 22.
//!
//! Residual: `(end_angle - start_angle) - target`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual = sweep - target.
pub fn residuals(sketch: &Sketch, arc: EntityId, target: f64, out: &mut [f64]) {
    let Ok(a) = sketch.arc_at(arc) else {
        out[0] = 0.0;
        return;
    };
    out[0] = a.sweep(&sketch.vars) - target;
}

/// Jacobian: 2 entries.
pub fn jacobian(
    sketch: &Sketch,
    arc: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Ok(a) = sketch.arc_at(arc) else { return };
    triplets.push((0, a.start_angle_var, -1.0));
    triplets.push((0, a.end_angle_var, 1.0));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_sweep_matches() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let arc = s.add_arc(c, 1.0, 0.0, std::f64::consts::PI).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, arc, std::f64::consts::PI, &mut out);
        assert!(out[0].abs() < 1e-12);
    }
}
