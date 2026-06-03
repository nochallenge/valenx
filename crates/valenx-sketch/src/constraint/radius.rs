//! Radius constraint — a circle's or arc's radius equals a target.
//!
//! Residual: r = radius - target.
//! Jacobian: dr/d(radius_var) = +1 (single entry).

use crate::geom::{Entity, EntityId};
use crate::sketch::Sketch;

/// Residual r = radius - target.
pub fn residuals(sketch: &Sketch, circle_or_arc: EntityId, target: f64, out: &mut [f64]) {
    let radius = match radius_of(sketch, circle_or_arc) {
        Some(r) => r,
        None => {
            out[0] = 0.0;
            return;
        }
    };
    out[0] = radius - target;
}

/// Jacobian: dr/d(radius_var) = +1.
pub fn jacobian(
    sketch: &Sketch,
    circle_or_arc: EntityId,
    _target: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let Some(rv) = radius_var_of(sketch, circle_or_arc) else {
        return;
    };
    triplets.push((0, rv, 1.0));
}

/// Read the current radius value from either a Circle or an Arc.
fn radius_of(sketch: &Sketch, id: EntityId) -> Option<f64> {
    let entity = sketch.entities.get(id.0.wrapping_sub(1))?;
    match entity {
        Entity::Circle(c) => Some(c.radius(&sketch.vars)),
        Entity::Arc(a) => Some(a.radius(&sketch.vars)),
        _ => None,
    }
}

/// Get the radius variable index from either a Circle or an Arc.
fn radius_var_of(sketch: &Sketch, id: EntityId) -> Option<usize> {
    let entity = sketch.entities.get(id.0.wrapping_sub(1))?;
    match entity {
        Entity::Circle(c) => Some(c.radius_var),
        Entity::Arc(a) => Some(a.radius_var),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_when_radius_matches_target() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let circle = s.add_circle(c, 2.5).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, circle, 2.5, &mut out);
        assert!(out[0].abs() < 1e-12);
    }

    #[test]
    fn residual_equals_actual_minus_target() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let circle = s.add_circle(c, 5.0).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, circle, 2.0, &mut out);
        assert!((out[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn jacobian_has_single_entry_with_value_one() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let circle = s.add_circle(c, 2.5).unwrap();
        let mut t = Vec::new();
        jacobian(&s, circle, 2.5, &mut t);
        assert_eq!(t.len(), 1);
        // c.x = 0, c.y = 1, radius_var = 2
        assert_eq!(t[0], (0, 2, 1.0));
    }

    #[test]
    fn works_on_arc_radius_too() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let arc = s.add_arc(c, 4.0, 0.0, std::f64::consts::PI).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, arc, 1.5, &mut out);
        assert!((out[0] - 2.5).abs() < 1e-12);

        let mut t = Vec::new();
        jacobian(&s, arc, 1.5, &mut t);
        assert_eq!(t.len(), 1);
        // c.x=0, c.y=1, radius_var=2, start_angle=3, end_angle=4
        assert_eq!(t[0], (0, 2, 1.0));
    }
}
