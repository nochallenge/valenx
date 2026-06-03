//! Move entities in place by translation.
//!
//! Phase 12E Task 42. Mutates `sketch.vars` to add `(dx, dy)` to every
//! point referenced (directly or via container entities).

use std::collections::HashSet;

use crate::geom::{Entity, EntityId};
use crate::sketch::Sketch;

/// Translate the listed entities by `delta = (dx, dy)`. Bumps the
/// underlying x/y variable values; constraints will redrive on the
/// next solve.
pub fn translate(sketch: &mut Sketch, entities: &[EntityId], delta: (f64, f64)) {
    let (dx, dy) = delta;
    // Collect unique x_var / y_var indices so a line whose endpoints
    // are shared with another entity doesn't get translated twice.
    let mut xs: HashSet<usize> = HashSet::new();
    let mut ys: HashSet<usize> = HashSet::new();
    for id in entities {
        let entity = match sketch.entities.get(id.0.wrapping_sub(1)) {
            Some(e) => e,
            None => continue,
        };
        match entity {
            Entity::Point(p) => {
                xs.insert(p.x_var);
                ys.insert(p.y_var);
            }
            Entity::Line(l) => {
                xs.insert(l.start.x_var);
                ys.insert(l.start.y_var);
                xs.insert(l.end.x_var);
                ys.insert(l.end.y_var);
            }
            Entity::Circle(c) => {
                xs.insert(c.center.x_var);
                ys.insert(c.center.y_var);
            }
            Entity::Arc(a) => {
                xs.insert(a.center.x_var);
                ys.insert(a.center.y_var);
            }
            Entity::BSpline(b) => {
                for cp in &b.control_points {
                    xs.insert(cp.x_var);
                    ys.insert(cp.y_var);
                }
            }
            Entity::Ellipse(e) => {
                xs.insert(e.center.x_var);
                ys.insert(e.center.y_var);
            }
            Entity::EllipticalArc(ea) => {
                xs.insert(ea.ellipse.center.x_var);
                ys.insert(ea.ellipse.center.y_var);
            }
        }
    }
    for x in xs {
        sketch.vars[x] += dx;
    }
    for y in ys {
        sketch.vars[y] += dy;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translate_moves_point_in_place() {
        let mut s = Sketch::new();
        let p = s.add_point(1.0, 2.0);
        translate(&mut s, &[p], (10.0, 20.0));
        let pp = s.point_at(p).unwrap();
        let (x, y) = pp.read(&s.vars);
        assert!((x - 11.0).abs() < 1e-12);
        assert!((y - 22.0).abs() < 1e-12);
    }

    #[test]
    fn translate_shared_endpoint_only_once() {
        // Two lines sharing endpoint b — translating both should not
        // double-shift the shared point.
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(2.0, 0.0);
        let ab = s.add_line(a, b).unwrap();
        let bc = s.add_line(b, c).unwrap();
        translate(&mut s, &[ab, bc], (5.0, 0.0));
        let pb = s.point_at(b).unwrap();
        let (x, _) = pb.read(&s.vars);
        assert!((x - 6.0).abs() < 1e-12, "b shifted to {x}, expected 6.0");
    }
}
