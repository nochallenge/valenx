//! Rotate entities about a pivot.
//!
//! Phase 12E Task 43.

use std::collections::HashSet;

use crate::geom::{Entity, EntityId};
use crate::sketch::Sketch;

/// Rotate every point referenced by `entities` about `pivot` by
/// `angle` radians (CCW).
pub fn rotate(sketch: &mut Sketch, entities: &[EntityId], pivot: (f64, f64), angle: f64) {
    let (px, py) = pivot;
    let c = angle.cos();
    let s = angle.sin();
    let mut points: HashSet<(usize, usize)> = HashSet::new();
    for id in entities {
        let entity = match sketch.entities.get(id.0.wrapping_sub(1)) {
            Some(e) => e,
            None => continue,
        };
        match entity {
            Entity::Point(p) => {
                points.insert((p.x_var, p.y_var));
            }
            Entity::Line(l) => {
                points.insert((l.start.x_var, l.start.y_var));
                points.insert((l.end.x_var, l.end.y_var));
            }
            Entity::Circle(circle) => {
                points.insert((circle.center.x_var, circle.center.y_var));
            }
            Entity::Arc(a) => {
                points.insert((a.center.x_var, a.center.y_var));
            }
            Entity::BSpline(b) => {
                for cp in &b.control_points {
                    points.insert((cp.x_var, cp.y_var));
                }
            }
            Entity::Ellipse(e) => {
                points.insert((e.center.x_var, e.center.y_var));
            }
            Entity::EllipticalArc(ea) => {
                points.insert((ea.ellipse.center.x_var, ea.ellipse.center.y_var));
            }
        }
    }
    for (xv, yv) in points {
        let x = sketch.vars[xv];
        let y = sketch.vars[yv];
        let dx = x - px;
        let dy = y - py;
        let nx = px + c * dx - s * dy;
        let ny = py + s * dx + c * dy;
        sketch.vars[xv] = nx;
        sketch.vars[yv] = ny;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn rotate_90_deg_about_origin_moves_x_to_y() {
        let mut s = Sketch::new();
        let p = s.add_point(1.0, 0.0);
        rotate(&mut s, &[p], (0.0, 0.0), FRAC_PI_2);
        let pp = s.point_at(p).unwrap();
        let (x, y) = pp.read(&s.vars);
        assert!(x.abs() < 1e-10, "x={x}");
        assert!((y - 1.0).abs() < 1e-10, "y={y}");
    }

    #[test]
    fn rotate_about_self_is_noop() {
        let mut s = Sketch::new();
        let p = s.add_point(3.0, 4.0);
        rotate(&mut s, &[p], (3.0, 4.0), 1.234);
        let pp = s.point_at(p).unwrap();
        let (x, y) = pp.read(&s.vars);
        assert!((x - 3.0).abs() < 1e-12);
        assert!((y - 4.0).abs() < 1e-12);
    }
}
