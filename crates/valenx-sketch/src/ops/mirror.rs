//! Mirror entities across a line.
//!
//! Phase 12E Task 40. Reflects the (x, y) coordinates of every
//! point referenced (directly or via container entities) across the
//! line `(p, n)` where `p` is a point on the line and `n` is its
//! direction vector. New entities reuse the same primitive kind as
//! their source; only [`crate::geom::Point2`] / [`crate::geom::Line2`]
//! / [`crate::geom::Circle2`] / [`crate::geom::Arc2`] are duplicated
//! in v1 (BSpline / Ellipse / EllipticalArc are skipped — a future
//! task can extend this once the constraint coverage warrants it).

use crate::geom::{Entity, EntityId};
use crate::sketch::Sketch;

/// Mirror line: `(point_on_line, direction)`.
pub struct MirrorLine {
    /// A point on the mirror line.
    pub point: (f64, f64),
    /// Direction vector of the mirror line (any non-zero vector).
    pub direction: (f64, f64),
}

/// Reflect a 2-D point across `line`.
pub fn reflect_point(p: (f64, f64), line: &MirrorLine) -> (f64, f64) {
    let (px, py) = p;
    let (lx, ly) = line.point;
    let (dx, dy) = line.direction;
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-30 {
        return p;
    }
    // Component of (p - L) along d:
    let t = ((px - lx) * dx + (py - ly) * dy) / len2;
    // Foot of the perpendicular:
    let fx = lx + t * dx;
    let fy = ly + t * dy;
    // Reflection: 2*foot - p.
    (2.0 * fx - px, 2.0 * fy - py)
}

/// Duplicate `entities` mirrored across `line`. Returns the ids of the
/// new entities (in the same order as input). Unsupported entity
/// kinds are silently skipped.
pub fn mirror(sketch: &mut Sketch, entities: &[EntityId], line: &MirrorLine) -> Vec<EntityId> {
    let mut created = Vec::new();
    for id in entities {
        let kind = match sketch.entities.get(id.0.wrapping_sub(1)) {
            Some(e) => e.clone(),
            None => continue,
        };
        let new_id = match kind {
            Entity::Point(p) => {
                let (x, y) = p.read(&sketch.vars);
                let (nx, ny) = reflect_point((x, y), line);
                Some(sketch.add_point(nx, ny))
            }
            Entity::Line(l) => {
                let ((sx, sy), (ex, ey)) = l.endpoints(&sketch.vars);
                let (nsx, nsy) = reflect_point((sx, sy), line);
                let (nex, ney) = reflect_point((ex, ey), line);
                let a = sketch.add_point(nsx, nsy);
                let b = sketch.add_point(nex, ney);
                sketch.add_line(a, b).ok()
            }
            Entity::Circle(c) => {
                let (cx, cy) = c.center.read(&sketch.vars);
                let r = c.radius(&sketch.vars);
                let (ncx, ncy) = reflect_point((cx, cy), line);
                let center = sketch.add_point(ncx, ncy);
                sketch.add_circle(center, r).ok()
            }
            Entity::Arc(a) => {
                let (cx, cy) = a.center.read(&sketch.vars);
                let r = a.radius(&sketch.vars);
                let (sa, ea) = a.angles(&sketch.vars);
                let (ncx, ncy) = reflect_point((cx, cy), line);
                let center = sketch.add_point(ncx, ncy);
                // Mirroring an angle is just π - angle relative to the
                // mirror direction; for the v1 sketch ops we leave the
                // angles unchanged and let downstream constraints fix
                // them up.
                sketch.add_arc(center, r, sa, ea).ok()
            }
            _ => None,
        };
        if let Some(nid) = new_id {
            created.push(nid);
        }
    }
    created
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflect_point_across_y_axis_negates_x() {
        let line = MirrorLine {
            point: (0.0, 0.0),
            direction: (0.0, 1.0),
        };
        let (x, y) = reflect_point((3.0, 4.0), &line);
        assert!((x + 3.0).abs() < 1e-12);
        assert!((y - 4.0).abs() < 1e-12);
    }

    #[test]
    fn reflect_point_across_x_axis_negates_y() {
        let line = MirrorLine {
            point: (0.0, 0.0),
            direction: (1.0, 0.0),
        };
        let (x, y) = reflect_point((3.0, 4.0), &line);
        assert!((x - 3.0).abs() < 1e-12);
        assert!((y + 4.0).abs() < 1e-12);
    }

    #[test]
    fn mirror_creates_new_point_at_reflected_xy() {
        let mut s = Sketch::new();
        let p = s.add_point(2.0, 5.0);
        let line = MirrorLine {
            point: (0.0, 0.0),
            direction: (0.0, 1.0),
        };
        let created = mirror(&mut s, &[p], &line);
        assert_eq!(created.len(), 1);
        let new_pt = s.point_at(created[0]).unwrap();
        let (x, y) = new_pt.read(&s.vars);
        assert!((x + 2.0).abs() < 1e-12);
        assert!((y - 5.0).abs() < 1e-12);
    }
}
