//! Copy entities to a new translated position.
//!
//! Phase 12E Task 41.

use crate::geom::{Entity, EntityId};
use crate::sketch::Sketch;

/// Duplicate `entities` translated by `offset = (dx, dy)`. Returns
/// the new entity ids. Unsupported kinds (BSpline / Ellipse /
/// EllipticalArc) are skipped silently in v1.
pub fn copy(sketch: &mut Sketch, entities: &[EntityId], offset: (f64, f64)) -> Vec<EntityId> {
    let (dx, dy) = offset;
    let mut created = Vec::new();
    for id in entities {
        let kind = match sketch.entities.get(id.0.wrapping_sub(1)) {
            Some(e) => e.clone(),
            None => continue,
        };
        let new_id = match kind {
            Entity::Point(p) => {
                let (x, y) = p.read(&sketch.vars);
                Some(sketch.add_point(x + dx, y + dy))
            }
            Entity::Line(l) => {
                let ((sx, sy), (ex, ey)) = l.endpoints(&sketch.vars);
                let a = sketch.add_point(sx + dx, sy + dy);
                let b = sketch.add_point(ex + dx, ey + dy);
                sketch.add_line(a, b).ok()
            }
            Entity::Circle(c) => {
                let (cx, cy) = c.center.read(&sketch.vars);
                let r = c.radius(&sketch.vars);
                let center = sketch.add_point(cx + dx, cy + dy);
                sketch.add_circle(center, r).ok()
            }
            Entity::Arc(a) => {
                let (cx, cy) = a.center.read(&sketch.vars);
                let r = a.radius(&sketch.vars);
                let (sa, ea) = a.angles(&sketch.vars);
                let center = sketch.add_point(cx + dx, cy + dy);
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
    fn copy_translates_point() {
        let mut s = Sketch::new();
        let p = s.add_point(1.0, 2.0);
        let created = copy(&mut s, &[p], (10.0, 20.0));
        assert_eq!(created.len(), 1);
        let np = s.point_at(created[0]).unwrap();
        let (x, y) = np.read(&s.vars);
        assert!((x - 11.0).abs() < 1e-12);
        assert!((y - 22.0).abs() < 1e-12);
    }

    #[test]
    fn copy_preserves_circle_radius() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let cir = s.add_circle(c, 3.5).unwrap();
        let created = copy(&mut s, &[cir], (5.0, 0.0));
        assert_eq!(created.len(), 1);
        let nc = s.circle_at(created[0]).unwrap();
        assert!((nc.radius(&s.vars) - 3.5).abs() < 1e-12);
    }
}
