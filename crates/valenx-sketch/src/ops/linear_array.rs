//! N copies of entities translated along a direction.
//!
//! Phase 12E Task 45.

use super::copy;
use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Create `count - 1` translated copies along `direction` (unit vector
/// scaled by `spacing`). The original entities are *not* duplicated.
pub fn linear_array(
    sketch: &mut Sketch,
    entities: &[EntityId],
    direction: (f64, f64),
    count: usize,
    spacing: f64,
) -> Vec<EntityId> {
    if count < 2 {
        return Vec::new();
    }
    let (dx, dy) = direction;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-30 {
        return Vec::new();
    }
    let ux = dx / len;
    let uy = dy / len;
    let mut created = Vec::new();
    for i in 1..count {
        let off = ((ux * spacing) * i as f64, (uy * spacing) * i as f64);
        let new_ids = copy::copy(sketch, entities, off);
        created.extend(new_ids);
    }
    created
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_array_3_along_x_at_spacing_1() {
        let mut s = Sketch::new();
        let p = s.add_point(0.0, 0.0);
        let created = linear_array(&mut s, &[p], (1.0, 0.0), 3, 1.0);
        assert_eq!(created.len(), 2);
        let p1 = s.point_at(created[0]).unwrap();
        let p2 = s.point_at(created[1]).unwrap();
        let (x1, _) = p1.read(&s.vars);
        let (x2, _) = p2.read(&s.vars);
        assert!((x1 - 1.0).abs() < 1e-12);
        assert!((x2 - 2.0).abs() < 1e-12);
    }

    #[test]
    fn linear_array_normalises_direction_vector() {
        let mut s = Sketch::new();
        let p = s.add_point(0.0, 0.0);
        // (3, 4) has magnitude 5; spacing 5 should move (3, 4).
        let created = linear_array(&mut s, &[p], (3.0, 4.0), 2, 5.0);
        assert_eq!(created.len(), 1);
        let np = s.point_at(created[0]).unwrap();
        let (x, y) = np.read(&s.vars);
        assert!((x - 3.0).abs() < 1e-10);
        assert!((y - 4.0).abs() < 1e-10);
    }
}
