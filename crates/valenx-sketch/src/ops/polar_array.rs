//! N copies of entities rotated about a centre point.
//!
//! Phase 12E Task 44.

use super::{copy, rotate};
use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Create `count - 1` rotated copies (the original is included as the
/// "0th" copy and is *not* duplicated). Returns the new entity ids.
/// `total_angle` is the angle swept by all copies in radians.
pub fn polar_array(
    sketch: &mut Sketch,
    entities: &[EntityId],
    center: (f64, f64),
    count: usize,
    total_angle: f64,
) -> Vec<EntityId> {
    if count < 2 {
        return Vec::new();
    }
    let step = total_angle / (count - 1) as f64;
    let mut created = Vec::new();
    // Copy step at a time: copy entities, then rotate the new copies.
    for i in 1..count {
        // Phase 12 v1: copy the originals each step (preserves their
        // original positions), then rotate the new copies by i * step.
        let new_ids = copy::copy(sketch, entities, (0.0, 0.0));
        rotate::rotate(sketch, &new_ids, center, step * i as f64);
        created.extend(new_ids);
    }
    created
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn polar_array_with_count_2_creates_one_copy() {
        let mut s = Sketch::new();
        let p = s.add_point(1.0, 0.0);
        let created = polar_array(&mut s, &[p], (0.0, 0.0), 2, PI);
        assert_eq!(created.len(), 1);
        let np = s.point_at(created[0]).unwrap();
        let (x, y) = np.read(&s.vars);
        // rotated by step = PI / 1 = PI radians
        assert!((x + 1.0).abs() < 1e-10);
        assert!(y.abs() < 1e-10);
    }

    #[test]
    fn polar_array_quarter_at_90_deg_steps() {
        let mut s = Sketch::new();
        let p = s.add_point(1.0, 0.0);
        // Step = FRAC_PI_2 / 1 = FRAC_PI_2 — single rotation 90 deg.
        let created = polar_array(&mut s, &[p], (0.0, 0.0), 2, FRAC_PI_2);
        assert_eq!(created.len(), 1);
        let np = s.point_at(created[0]).unwrap();
        let (x, y) = np.read(&s.vars);
        assert!(x.abs() < 1e-10);
        assert!((y - 1.0).abs() < 1e-10);
    }
}
