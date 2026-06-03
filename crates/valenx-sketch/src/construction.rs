//! Construction-geometry flag for sketch entities.
//!
//! Phase 12C. A construction entity is visible in the sketch overlay
//! (rendered as a dashed light line) but is excluded from the
//! extruded profile. It is solver-visible — constraints can still
//! reference construction geometry and the solver still drives it.
//!
//! Storage: a parallel `Vec<bool>` on [`crate::sketch::Sketch`] sized
//! to match [`crate::sketch::Sketch::entities`]. Adding entities
//! pushes `false` (regular). [`crate::sketch::Sketch::toggle_construction`]
//! flips the flag for a given [`crate::geom::EntityId`].

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Returns `true` if the entity is flagged as construction geometry.
///
/// Returns `false` for unknown ids — the caller is expected to have
/// validated `id` against the sketch first.
pub fn is_construction(sketch: &Sketch, id: EntityId) -> bool {
    let idx = id.0.wrapping_sub(1);
    sketch.construction.get(idx).copied().unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_entities_are_regular_not_construction() {
        let mut s = Sketch::new();
        let p = s.add_point(0.0, 0.0);
        assert!(!is_construction(&s, p));
    }

    #[test]
    fn toggle_construction_flips_flag() {
        let mut s = Sketch::new();
        let p = s.add_point(0.0, 0.0);
        s.toggle_construction(p);
        assert!(is_construction(&s, p));
        s.toggle_construction(p);
        assert!(!is_construction(&s, p));
    }

    #[test]
    fn is_construction_for_unknown_id_returns_false() {
        let s = Sketch::new();
        assert!(!is_construction(&s, EntityId(99)));
    }
}
