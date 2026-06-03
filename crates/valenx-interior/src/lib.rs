//! # valenx-interior
//!
//! Sweet Home 3D-style interior-design workflow — Phase 62.
//!
//! Modules:
//! - [`room::Room`] / [`room::WallRef`] / [`room::compute_area`] —
//!   room data model + polygon area.
//! - [`furniture::Furniture`] — 12-kind catalog.
//! - [`furniture::to_solid`] — parametric box geometry per kind.
//! - [`furniture::Placement`] / [`Furniture::place`] — placement
//!   with position + rotation.
//! - [`panel::InteriorPanelState`] — room list + furniture palette.
//! - [`error`] — typed [`InteriorError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod furniture;
pub mod panel;
pub mod room;

pub use error::{ErrorCategory, InteriorError};
pub use furniture::{to_solid, Furniture, Placement, Solid};
pub use panel::InteriorPanelState;
pub use room::{compute_area, Room, WallRef};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Vector2, Vector3};

    #[test]
    fn catalog_has_twelve_kinds() {
        assert_eq!(Furniture::all().len(), 12);
    }

    #[test]
    fn all_kinds_have_distinct_names() {
        let names: std::collections::BTreeSet<&str> =
            Furniture::all().iter().map(|f| f.name()).collect();
        assert_eq!(names.len(), 12);
    }

    #[test]
    fn default_sizes_strictly_positive() {
        for f in Furniture::all() {
            let s = f.default_size();
            assert!(s.x > 0.0 && s.y > 0.0 && s.z > 0.0, "{f:?} has bad size");
        }
    }

    #[test]
    fn to_solid_tags_kind() {
        let s = to_solid(Furniture::Bed, Furniture::Bed.default_size());
        assert_eq!(s.kind, Furniture::Bed);
        assert_eq!(s.size, Furniture::Bed.default_size());
        assert_eq!(s.origin, Vector3::zeros());
    }

    #[test]
    fn place_carries_room_id() {
        let p = Furniture::Chair.place(Vector3::new(1.0, 2.0, 0.0), 0.5, "kitchen");
        assert_eq!(p.kind, Furniture::Chair);
        assert_eq!(p.room_id, "kitchen");
        assert_eq!(p.rotation_rad, 0.5);
    }

    #[test]
    fn area_unit_square() {
        let poly = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 0.0),
            Vector2::new(1.0, 1.0),
            Vector2::new(0.0, 1.0),
        ];
        assert!((compute_area(&poly) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn area_degenerate_zero() {
        let poly = vec![Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.0)];
        assert_eq!(compute_area(&poly), 0.0);
    }

    #[test]
    fn panel_add_room_rejects_duplicate() {
        let mut p = InteriorPanelState::new();
        p.add_room(Room::new("a", "Kitchen", 2.5)).unwrap();
        let r = p.add_room(Room::new("a", "Bath", 2.5));
        assert!(matches!(r, Err(InteriorError::BadParameter { .. })));
    }

    #[test]
    fn panel_click_to_place_requires_selection() {
        let mut p = InteriorPanelState::new();
        p.add_room(Room::new("k", "Kitchen", 2.5)).unwrap();
        let r = p.click_to_place(Vector3::zeros(), "k");
        assert!(matches!(r, Err(InteriorError::BadParameter { .. })));
    }

    #[test]
    fn panel_click_to_place_rejects_unknown_room() {
        let mut p = InteriorPanelState::new();
        p.select(Furniture::Chair);
        let r = p.click_to_place(Vector3::zeros(), "nope");
        assert!(matches!(r, Err(InteriorError::UnknownRoom(_))));
    }

    #[test]
    fn panel_happy_path_places_one_item() {
        let mut p = InteriorPanelState::new();
        p.add_room(Room::new("k", "Kitchen", 2.5)).unwrap();
        p.select(Furniture::Table);
        p.click_to_place(Vector3::new(1.0, 1.0, 0.0), "k").unwrap();
        assert_eq!(p.placements.len(), 1);
        assert_eq!(p.placements[0].kind, Furniture::Table);
    }
}
