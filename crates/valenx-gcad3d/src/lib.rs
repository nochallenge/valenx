//! # valenx-gcad3d
//!
//! gCAD3D parametric primitives — Phase 66.
//!
//! Modules:
//! - [`plane::Plane3d`] — `xy_at` / `xz_at` / `yz_at` constructors.
//! - [`line::Line3d`] — `between` / `parallel_to` / `perpendicular_at`.
//! - [`arc::Arc3d`] — `three_point` / `tangent_to`.
//! - [`surface::ruled`] — ruled-surface descriptor.
//! - [`text::extrude`] — fixed-pitch 3x5 block-letter extruded text.
//! - [`text::extrude_strokes`] — vector stroke font (Phase 66.5).
//! - [`panel::Gcad3dPanelState`] — UI envelope.
//! - [`error`] — typed [`Gcad3dError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod arc;
pub mod error;
pub mod line;
pub mod panel;
pub mod plane;
pub mod surface;
pub mod text;

pub use arc::Arc3d;
pub use error::{ErrorCategory, Gcad3dError};
pub use line::Line3d;
pub use panel::Gcad3dPanelState;
pub use plane::Plane3d;
pub use surface::{ruled, RuledSurface};
pub use text::{
    extrude as extrude_text, extrude_strokes as extrude_text_strokes, Glyph, StrokeGlyph,
    StrokeText, TextSolid,
};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn plane_constructors_correct_normal() {
        assert_eq!(Plane3d::xy_at(1.0).normal, Vector3::z());
        assert_eq!(Plane3d::xz_at(2.0).normal, Vector3::y());
        assert_eq!(Plane3d::yz_at(3.0).normal, Vector3::x());
    }

    #[test]
    fn line_between_rejects_coincident() {
        let r = Line3d::between(Vector3::zeros(), Vector3::zeros());
        assert!(matches!(r, Err(Gcad3dError::Degenerate(_))));
    }

    #[test]
    fn line_between_normalises_direction() {
        let l = Line3d::between(Vector3::zeros(), Vector3::new(3.0, 4.0, 0.0)).unwrap();
        assert!((l.direction.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn line_parallel_to_keeps_direction() {
        let l = Line3d::between(Vector3::zeros(), Vector3::x()).unwrap();
        let p = Line3d::parallel_to(&l, 0.5).unwrap();
        assert!((p.direction.norm() - 1.0).abs() < 1e-12);
        assert!((p.direction - l.direction).norm() < 1e-12);
    }

    #[test]
    fn line_perpendicular_at_rejects_on_line() {
        let l = Line3d::between(Vector3::zeros(), Vector3::x()).unwrap();
        let r = Line3d::perpendicular_at(&l, Vector3::new(1.0, 0.0, 0.0));
        assert!(matches!(r, Err(Gcad3dError::Degenerate(_))));
    }

    #[test]
    fn arc_three_point_unit_circle_in_xy() {
        let p1 = Vector3::new(1.0, 0.0, 0.0);
        let p2 = Vector3::new(0.0, 1.0, 0.0);
        let p3 = Vector3::new(-1.0, 0.0, 0.0);
        let a = Arc3d::three_point(p1, p2, p3).unwrap();
        assert!((a.radius - 1.0).abs() < 1e-9);
        assert!((a.centre - Vector3::zeros()).norm() < 1e-9);
    }

    #[test]
    fn arc_three_point_rejects_colinear() {
        let p1 = Vector3::zeros();
        let p2 = Vector3::x();
        let p3 = Vector3::new(2.0, 0.0, 0.0);
        assert!(matches!(
            Arc3d::three_point(p1, p2, p3),
            Err(Gcad3dError::Degenerate(_))
        ));
    }

    #[test]
    fn arc_tangent_to_radius_correct() {
        let l = Line3d::between(Vector3::zeros(), Vector3::x()).unwrap();
        let a = Arc3d::tangent_to(&l, 2.0).unwrap();
        assert_eq!(a.radius, 2.0);
    }

    #[test]
    fn arc_tangent_rejects_negative_radius() {
        let l = Line3d::between(Vector3::zeros(), Vector3::x()).unwrap();
        assert!(matches!(
            Arc3d::tangent_to(&l, -1.0),
            Err(Gcad3dError::BadParameter { .. })
        ));
    }

    #[test]
    fn surface_ruled_descriptor() {
        let s = surface::ruled("c1", "c2");
        assert_eq!(s.curve1, "c1");
        assert_eq!(s.curve2, "c2");
    }

    #[test]
    fn text_extrude_one_letter_has_cells() {
        let t = text::extrude("A", 10.0, 1.0).unwrap();
        assert_eq!(t.glyphs.len(), 1);
        assert!(!t.glyphs[0].cells.is_empty());
    }

    #[test]
    fn text_extrude_rejects_bad_font_size() {
        assert!(matches!(
            text::extrude("A", 0.0, 1.0),
            Err(Gcad3dError::BadParameter { .. })
        ));
    }

    #[test]
    fn text_extrude_rejects_unsupported_char() {
        assert!(matches!(
            text::extrude("@", 10.0, 1.0),
            Err(Gcad3dError::UnsupportedChar('@'))
        ));
    }

    #[test]
    fn text_extrude_space_advances_cursor() {
        let t = text::extrude("A A", 10.0, 1.0).unwrap();
        // Two letters extracted; the space is gap-only.
        assert_eq!(t.glyphs.len(), 2);
    }

    #[test]
    fn panel_add_arc_records_status() {
        let mut p = Gcad3dPanelState::new();
        let a = Arc3d::three_point(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
        )
        .unwrap();
        p.add_arc(a);
        assert_eq!(p.arcs.len(), 1);
        assert!(p.last_status.is_some());
    }
}
