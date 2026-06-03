//! # valenx-vector-graphics
//!
//! Inkscape-style SVG vector primitives — Phase 68.
//!
//! Modules:
//! - [`entity::VectorEntity`] / [`entity::PathSegment`] — six entity
//!   kinds + the SVG `d` mini-language.
//! - [`svg::to_svg`] / [`svg::from_svg`] — round-trip writer + minimal
//!   reader.
//! - [`path::bbox`] / [`path::length`] — bounding box + arc-length.
//! - [`panel::VectorPanelState`] — UI envelope.
//! - [`error`] — typed [`VectorError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod entity;
pub mod error;
pub mod panel;
pub mod path;
pub mod svg;

pub use entity::{PathSegment, VectorEntity};
pub use error::{ErrorCategory, VectorError};
pub use panel::VectorPanelState;
pub use svg::{from_svg, to_svg};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector2;

    fn sample_entities() -> Vec<VectorEntity> {
        vec![
            VectorEntity::Line {
                a: Vector2::new(0.0, 0.0),
                b: Vector2::new(10.0, 10.0),
            },
            VectorEntity::Rect {
                origin: Vector2::new(1.0, 1.0),
                size: Vector2::new(5.0, 3.0),
            },
            VectorEntity::Ellipse {
                centre: Vector2::new(5.0, 5.0),
                rx: 2.0,
                ry: 1.5,
            },
            VectorEntity::Polygon(vec![
                Vector2::new(0.0, 0.0),
                Vector2::new(1.0, 0.0),
                Vector2::new(0.5, 1.0),
            ]),
            VectorEntity::Text {
                anchor: Vector2::new(2.0, 2.0),
                font_size: 12.0,
                text: "hi".into(),
            },
            VectorEntity::Path(vec![
                PathSegment::MoveTo(Vector2::zeros()),
                PathSegment::LineTo(Vector2::new(1.0, 1.0)),
                PathSegment::CurveTo {
                    c1: Vector2::new(2.0, 2.0),
                    c2: Vector2::new(3.0, 0.0),
                    end: Vector2::new(4.0, 4.0),
                },
                PathSegment::Close,
            ]),
        ]
    }

    #[test]
    fn to_svg_contains_each_kind() {
        let s = to_svg(&sample_entities());
        assert!(s.contains("<line"));
        assert!(s.contains("<rect"));
        assert!(s.contains("<ellipse"));
        assert!(s.contains("<polygon"));
        assert!(s.contains("<text"));
        assert!(s.contains("<path"));
    }

    #[test]
    fn round_trip_count_preserved() {
        let e = sample_entities();
        let s = to_svg(&e);
        let back = from_svg(&s).expect("parse ok");
        assert_eq!(back.len(), e.len());
    }

    #[test]
    fn round_trip_preserves_line() {
        let e = vec![VectorEntity::Line {
            a: Vector2::new(1.0, 2.0),
            b: Vector2::new(3.0, 4.0),
        }];
        let s = to_svg(&e);
        let back = from_svg(&s).expect("parse ok");
        assert_eq!(back, e);
    }

    #[test]
    fn round_trip_preserves_path_subset() {
        let e = vec![VectorEntity::Path(vec![
            PathSegment::MoveTo(Vector2::new(0.0, 0.0)),
            PathSegment::LineTo(Vector2::new(1.0, 0.0)),
            PathSegment::Close,
        ])];
        let s = to_svg(&e);
        let back = from_svg(&s).expect("parse ok");
        assert_eq!(back, e);
    }

    #[test]
    fn bbox_simple_path() {
        let p = vec![
            PathSegment::MoveTo(Vector2::new(0.0, 0.0)),
            PathSegment::LineTo(Vector2::new(10.0, 5.0)),
        ];
        let (lo, hi) = path::bbox(&p);
        assert_eq!(lo, Vector2::zeros());
        assert_eq!(hi, Vector2::new(10.0, 5.0));
    }

    #[test]
    fn length_straight_segment() {
        let p = vec![
            PathSegment::MoveTo(Vector2::zeros()),
            PathSegment::LineTo(Vector2::new(3.0, 4.0)),
        ];
        let l = path::length(&p);
        assert!((l - 5.0).abs() < 1e-9);
    }

    #[test]
    fn length_close_returns_to_start() {
        let p = vec![
            PathSegment::MoveTo(Vector2::zeros()),
            PathSegment::LineTo(Vector2::new(3.0, 0.0)),
            PathSegment::LineTo(Vector2::new(3.0, 4.0)),
            PathSegment::Close,
        ];
        let l = path::length(&p);
        // 3 + 4 + 5 = 12.
        assert!((l - 12.0).abs() < 1e-9);
    }

    #[test]
    fn parse_rejects_unterminated_tag() {
        let r = from_svg("<line x1=\"0\"");
        assert!(matches!(r, Err(VectorError::Parse { .. })));
    }

    #[test]
    fn panel_set_entities_records_status() {
        let mut p = VectorPanelState::new();
        p.set_entities(sample_entities());
        assert_eq!(p.entities.len(), 6);
        assert!(p.last_status.is_some());
    }
}
