//! # valenx-librecad-2d
//!
//! LibreCAD-style 2D drafting workbench — Phase 60. Pure 2D, full DXF
//! round-trip (LINE / CIRCLE / ARC / LWPOLYLINE / SPLINE / TEXT /
//! MTEXT / DIMENSION / HATCH / INSERT) plus layers + blocks.
//!
//! # Why a separate crate
//!
//! Phase 5's `valenx-draft` already does DXF *writing* for 2D
//! projections of 3D parts. This crate is the **2D-only** workbench
//! that consumes DXF as a first-class input format and exposes a 2D
//! drafting data model independent of the 3D pipeline.
//!
//! Module map:
//! - [`drawing`] — `Drawing2D`, `Layer`, `Block`, `Entity2D`.
//! - [`dxf`] — `read_full`, `write_full`, `parse`, `serialise`.
//! - [`panel`] — UI panel state.
//! - [`persist`] — RON envelope (workbench-internal format).
//! - [`error`] — typed [`LibreCadError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod drawing;
pub mod dxf;
pub mod error;
pub mod geometry;
pub mod panel;
pub mod persist;

pub use drawing::{Block, Drawing2D, Entity2D, Layer};
pub use error::{ErrorCategory, LibreCadError};
pub use geometry::{bounding_box_diagonal_2d, polygon_area, polyline_length};
pub use panel::LibreCadPanelState;
pub use persist::{from_ron_str, to_ron_string, PanelFile, VERSION};

#[cfg(test)]
mod tests {
    use super::*;

    fn complex_drawing() -> Drawing2D {
        let mut d = Drawing2D::new();
        d.layers.push(Layer {
            name: "WALLS".into(),
            color: 1,
            linetype: "CONTINUOUS".into(),
            visible: true,
        });
        d.entities.push(Entity2D::Line {
            layer: "WALLS".into(),
            a: [0.0, 0.0],
            b: [10.0, 0.0],
        });
        d.entities.push(Entity2D::Circle {
            layer: "0".into(),
            centre: [5.0, 5.0],
            radius: 2.5,
        });
        d.entities.push(Entity2D::Arc {
            layer: "0".into(),
            centre: [0.0, 0.0],
            radius: 5.0,
            start_angle_deg: 0.0,
            end_angle_deg: 90.0,
        });
        d.entities.push(Entity2D::Polyline {
            layer: "0".into(),
            vertices: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            closed: true,
        });
        d.entities.push(Entity2D::Text {
            layer: "0".into(),
            position: [3.0, 3.0],
            height: 0.5,
            text: "HELLO".into(),
        });
        d.entities.push(Entity2D::Insert {
            layer: "0".into(),
            block: "BLK1".into(),
            position: [10.0, 10.0],
            scale: 1.0,
            rotation_deg: 45.0,
        });
        d.blocks.push(Block {
            name: "BLK1".into(),
            origin: [0.0, 0.0],
            entities: vec![Entity2D::Line {
                layer: "0".into(),
                a: [0.0, 0.0],
                b: [1.0, 1.0],
            }],
        });
        d
    }

    #[test]
    fn empty_drawing_has_default_layer() {
        let d = Drawing2D::new();
        assert_eq!(d.layers.len(), 1);
        assert_eq!(d.layers[0].name, "0");
    }

    #[test]
    fn entity_kind_dispatches_correctly() {
        let e = Entity2D::Line {
            layer: "0".into(),
            a: [0.0, 0.0],
            b: [1.0, 1.0],
        };
        assert_eq!(e.kind(), "LINE");
        assert_eq!(e.layer(), "0");
    }

    #[test]
    fn dxf_round_trip_preserves_entity_count() {
        let d = complex_drawing();
        let s = dxf::serialise(&d);
        let d2 = dxf::parse(&s).expect("parse ok");
        assert_eq!(d.entities.len(), d2.entities.len());
        assert_eq!(d.blocks.len(), d2.blocks.len());
    }

    #[test]
    fn dxf_round_trip_preserves_line_endpoints() {
        let mut d = Drawing2D::new();
        d.entities.push(Entity2D::Line {
            layer: "0".into(),
            a: [1.0, 2.0],
            b: [3.0, 4.0],
        });
        let s = dxf::serialise(&d);
        let d2 = dxf::parse(&s).expect("parse ok");
        match &d2.entities[0] {
            Entity2D::Line { a, b, .. } => {
                assert_eq!(*a, [1.0, 2.0]);
                assert_eq!(*b, [3.0, 4.0]);
            }
            other => panic!("expected Line, got {other:?}"),
        }
    }

    #[test]
    fn dxf_round_trip_preserves_dimension_points() {
        // Regression: the serialiser wrote `text_pos` on group codes 10/20 —
        // the same codes the parser reads into definition point `a` — so on
        // reload the second 10/20 clobbered `a`. `text_pos` now uses 13/23.
        let mut d = Drawing2D::new();
        d.entities.push(Entity2D::Dimension {
            layer: "0".into(),
            a: [1.0, 2.0],
            b: [3.0, 4.0],
            text_pos: [9.0, 8.0],
            text: "5.0".into(),
        });
        let s = dxf::serialise(&d);
        let d2 = dxf::parse(&s).expect("parse ok");
        match &d2.entities[0] {
            Entity2D::Dimension {
                a, b, text_pos, text, ..
            } => {
                assert_eq!(*a, [1.0, 2.0], "definition point a must survive reload");
                assert_eq!(*b, [3.0, 4.0]);
                assert_eq!(*text_pos, [9.0, 8.0], "text midpoint must round-trip");
                assert_eq!(text, "5.0");
            }
            other => panic!("expected Dimension, got {other:?}"),
        }
    }

    #[test]
    fn dxf_round_trip_preserves_polyline_vertices() {
        let mut d = Drawing2D::new();
        d.entities.push(Entity2D::Polyline {
            layer: "0".into(),
            vertices: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            closed: true,
        });
        let s = dxf::serialise(&d);
        let d2 = dxf::parse(&s).expect("parse ok");
        match &d2.entities[0] {
            Entity2D::Polyline { vertices, closed, .. } => {
                assert_eq!(vertices.len(), 3);
                assert!(*closed);
            }
            other => panic!("expected Polyline, got {other:?}"),
        }
    }

    #[test]
    fn dxf_parse_rejects_garbage_group_code() {
        let bad = "not_a_number\nLINE\n";
        let res = dxf::parse(bad);
        assert!(matches!(res, Err(LibreCadError::DxfParse { .. })));
    }

    #[test]
    fn ron_round_trip_preserves_drawing() {
        let d = complex_drawing();
        let s = to_ron_string(&d).unwrap();
        let f = from_ron_str(&s).unwrap();
        assert_eq!(f.drawing.entities.len(), d.entities.len());
        assert_eq!(f.drawing.blocks.len(), d.blocks.len());
        assert_eq!(f.drawing.layers.len(), d.layers.len());
    }
}
