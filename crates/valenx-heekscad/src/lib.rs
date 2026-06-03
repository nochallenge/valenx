//! # valenx-heekscad
//!
//! HeeksCAD primitives + CAM — Phase 69.
//!
//! Modules:
//! - [`drawing::Drawing`] / [`drawing::Layer`] / [`drawing::HeeksObject`]
//!   — object tree with layer-aware grouping (Sketch / Pad / Pocket /
//!   Drill primitives).
//! - [`cam::pocket_op`] / [`cam::profile_op`] / [`cam::drill_op`] —
//!   HeeksCNC-style toolpath ops.
//! - [`nc_export::write_heeks`] / [`nc_export::write_heeks_string`]
//!   — `.nc` writer matching HeeksCAD's post-processor.
//! - [`persist::to_ron_string`] / [`persist::from_ron_str`] —
//!   workbench-internal RON envelope.
//! - [`panel::HeeksCadPanelState`] — UI envelope.
//! - [`error`] — typed [`HeeksCadError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cam;
pub mod drawing;
pub mod error;
pub mod nc_export;
pub mod panel;
pub mod persist;

pub use cam::{drill_op, pocket_op, profile_op, Move, Tool, Toolpath};
pub use drawing::{Drawing, HeeksObject, Layer, SketchEntity};
pub use error::{ErrorCategory, HeeksCadError};
pub use nc_export::{write_heeks, write_heeks_string};
pub use panel::HeeksCadPanelState;
pub use persist::{from_ron_str, to_ron_string, PanelFile, VERSION};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn tool() -> Tool {
        Tool {
            diameter: 3.0,
            plunge_rate: 100.0,
            feed_rate: 600.0,
        }
    }

    fn unit_square() -> Vec<[f64; 2]> {
        vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
    }

    fn complex_drawing() -> Drawing {
        let mut d = Drawing::new();
        d.layers.push(Layer::new("CAM"));
        d.layer_mut("0").unwrap().objects.push(HeeksObject::Sketch {
            name: "s1".into(),
            plane_normal: Vector3::z(),
            entities: vec![
                SketchEntity::Line {
                    a: [0.0, 0.0],
                    b: [10.0, 0.0],
                },
                SketchEntity::Arc {
                    centre: [5.0, 0.0],
                    radius: 5.0,
                    start: 0.0,
                    sweep: std::f64::consts::PI,
                },
            ],
        });
        d.layer_mut("CAM").unwrap().objects.push(HeeksObject::Pad {
            name: "pad1".into(),
            sketch_ref: "s1".into(),
            height: 5.0,
        });
        d
    }

    #[test]
    fn drawing_default_layer() {
        let d = Drawing::new();
        assert_eq!(d.layers.len(), 1);
        assert_eq!(d.layers[0].name, "0");
    }

    #[test]
    fn drawing_total_objects() {
        let d = complex_drawing();
        assert_eq!(d.total_objects(), 2);
    }

    #[test]
    fn pocket_op_rejects_short_boundary() {
        let r = pocket_op(&[[0.0, 0.0], [1.0, 0.0]], 1.0, tool());
        assert!(matches!(r, Err(HeeksCadError::BadParameter { .. })));
    }

    #[test]
    fn pocket_op_rejects_zero_depth() {
        let r = pocket_op(&unit_square(), 0.0, tool());
        assert!(matches!(r, Err(HeeksCadError::BadParameter { .. })));
    }

    #[test]
    fn pocket_op_emits_plunge_and_retract() {
        let tp = pocket_op(&unit_square(), 2.0, tool()).unwrap();
        assert_eq!(tp.op_name, "pocket");
        assert!(tp.moves.iter().any(|m| matches!(m, Move::Plunge { .. })));
        assert!(tp.moves.iter().any(|m| matches!(m, Move::Retract { .. })));
    }

    #[test]
    fn profile_op_closes_polygon() {
        let tp = profile_op(&unit_square(), 1.0, tool()).unwrap();
        // Last feed move returns to start.
        let last_feed = tp.moves.iter().filter_map(|m| match m {
            Move::Feed { x, y } => Some([*x, *y]),
            _ => None,
        });
        let pts: Vec<[f64; 2]> = last_feed.collect();
        assert_eq!(pts.last(), Some(&[0.0, 0.0]));
    }

    #[test]
    fn drill_op_rejects_empty_positions() {
        let r = drill_op(&[], 5.0, tool());
        assert!(matches!(r, Err(HeeksCadError::BadParameter { .. })));
    }

    #[test]
    fn drill_op_one_rapid_per_hole() {
        let positions = vec![[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]];
        let tp = drill_op(&positions, 3.0, tool()).unwrap();
        let rapids = tp
            .moves
            .iter()
            .filter(|m| matches!(m, Move::Rapid { .. }))
            .count();
        assert_eq!(rapids, 3);
    }

    #[test]
    fn nc_writer_includes_header() {
        let tp = profile_op(&unit_square(), 1.0, tool()).unwrap();
        let s = write_heeks_string(&tp);
        assert!(s.contains("valenx-heekscad"));
        assert!(s.starts_with(";"));
        assert!(s.contains("M30"));
    }

    #[test]
    fn nc_writer_g0_g1_correctly() {
        let tp = profile_op(&unit_square(), 1.0, tool()).unwrap();
        let s = write_heeks_string(&tp);
        assert!(s.contains("G0 X"));
        assert!(s.contains("G1 X"));
    }

    #[test]
    fn ron_round_trip_preserves_drawing() {
        let d = complex_drawing();
        let s = to_ron_string(&d).unwrap();
        let f = from_ron_str(&s).unwrap();
        assert_eq!(f.drawing.layers.len(), d.layers.len());
        assert_eq!(f.drawing.total_objects(), d.total_objects());
    }

    #[test]
    fn panel_add_toolpath_records_status() {
        let mut p = HeeksCadPanelState::new();
        let tp = profile_op(&unit_square(), 1.0, tool()).unwrap();
        p.add_toolpath(tp);
        assert_eq!(p.toolpaths.len(), 1);
        assert!(p.last_status.is_some());
    }
}
