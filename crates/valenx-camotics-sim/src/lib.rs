//! # valenx-camotics-sim
//!
//! CAMotics-style **animated** material-removal simulation — Phase 56.
//!
//! Wraps the Phase 17 voxel grid + cut algorithm in an [`animation::Animation`]
//! type that emits per-frame meshes plus per-frame metadata (cut volume,
//! tool position, material-removal rate).
//!
//! # Design
//!
//! [CAMotics](https://camotics.org) renders machining as an animation:
//! show the stock at `t=0`, the finished part at `t=1`, and interpolate
//! between by replaying the toolpath. We mirror that flow:
//!
//! - [`animation::Animation`] owns stock + toolpath + frame budget.
//! - [`animation::Animation::frame(t)`] returns the **faceted** mesh
//!   of the stock after replaying the first `round(t * (n_moves - 1))`
//!   moves.
//! - [`animation::Animation::frame_smooth(t)`] returns the same frame
//!   extracted with **Surface Nets** — a smooth, dual-contoured
//!   surface that rounds off the voxel stair-stepping (Phase 56.5).
//! - [`animation::Animation::frames()`] yields every (faceted) frame
//!   for sequential playback.
//! - [`report::FrameMetadata`] records cut volume + MRR + tool
//!   position per frame; [`report::MaterialRemovalReport`] aggregates.
//! - [`panel::CamoticsPanelState`] wires everything for the
//!   workbench UI envelope (no FileDialog, no live windowing — that
//!   lives in `valenx-app`).
//!
//! ## Why re-run from scratch per frame
//!
//! Voxel state isn't snapshottable without a deep clone, and Phase 17's
//! [`valenx_cam::voxel::Voxel::cut_segment`] is destructive. The v1
//! simulator pays the `n_frames * total_moves` cost for clarity. A
//! delta-cache is future work.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod animation;
pub mod error;
pub mod panel;
pub mod persist;
pub mod report;

pub use animation::Animation;
pub use error::{CamoticsError, ErrorCategory};
pub use panel::CamoticsPanelState;
pub use persist::{from_ron_str, to_ron_string, PanelFile, VERSION};
pub use report::{FrameMetadata, MaterialRemovalReport};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_cam::stock::Stock;
    use valenx_cam::toolpath::{Move, MoveKind, Toolpath};

    fn demo_animation() -> Animation {
        let stock = Stock::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(20.0, 20.0, 10.0),
            "aluminum",
        )
        .expect("stock ok");
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 5.0, 5.0), 0.0));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(20.0, 5.0, 5.0),
            300.0,
        ));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(20.0, 15.0, 5.0),
            300.0,
        ));
        Animation::new(stock, tp, 1.5, 5, (20, 20, 10)).expect("anim ok")
    }

    #[test]
    fn frame_zero_is_full_stock() {
        let a = demo_animation();
        let m = a.frame(0.0).unwrap();
        assert!(m.total_elements() > 0, "frame 0 mesh is empty");
    }

    #[test]
    fn frame_smooth_shares_vertices_unlike_faceted() {
        // After the toolpath cuts the stock, the Surface-Nets frame
        // must be non-empty. For a binary occupancy grid it carries the
        // *same* triangle count as the faceted frame — naive Surface
        // Nets emits one quad per sign-changing grid edge, in 1:1
        // correspondence with the faceted extractor's one quad per
        // boundary voxel face. Its real win is **vertex sharing**: one
        // vertex per boundary cell shared by up to 4 quads, where the
        // faceted extractor pushes 4 fresh un-welded nodes per quad —
        // so the smooth frame carries far fewer vertices.
        let a = demo_animation();
        let faceted = a.frame(1.0).unwrap();
        let smooth = a.frame_smooth(1.0).unwrap();
        assert!(
            smooth.total_elements() > 0,
            "smooth frame should not be empty"
        );
        assert_eq!(
            smooth.total_elements(),
            faceted.total_elements(),
            "binary-grid surface nets emits one quad per sign-changing edge"
        );
        assert!(
            smooth.nodes.len() < faceted.nodes.len(),
            "surface-nets frame should share vertices: smooth {} vs faceted {}",
            smooth.nodes.len(),
            faceted.nodes.len()
        );
    }

    #[test]
    fn frames_returns_n_frames() {
        let a = demo_animation();
        let fs = a.frames().unwrap();
        assert_eq!(fs.len(), 5);
    }

    #[test]
    fn report_has_monotone_cut_volume() {
        let a = demo_animation();
        let r = a.material_removal_report().unwrap();
        for window in r.frames.windows(2) {
            assert!(
                window[1].cut_volume_mm3 >= window[0].cut_volume_mm3 - 1e-9,
                "cut volume regressed: {} -> {}",
                window[0].cut_volume_mm3,
                window[1].cut_volume_mm3
            );
        }
    }

    #[test]
    fn invalid_inputs_rejected() {
        let stock = Stock::new(Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), "alu").unwrap();
        let tp = Toolpath::new();
        assert!(matches!(
            Animation::new(stock.clone(), tp.clone(), 0.0, 5, (10, 10, 10)),
            Err(CamoticsError::BadParameter {
                name: "tool_radius_mm",
                ..
            })
        ));
        assert!(matches!(
            Animation::new(stock.clone(), tp.clone(), 1.0, 1, (10, 10, 10)),
            Err(CamoticsError::BadParameter {
                name: "n_frames",
                ..
            })
        ));
        assert!(matches!(
            Animation::new(stock, tp, 1.0, 5, (0, 10, 10)),
            Err(CamoticsError::BadParameter {
                name: "voxel_resolution",
                ..
            })
        ));
    }

    #[test]
    fn frame_out_of_range_errors() {
        let a = demo_animation();
        assert!(matches!(
            a.frame_at(5),
            Err(CamoticsError::FrameOutOfRange(5, 5))
        ));
    }

    #[test]
    fn panel_round_trip() {
        let a = demo_animation();
        let s = to_ron_string(&Some(a.clone()), 0.5).unwrap();
        let f = from_ron_str(&s).unwrap();
        assert_eq!(f.current_t, 0.5);
        assert!(f.animation.is_some());
    }
}
