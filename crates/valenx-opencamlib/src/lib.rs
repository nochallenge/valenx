//! # valenx-opencamlib
//!
//! Pure-Rust port of [OpenCamLib]'s most-used CAM algorithms — Phase 57.
//!
//! [OpenCamLib]: https://github.com/aewallin/opencamlib
//!
//! # Scope
//!
//! - [`cutter::DropCutter`] — find the Z at which a cylindrical tool
//!   just touches a triangle soup above `(x, y)`.
//! - [`cutter::AdaptiveDropCutter`] — drop-cutter sampling on a grid
//!   backed by an [`octree::Octree`] for skip-empty optimisation.
//! - [`cutter::WaterlinePathPlanner`] — emit constant-Z waterlines
//!   (v1: raw intersection points, no edge stitching).
//! - [`cutter::PushCutter`] — push a tool laterally; record the
//!   resulting Z trail.
//! - [`cutter::EdgeCutter`] — pure edge contact, no face-interior
//!   interaction.
//! - [`octree::Octree`] and [`aabb_tree::AabbTree`] — spatial indexes
//!   used by the cutters.
//! - [`triangle::Triangle`] — flat triangle primitive +
//!   `from_valenx_mesh` adapter.
//!
//! # v1 limitations
//!
//! - Tool model is a cylinder with no shank / spindle clearance.
//! - WaterlinePathPlanner emits unordered intersection points (no
//!   loop-closing); downstream CAM passes need to stitch.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aabb_tree;
pub mod cutter;
pub mod error;
pub mod octree;
pub mod panel;
pub mod triangle;

pub use aabb_tree::{AabbTree, Node};
pub use cutter::{
    AdaptiveDropCutter, DropCutter, EdgeCutter, PushCutter, Tool, WaterlinePathPlanner,
};
pub use error::{ErrorCategory, OpencamlibError};
pub use octree::{OctNode, Octree};
pub use panel::OpenCamLibPanelState;
pub use triangle::{from_valenx_mesh, Triangle};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn flat_triangle_at_z(z: f64) -> Vec<Triangle> {
        vec![
            Triangle::new(
                Vector3::new(0.0, 0.0, z),
                Vector3::new(10.0, 0.0, z),
                Vector3::new(0.0, 10.0, z),
            ),
            Triangle::new(
                Vector3::new(10.0, 0.0, z),
                Vector3::new(10.0, 10.0, z),
                Vector3::new(0.0, 10.0, z),
            ),
        ]
    }

    #[test]
    fn drop_cutter_on_flat_surface() {
        let tris = flat_triangle_at_z(5.0);
        let dc = DropCutter::new(&tris);
        let z = dc.drop(Tool::new(1.0, 5.0), (3.0, 3.0));
        assert!((z - 5.0).abs() < 1e-9, "z = {z}");
    }

    #[test]
    fn drop_cutter_misses_outside_surface() {
        let tris = flat_triangle_at_z(5.0);
        let dc = DropCutter::new(&tris);
        let z = dc.drop(Tool::new(1.0, 5.0), (20.0, 20.0));
        assert_eq!(z, f64::NEG_INFINITY);
    }

    #[test]
    fn aabb_tree_xy_query_returns_overlapping_tris() {
        let tris = flat_triangle_at_z(0.0);
        let tree = AabbTree::new(&tris);
        let hits = tree.xy_query(5.0, 5.0);
        assert!(!hits.is_empty(), "no triangles found at centre");
    }

    #[test]
    fn octree_xy_query_returns_overlapping_tris() {
        let tris = flat_triangle_at_z(0.0);
        let tree = Octree::new(&tris);
        let hits = tree.xy_query(5.0, 5.0);
        assert!(!hits.is_empty(), "octree returned no candidates");
    }

    #[test]
    fn adaptive_drop_grid_samples_full_region() {
        let tris = flat_triangle_at_z(2.0);
        let adc = AdaptiveDropCutter::new(&tris);
        let grid = adc.drop_grid(Tool::new(0.5, 5.0), ((0.0, 10.0), (0.0, 10.0)), 5.0);
        assert_eq!(grid.len(), 9);
        for (_, z) in grid {
            assert!((z - 2.0).abs() < 1e-9, "z={z}");
        }
    }

    #[test]
    fn waterline_intersects_inclined_triangle() {
        let tri = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(10.0, 0.0, 10.0),
            Vector3::new(0.0, 10.0, 0.0),
        );
        let planner = WaterlinePathPlanner::new(std::slice::from_ref(&tri));
        let pts = planner.waterline_at(5.0);
        assert!(!pts.is_empty(), "no waterline points at Z=5");
    }

    #[test]
    fn push_cutter_emits_trail() {
        let tris = flat_triangle_at_z(3.0);
        let pc = PushCutter::new(&tris);
        let trail = pc.push(Tool::new(0.5, 5.0), (0.5, 5.0), (1.0, 0.0), 8.0, 1.0);
        assert_eq!(trail.len(), 9);
        for p in trail {
            assert!((p.z - 3.0).abs() < 1e-9, "z={}", p.z);
        }
    }

    #[test]
    fn edge_cutter_finds_edge() {
        let tri = Triangle::new(
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(10.0, 0.0, 1.0),
            Vector3::new(5.0, 10.0, 1.0),
        );
        let ec = EdgeCutter::new(std::slice::from_ref(&tri));
        let z = ec.edge_only((5.0, 0.0), 0.1);
        assert!((z - 1.0).abs() < 1e-9);
        let z2 = ec.edge_only((-5.0, -5.0), 0.1);
        assert_eq!(z2, f64::NEG_INFINITY);
    }
}
