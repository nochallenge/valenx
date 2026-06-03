//! # valenx-salome-bridge
//!
//! Salome-style integrated CAD + meshing + FEM platform bridge —
//! Phase 63. Salome unifies modeling, meshing, and analysis in a
//! single environment; this crate provides the **Study** DAG and
//! the three module facades (geom / mesh / analysis) that the rest
//! of the workspace plugs into.
//!
//! Modules:
//! - [`study::Study`] / [`study::StudyNode`] / [`study::NodeKind`]
//!   — central data model with add/depends/rebuild ops.
//! - [`mesh_module::mesh_solid`] — bridge planner for the meshing
//!   adapters.
//! - [`geom_module::plan`] / [`geom_module::GeomOp`] — bridge
//!   planner for the CAD ops.
//! - [`analysis_module::plan`] / [`analysis_module::Analysis`] —
//!   bridge planner for the FEM solvers.
//! - [`panel::SalomePanelState`] — UI envelope (Study browser).
//! - [`error`] — typed [`SalomeError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod analysis_module;
pub mod error;
pub mod geom_module;
pub mod mesh_module;
pub mod panel;
pub mod study;

pub use analysis_module::{Analysis, Solver};
pub use error::{ErrorCategory, SalomeError};
pub use geom_module::GeomOp;
pub use mesh_module::MeshParams;
pub use panel::SalomePanelState;
pub use study::{NodeId, NodeKind, Study, StudyNode};

#[cfg(test)]
mod tests {
    use super::*;

    fn build_study() -> (Study, NodeId, NodeId, NodeId, NodeId) {
        let mut s = Study::new();
        let g = s.add_solid("box1", vec![]);
        let m = s.add_mesh("mesh1", vec![g]);
        let a = s.add_fem("static1", vec![m]);
        let r = s.add_result("disp1", vec![a]);
        (s, g, m, a, r)
    }

    #[test]
    fn add_returns_monotonic_ids() {
        let (s, g, m, a, r) = build_study();
        assert!(g < m && m < a && a < r);
        assert_eq!(s.nodes.len(), 4);
    }

    #[test]
    fn depends_returns_direct_parents() {
        let (s, g, m, a, _) = build_study();
        assert_eq!(s.depends(g).unwrap(), vec![]);
        assert_eq!(s.depends(m).unwrap(), vec![g]);
        assert_eq!(s.depends(a).unwrap(), vec![m]);
    }

    #[test]
    fn unknown_node_errors() {
        let (s, _, _, _, _) = build_study();
        assert!(matches!(s.depends(999), Err(SalomeError::UnknownNode(999))));
    }

    #[test]
    fn rebuild_runs_topo_order() {
        let (mut s, g, m, a, r) = build_study();
        let order = s.rebuild(g).unwrap();
        assert_eq!(order.first(), Some(&g));
        assert_eq!(order.last(), Some(&r));
        // m before a, a before r.
        let pm = order.iter().position(|x| *x == m).unwrap();
        let pa = order.iter().position(|x| *x == a).unwrap();
        let pr = order.iter().position(|x| *x == r).unwrap();
        assert!(pm < pa && pa < pr);
    }

    #[test]
    fn rebuild_marks_up_to_date() {
        let (mut s, g, _, _, _) = build_study();
        s.rebuild(g).unwrap();
        for n in &s.nodes {
            assert!(n.up_to_date);
        }
    }

    #[test]
    fn mesh_module_rejects_bad_h() {
        let r = mesh_module::mesh_solid("s", &MeshParams { h: -1.0, order: 1 });
        assert!(matches!(r, Err(SalomeError::BadParameter { .. })));
    }

    #[test]
    fn mesh_module_rejects_bad_order() {
        let r = mesh_module::mesh_solid("s", &MeshParams { h: 1.0, order: 3 });
        assert!(matches!(r, Err(SalomeError::BadParameter { .. })));
    }

    #[test]
    fn geom_module_box_rejects_zero() {
        let r = geom_module::plan(
            &GeomOp::Box {
                lx: 1.0,
                ly: 0.0,
                lz: 1.0,
            },
            "b",
        );
        assert!(matches!(r, Err(SalomeError::BadParameter { .. })));
    }

    #[test]
    fn analysis_module_rejects_empty_mesh() {
        let r = analysis_module::plan(
            &Analysis {
                mesh_name: String::new(),
                load_case: "lc".into(),
                solver: Solver::LinearElastic,
            },
            "r",
        );
        assert!(matches!(r, Err(SalomeError::BadParameter { .. })));
    }

    #[test]
    fn panel_select_updates_status() {
        let mut p = SalomePanelState::new();
        p.select(Some(7));
        assert_eq!(p.selection, Some(7));
        assert!(p.last_status.is_some());
    }
}
