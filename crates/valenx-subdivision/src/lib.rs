//! # valenx-subdivision
//!
//! Wings 3D-style subdivision modeling — Phase 61.
//!
//! Modules:
//! - [`mesh::SubdivMesh`] — polygon mesh primitive.
//! - [`catmull_clark::subdivide`] — quad-output scheme for any
//!   polygon mesh.
//! - [`loop_subdiv::subdivide`] — triangle-only refinement scheme.
//! - [`crease::set_crease`] / [`crease::CreasedMesh`] — semi-sharp
//!   crease annotation.
//! - [`catmull_clark_creased::subdivide`] — crease-aware
//!   Catmull-Clark.
//! - [`extrude::extrude_face`] — Wings-style face extrude.
//! - [`panel::SubdivPanelState`] — UI panel envelope.
//! - [`error`] — typed [`SubdivError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod catmull_clark;
pub mod catmull_clark_creased;
pub mod crease;
pub mod error;
pub mod extrude;
pub mod loop_subdiv;
pub mod mesh;
pub mod panel;

pub use crease::{set_crease, CreasedMesh};
pub use error::{ErrorCategory, SubdivError};
pub use extrude::extrude_face;
pub use mesh::SubdivMesh;
pub use panel::{Scheme, SubdivPanelState};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cube_catmull_clark_one_iter_quads_only() {
        let m = SubdivMesh::unit_cube();
        let s = catmull_clark::subdivide(&m, 1);
        // 6 quads → 24 sub-quads.
        assert_eq!(s.n_faces(), 24);
        for f in &s.faces {
            assert_eq!(f.len(), 4);
        }
    }

    #[test]
    fn cube_catmull_clark_two_iter_quad_count() {
        let m = SubdivMesh::unit_cube();
        let s = catmull_clark::subdivide(&m, 2);
        // 24 → 96 quads after second iteration.
        assert_eq!(s.n_faces(), 96);
    }

    #[test]
    fn tetra_loop_one_iter_face_count() {
        let m = SubdivMesh::tetrahedron();
        let s = loop_subdiv::subdivide(&m, 1).expect("loop ok");
        // 4 triangles → 16 triangles.
        assert_eq!(s.n_faces(), 16);
        for f in &s.faces {
            assert_eq!(f.len(), 3);
        }
    }

    #[test]
    fn loop_rejects_non_triangle() {
        let m = SubdivMesh::unit_cube();
        let r = loop_subdiv::subdivide(&m, 1);
        assert!(matches!(r, Err(SubdivError::Topology(_))));
    }

    #[test]
    fn catmull_clark_skips_malformed_mesh_without_panic() {
        // A face index past `vertices` (SubdivMesh has public fields). The
        // infallible Catmull-Clark entry must return rather than panic.
        let mut m = SubdivMesh::unit_cube();
        m.faces.push(vec![0, 1, 99]); // 99 is past the 8 cube vertices
        let _ = catmull_clark::subdivide(&m, 1);
    }

    #[test]
    fn loop_rejects_out_of_range_index() {
        // Loop returns Result, so a face index past `vertices` is a clean error.
        let mut m = SubdivMesh::tetrahedron();
        m.faces.push(vec![0, 1, 99]); // arity 3 passes the arity check; index OOB
        let r = loop_subdiv::subdivide(&m, 1);
        assert!(matches!(r, Err(SubdivError::IndexOutOfRange { .. })));
    }

    #[test]
    fn crease_round_trip() {
        let m = SubdivMesh::unit_cube();
        let c = set_crease(&m, &[(0, 1), (1, 2)], 2.0).expect("ok");
        assert_eq!(c.sharpness(0, 1), 2.0);
        assert_eq!(c.sharpness(1, 0), 2.0); // order-independent
        assert_eq!(c.sharpness(0, 2), 0.0);
    }

    #[test]
    fn crease_rejects_bad_vertex() {
        let m = SubdivMesh::unit_cube();
        let r = set_crease(&m, &[(0, 99)], 1.0);
        assert!(matches!(r, Err(SubdivError::IndexOutOfRange { .. })));
    }

    #[test]
    fn crease_rejects_negative_sharpness() {
        let m = SubdivMesh::unit_cube();
        let r = set_crease(&m, &[(0, 1)], -1.0);
        assert!(matches!(r, Err(SubdivError::BadParameter { .. })));
    }

    #[test]
    fn catmull_clark_creased_keeps_quad_count() {
        let m = SubdivMesh::unit_cube();
        let c = set_crease(&m, &[(0, 1)], 2.0).expect("ok");
        let s = catmull_clark_creased::subdivide(&c, 1);
        assert_eq!(s.n_faces(), 24);
    }

    #[test]
    fn extrude_face_doubles_face_count_plus_walls() {
        let m = SubdivMesh::unit_cube();
        // Cube top face (index 1). 6 faces → 6 + 4 wall quads = 10.
        let s = extrude_face(&m, 1, 0.5).expect("ok");
        assert_eq!(s.n_faces(), 10);
        assert_eq!(s.n_verts(), 8 + 4);
    }

    #[test]
    fn extrude_face_rejects_bad_id() {
        let m = SubdivMesh::unit_cube();
        let r = extrude_face(&m, 99, 0.5);
        assert!(matches!(r, Err(SubdivError::IndexOutOfRange { .. })));
    }

    #[test]
    fn extrude_face_rejects_nan() {
        let m = SubdivMesh::unit_cube();
        let r = extrude_face(&m, 0, f64::NAN);
        assert!(matches!(r, Err(SubdivError::BadParameter { .. })));
    }

    #[test]
    fn panel_default_scheme_catmull_clark() {
        let p = SubdivPanelState::default();
        assert_eq!(p.scheme, Scheme::CatmullClark);
        assert_eq!(p.iterations, 2);
        assert!(p.menu_label.contains("Subdivision"));
    }

    #[test]
    fn panel_set_source_records_status() {
        let mut p = SubdivPanelState::new();
        p.set_source(SubdivMesh::unit_cube());
        assert!(p.last_status.is_some());
        assert!(p.last_error.is_none());
        assert_eq!(p.result.n_verts(), 8);
    }
}
