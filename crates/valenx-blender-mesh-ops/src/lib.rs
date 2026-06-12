//! # valenx-blender-mesh-ops
//!
//! Blender precision-modeling subset — Phase 65.
//!
//! Modules:
//! - [`mesh::Mesh`] — polygon mesh primitive.
//! - [`extrude::region`] — extrude a set of faces.
//! - [`bevel::edges`] — bevel one or more edges.
//! - [`inset::faces`] — inset one or more faces.
//! - [`loop_cut::insert`] — insert mid-edge vertices along an edge
//!   loop.
//! - [`bridge::edge_loops`] — connect two same-length edge loops
//!   with a quad strip.
//! - [`bool_modifier::union`] / `diff` / `intersect` — mesh-domain
//!   boolean modifier (real co-refinement CSG via `valenx-cgal-port`;
//!   see module docs for precision limits).
//! - [`solidify::shell`] — Solidify modifier (Blender naming for
//!   shell).
//! - [`panel::BlenderOpPanelState`] — UI envelope with op palette.
//! - [`error`] — typed [`BlenderOpError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bevel;
pub mod bool_modifier;
pub mod bridge;
pub mod error;
pub mod extrude;
pub mod inset;
pub mod loop_cut;
pub mod mesh;
pub mod panel;
pub mod solidify;

pub use error::{BlenderOpError, ErrorCategory};
pub use mesh::Mesh;
pub use panel::{BlenderOp, BlenderOpPanelState};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn extrude_region_top_face_grows_mesh() {
        let m = Mesh::unit_cube();
        let out = extrude::region(&m, &[1], Vector3::new(0.0, 0.0, 0.5)).unwrap();
        // 8 + 4 verts (top face is a quad)
        assert_eq!(out.n_verts(), 12);
        // 6 + 4 wall faces.
        assert_eq!(out.n_faces(), 10);
    }

    #[test]
    fn extrude_region_rejects_empty() {
        let m = Mesh::unit_cube();
        let r = extrude::region(&m, &[], Vector3::z());
        assert!(matches!(r, Err(BlenderOpError::BadParameter { .. })));
    }

    #[test]
    fn bevel_grows_vertex_count() {
        let m = Mesh::unit_cube();
        let r = bevel::edges(&m, &[(0, 1)], 0.1, 3).unwrap();
        // 3+1 new bevel-row vertices per edge.
        assert_eq!(r.n_verts(), 8 + 4);
    }

    #[test]
    fn bevel_rejects_negative_distance() {
        let m = Mesh::unit_cube();
        let r = bevel::edges(&m, &[(0, 1)], -0.1, 2);
        assert!(matches!(r, Err(BlenderOpError::BadParameter { .. })));
    }

    #[test]
    fn bevel_face_less_edge_stays_finite() {
        // An edge whose endpoints belong to no face leaves both endpoint
        // vertex normals at zero; pre-fix the unguarded `.normalize()` of
        // the zero interpolant produced NaN vertices that poison the mesh.
        // The bevel must instead emit finite positions (no offset where
        // the normal is undefined).
        let mut m = Mesh::unit_cube();
        m.faces.clear(); // all 8 vertices now belong to no face
        let (v0, v1) = (m.vertices[0], m.vertices[1]);
        let base = m.vertices.len(); // bevel rows are appended after the originals
        let r = bevel::edges(&m, &[(0, 1)], 0.1, 3).unwrap();
        for v in &r.vertices {
            assert!(
                v.x.is_finite() && v.y.is_finite() && v.z.is_finite(),
                "bevel emitted a non-finite vertex {v:?}"
            );
        }
        // With the degenerate (zero) normal the fallback is no offset, so
        // every bevel-row vertex lies exactly on the edge -- its ends
        // coincide with the edge endpoints. This pins the on-edge fallback
        // semantics against drift (e.g. a future change to a nonzero
        // fallback would push these off the edge).
        assert!(
            (r.vertices[base] - v0).norm() < 1e-12,
            "bevel row should start on the edge"
        );
        assert!(
            (r.vertices[base + 3] - v1).norm() < 1e-12,
            "bevel row should end on the edge"
        );
    }

    #[test]
    fn inset_top_face_doubles_face_count_minus_one() {
        let m = Mesh::unit_cube();
        let r = inset::faces(&m, &[1], 0.1).unwrap();
        // 6 + 4 ring quads = 10. Vertices: 8 + 4 = 12.
        assert_eq!(r.n_faces(), 10);
        assert_eq!(r.n_verts(), 12);
    }

    #[test]
    fn loop_cut_inserts_n_verts_per_edge() {
        let m = Mesh::unit_cube();
        let r = loop_cut::insert(&m, &[(0, 1), (1, 2)], 2).unwrap();
        assert_eq!(r.n_verts(), 8 + 2 * 2);
    }

    #[test]
    fn bridge_loops_rejects_mismatched_lengths() {
        let m = Mesh::unit_cube();
        let r = bridge::edge_loops(&m, &[0, 1, 2], &[4, 5]);
        assert!(matches!(r, Err(BlenderOpError::Topology(_))));
    }

    #[test]
    fn bridge_loops_adds_quads() {
        let m = Mesh::unit_cube();
        let r = bridge::edge_loops(&m, &[0, 1, 2, 3], &[4, 5, 6, 7]).unwrap();
        assert_eq!(r.n_faces(), 6 + 4);
    }

    #[test]
    fn boolean_union_of_disjoint_cubes_keeps_both() {
        // 65.5: real co-refinement CSG. Two disjoint unit cubes union
        // to two welded cubes — 16 verts, 24 triangles.
        let a = Mesh::unit_cube();
        let b = {
            let mut m = Mesh::unit_cube();
            for v in &mut m.vertices {
                v.x += 5.0;
            }
            m
        };
        let r = bool_modifier::union(&a, &b);
        assert_eq!(r.n_verts(), 16);
        assert_eq!(r.n_faces(), 24);
    }

    #[test]
    fn boolean_difference_of_disjoint_cubes_keeps_a() {
        // Subtracting a far-away cube leaves A untouched (24 tris).
        let a = Mesh::unit_cube();
        let b = {
            let mut m = Mesh::unit_cube();
            for v in &mut m.vertices {
                v.x += 5.0;
            }
            m
        };
        let r = bool_modifier::diff(&a, &b);
        assert_eq!(r.n_verts(), 8);
        assert_eq!(r.n_faces(), 12);
    }

    #[test]
    fn solidify_doubles_face_count_closed() {
        let m = Mesh::unit_cube();
        let r = solidify::shell(&m, 0.05).unwrap();
        // Cube is closed -> no boundary stitching; double the faces.
        assert_eq!(r.n_faces(), 12);
        assert_eq!(r.n_verts(), 16);
    }

    #[test]
    fn solidify_rejects_negative_thickness() {
        let m = Mesh::unit_cube();
        let r = solidify::shell(&m, -0.1);
        assert!(matches!(r, Err(BlenderOpError::BadParameter { .. })));
    }

    #[test]
    fn panel_palette_has_seven() {
        assert_eq!(BlenderOp::all().len(), 7);
    }

    #[test]
    fn panel_select_records_status() {
        let mut p = BlenderOpPanelState::new();
        p.select(BlenderOp::Extrude);
        assert_eq!(p.selected, Some(BlenderOp::Extrude));
        assert!(p.last_status.is_some());
    }
}
