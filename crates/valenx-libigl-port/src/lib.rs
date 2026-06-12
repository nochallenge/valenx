//! # valenx-libigl-port
//!
//! Pure-Rust port of [libigl] discrete-geometry-processing algorithms —
//! Phase 59.
//!
//! [libigl]: https://github.com/libigl/libigl
//!
//! Modules:
//! - [`triangle::TriMesh`] — owned `(vertices, triangles)` pair +
//!   `vertex_one_ring` adjacency.
//! - [`laplacian`] — discrete differential operators (cotangent
//!   Laplacian, lumped mass matrix, dense symmetric solver) — the
//!   substrate the parameterisation / geodesics algorithms rest on.
//! - [`param::lscm`] / [`param::arap`] — UV parameterisation: real
//!   Least-Squares Conformal Mapping and As-Rigid-As-Possible
//!   local/global solve.
//! - [`deform::laplacian_smooth`] — uniform Laplacian smoothing.
//! - [`deform::biharmonic`] — handle-based deformation (v1: BFS-hop
//!   weighted offset; will swap for bi-Laplacian solve).
//! - [`field::heat_geodesics`] — geodesic distance via the real heat
//!   method (Crane et al.): heat flow + Poisson reconstruction.
//! - [`cut::random_cuts`] — deterministic random edge cuts for
//!   unfolding.
//! - [`pca::shape_descriptor`] — eigenvalues of vertex covariance.
//! - [`panel::LibiglPanelState`] — workbench-panel envelope.
//! - [`error`] — typed [`LibiglError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![allow(clippy::needless_range_loop)]

pub mod cut;
pub mod deform;
pub mod error;
pub mod field;
pub mod laplacian;
pub mod panel;
pub mod param;
pub mod pca;
pub mod triangle;

pub use cut::random_cuts;
pub use deform::{biharmonic, laplacian_smooth};
pub use error::{ErrorCategory, LibiglError};
pub use field::heat_geodesics;
pub use panel::LibiglPanelState;
pub use param::{arap, lscm};
pub use pca::shape_descriptor;
pub use triangle::TriMesh;

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn unit_quad_mesh() -> TriMesh {
        TriMesh {
            vertices: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3]],
        }
    }

    #[test]
    fn lscm_returns_one_uv_per_vertex() {
        let m = unit_quad_mesh();
        let uvs = lscm(&m, &[]).unwrap();
        assert_eq!(uvs.len(), m.vertices.len());
    }

    #[test]
    fn arap_returns_one_uv_per_vertex() {
        let m = unit_quad_mesh();
        let uvs = arap(&m).unwrap();
        assert_eq!(uvs.len(), m.vertices.len());
    }

    #[test]
    fn laplacian_smooth_preserves_count() {
        let m = unit_quad_mesh();
        let smoothed = laplacian_smooth(&m, 5, 0.5).unwrap();
        assert_eq!(smoothed.n_verts(), m.n_verts());
        assert_eq!(smoothed.n_tris(), m.n_tris());
    }

    #[test]
    fn laplacian_smooth_bad_lambda() {
        let m = unit_quad_mesh();
        assert!(matches!(
            laplacian_smooth(&m, 1, 2.0),
            Err(LibiglError::BadParameter { name: "lambda", .. })
        ));
    }

    #[test]
    fn biharmonic_pins_handle_exactly() {
        let m = unit_quad_mesh();
        let target = Vector3::new(5.0, 5.0, 5.0);
        let deformed = biharmonic(&m, &[(0, target)]).unwrap();
        let p = deformed.vertices[0];
        assert!((p - target).norm() < 1e-9);
    }

    #[test]
    fn biharmonic_no_handles_errors() {
        let m = unit_quad_mesh();
        assert!(matches!(
            biharmonic(&m, &[]),
            Err(LibiglError::BadParameter {
                name: "handles",
                ..
            })
        ));
    }

    #[test]
    fn heat_geodesics_source_distance_is_zero() {
        let m = unit_quad_mesh();
        let d = heat_geodesics(&m, 0).unwrap();
        assert_eq!(d.len(), m.n_verts());
        assert_eq!(d[0], 0.0);
    }

    #[test]
    fn heat_geodesics_bad_source() {
        let m = unit_quad_mesh();
        assert!(matches!(
            heat_geodesics(&m, 99),
            Err(LibiglError::BadParameter {
                name: "source_vertex",
                ..
            })
        ));
    }

    #[test]
    fn random_cuts_returns_requested_count() {
        let m = unit_quad_mesh();
        let cuts = random_cuts(&m, 3).unwrap();
        assert_eq!(cuts.len(), 3);
    }

    #[test]
    fn random_cuts_deterministic() {
        let m = unit_quad_mesh();
        let a = random_cuts(&m, 3).unwrap();
        let b = random_cuts(&m, 3).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn shape_descriptor_returns_three_eigenvalues() {
        let m = unit_quad_mesh();
        let v = shape_descriptor(&m).unwrap();
        assert_eq!(v.len(), 3);
        assert!(v[0] >= v[1] && v[1] >= v[2]);
    }
}
