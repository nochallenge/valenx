//! valenx-cgal-port — pure-Rust port of CGAL geometry algorithms subset
//! (Phase 58).
//!
//! Modules:
//! - [`delaunay`] — 2D Delaunay (Bowyer-Watson).
//! - [`convex_hull`] — 2D Graham scan + 3D incremental hull.
//! - [`alpha_shape`] — 2D alpha-shape boundary (Delaunay filter + edge stitch).
//! - [`voronoi`] — 2D Voronoi via Delaunay dual.
//! - [`mesh_boolean`] — co-refinement triangle-soup CSG (real
//!   union / difference / intersection; see the module docs for the
//!   float-arithmetic precision limits).
//! - [`aabb_tree`] — 3D AABB tree for ray queries.
//! - [`panel`] — UI panel state envelope.
//! - [`error`] — typed [`CgalError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aabb_tree;
pub mod alpha_shape;
pub mod convex_hull;
pub mod delaunay;
pub mod error;
pub mod mesh_boolean;
pub mod panel;
pub mod voronoi;

pub use aabb_tree::{AabbTree, Node, Triangle3, TriangleId};
pub use alpha_shape::alpha_shape_2d;
pub use convex_hull::{hull_2d, hull_3d};
pub use delaunay::triangulate_2d;
pub use error::{CgalError, ErrorCategory};
pub use mesh_boolean::{
    difference, intersection, mesh3_from_indexed, mesh3_to_indexed, union, Mesh3,
};
pub use panel::CgalPanelState;
pub use voronoi::{voronoi_2d, VoronoiCell};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn hull_2d_of_square_returns_corners() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]];
        let hull = hull_2d(&pts).unwrap();
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn hull_2d_not_enough_points() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0]];
        assert!(matches!(
            hull_2d(&pts),
            Err(CgalError::NotEnoughPoints {
                needed: 3,
                given: 2
            })
        ));
    }

    #[test]
    fn delaunay_triangulates_square() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let tris = triangulate_2d(&pts).unwrap();
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn delaunay_rejects_collinear() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]];
        assert!(matches!(triangulate_2d(&pts), Err(CgalError::Degenerate)));
    }

    #[test]
    fn alpha_shape_returns_boundary() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]];
        let bnd = alpha_shape_2d(&pts, 10.0).unwrap();
        assert!(!bnd.is_empty());
    }

    #[test]
    fn alpha_shape_bad_alpha() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        assert!(matches!(
            alpha_shape_2d(&pts, 0.0),
            Err(CgalError::BadParameter { name: "alpha", .. })
        ));
    }

    #[test]
    fn voronoi_returns_one_cell_per_site() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let cells = voronoi_2d(&pts).unwrap();
        assert_eq!(cells.len(), 3);
        for c in &cells {
            assert!(!c.vertices.is_empty(), "cell at {:?} empty", c.site);
        }
    }

    #[test]
    fn hull_3d_of_tetrahedron() {
        let pts = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let faces = hull_3d(&pts).unwrap();
        assert_eq!(faces.len(), 4);
    }

    #[test]
    fn aabb_tree_returns_ray_hits() {
        let tris = vec![
            Triangle3 {
                v: [
                    Vector3::new(0.0, 0.0, 0.0),
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(0.0, 1.0, 0.0),
                ],
            },
            Triangle3 {
                v: [
                    Vector3::new(2.0, 2.0, 0.0),
                    Vector3::new(3.0, 2.0, 0.0),
                    Vector3::new(2.0, 3.0, 0.0),
                ],
            },
        ];
        let tree = AabbTree::new(tris);
        let hits = tree.intersect_ray(Vector3::new(0.3, 0.3, 5.0), Vector3::new(0.0, 0.0, -1.0));
        assert!(hits.iter().any(|t| t.0 == 0));
    }

    #[test]
    fn mesh_boolean_union_concatenates() {
        let m1 = Mesh3 {
            triangles: vec![Triangle3 {
                v: [
                    Vector3::zeros(),
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(0.0, 1.0, 0.0),
                ],
            }],
        };
        let m2 = Mesh3 {
            triangles: vec![Triangle3 {
                v: [
                    Vector3::new(2.0, 2.0, 0.0),
                    Vector3::new(3.0, 2.0, 0.0),
                    Vector3::new(2.0, 3.0, 0.0),
                ],
            }],
        };
        let u = union(&m1, &m2);
        assert_eq!(u.len(), 2);
    }
}
