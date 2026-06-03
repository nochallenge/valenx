//! Mesh-domain boolean modifier — Phase 65 / 65.5.
//!
//! Blender's Boolean modifier carves one mesh by another. The v1 here
//! used to concatenate the operands; **65.5 wires the real
//! co-refinement CSG** from [`valenx_cgal_port::mesh_boolean`].
//!
//! Pipeline:
//!
//! 1. Fan-triangulate every polygon of both operand meshes into a
//!    triangle soup ([`valenx_cgal_port::Mesh3`]).
//! 2. Run the real boolean — co-refines the two soups along their
//!    intersection polylines, classifies each facet inside / outside
//!    the other solid, and selects per the requested op.
//! 3. Weld the result back into an indexed polygon [`Mesh`].
//!
//! The precision limits of the underlying kernel apply here too:
//! float arithmetic (robust for well-conditioned closed-manifold
//! inputs, not an exact kernel) and partial coplanar overlap. See the
//! `valenx-cgal-port` module docs for the full caveat list.

use nalgebra::Vector3;
use valenx_cgal_port::mesh_boolean::{self, Mesh3};

use crate::mesh::Mesh;

/// Boolean op kind.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Op {
    /// A union B.
    Union,
    /// A difference B (a - b).
    Difference,
    /// A intersection B.
    Intersection,
}

/// Combine two meshes per `op` with the real co-refinement CSG.
///
/// Both operands should be closed manifolds — the inside/outside
/// classification step is undefined for meshes with boundary. The
/// result is a fresh polygon [`Mesh`] of triangles (the CSG works in
/// the triangle domain; n-gon faces are fan-triangulated on input and
/// not re-merged on output).
pub fn boolean(a: &Mesh, b: &Mesh, op: Op) -> Mesh {
    let soup_a = mesh_to_soup(a);
    let soup_b = mesh_to_soup(b);
    let result = match op {
        Op::Union => mesh_boolean::union(&soup_a, &soup_b),
        Op::Difference => mesh_boolean::difference(&soup_a, &soup_b),
        Op::Intersection => mesh_boolean::intersection(&soup_a, &soup_b),
    };
    soup_to_mesh(&result)
}

/// Fan-triangulate a polygon mesh into a [`Mesh3`] triangle soup.
fn mesh_to_soup(m: &Mesh) -> Mesh3 {
    let mut soup = Mesh3::new();
    for face in &m.faces {
        if face.len() < 3 {
            continue;
        }
        // Fan around the first vertex.
        for k in 1..face.len() - 1 {
            let (i0, i1, i2) = (face[0], face[k], face[k + 1]);
            if i0 >= m.vertices.len() || i1 >= m.vertices.len() || i2 >= m.vertices.len() {
                continue;
            }
            soup.triangles.push(valenx_cgal_port::Triangle3 {
                v: [m.vertices[i0], m.vertices[i1], m.vertices[i2]],
            });
        }
    }
    soup
}

/// Weld a triangle soup back into an indexed polygon [`Mesh`].
fn soup_to_mesh(soup: &Mesh3) -> Mesh {
    // Weld within a tolerance scaled off the geometry's extent so the
    // dedup behaves consistently regardless of model size.
    let tol = weld_tolerance(soup);
    let (verts, faces) = mesh_boolean::mesh3_to_indexed(soup, tol);
    Mesh {
        vertices: verts,
        faces: faces.iter().map(|f| vec![f[0], f[1], f[2]]).collect(),
    }
}

/// Pick a vertex-weld tolerance: 1e-6 of the bounding-box diagonal,
/// clamped to a sane floor.
fn weld_tolerance(soup: &Mesh3) -> f64 {
    if soup.is_empty() {
        return 1e-9;
    }
    let mut lo = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for t in &soup.triangles {
        for p in &t.v {
            lo = Vector3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vector3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    let diag = (hi - lo).norm();
    (diag * 1e-6).max(1e-9)
}

/// Convenience wrapper — union.
pub fn union(a: &Mesh, b: &Mesh) -> Mesh {
    boolean(a, b, Op::Union)
}

/// Convenience wrapper — difference.
pub fn diff(a: &Mesh, b: &Mesh) -> Mesh {
    boolean(a, b, Op::Difference)
}

/// Convenience wrapper — intersection.
pub fn intersect(a: &Mesh, b: &Mesh) -> Mesh {
    boolean(a, b, Op::Intersection)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Axis-aligned box [lo, hi] as a 6-quad polygon mesh.
    fn box_mesh(lo: Vector3<f64>, hi: Vector3<f64>) -> Mesh {
        let v = vec![
            Vector3::new(lo.x, lo.y, lo.z),
            Vector3::new(hi.x, lo.y, lo.z),
            Vector3::new(hi.x, hi.y, lo.z),
            Vector3::new(lo.x, hi.y, lo.z),
            Vector3::new(lo.x, lo.y, hi.z),
            Vector3::new(hi.x, lo.y, hi.z),
            Vector3::new(hi.x, hi.y, hi.z),
            Vector3::new(lo.x, hi.y, hi.z),
        ];
        let faces = vec![
            vec![0, 3, 2, 1], // bottom
            vec![4, 5, 6, 7], // top
            vec![0, 1, 5, 4], // front
            vec![1, 2, 6, 5], // right
            vec![2, 3, 7, 6], // back
            vec![3, 0, 4, 7], // left
        ];
        Mesh { vertices: v, faces }
    }

    #[test]
    fn disjoint_union_keeps_both_solids() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = box_mesh(Vector3::new(5.0, 5.0, 5.0), Vector3::new(6.0, 6.0, 6.0));
        let u = union(&a, &b);
        // Two welded cubes → 16 unique vertices, 24 triangles.
        assert_eq!(u.n_verts(), 16);
        assert_eq!(u.n_faces(), 24);
    }

    #[test]
    fn disjoint_intersection_is_empty() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = box_mesh(Vector3::new(5.0, 5.0, 5.0), Vector3::new(6.0, 6.0, 6.0));
        let i = intersect(&a, &b);
        assert_eq!(i.n_faces(), 0, "disjoint solids share no volume");
    }

    #[test]
    fn intersection_lens_lies_in_shared_box() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        let b = box_mesh(Vector3::new(1.0, 1.0, 1.0), Vector3::new(3.0, 3.0, 3.0));
        let i = intersect(&a, &b);
        assert!(i.n_faces() > 0, "overlapping boxes share a unit cube");
        for p in &i.vertices {
            assert!(p.x >= 1.0 - 1e-6 && p.x <= 2.0 + 1e-6, "x={}", p.x);
            assert!(p.y >= 1.0 - 1e-6 && p.y <= 2.0 + 1e-6, "y={}", p.y);
            assert!(p.z >= 1.0 - 1e-6 && p.z <= 2.0 + 1e-6, "z={}", p.z);
        }
    }

    #[test]
    fn difference_keeps_a_outside_b() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(4.0, 4.0, 4.0));
        let b = box_mesh(Vector3::new(1.0, 1.0, -1.0), Vector3::new(2.0, 2.0, 5.0));
        let d = diff(&a, &b);
        // The carved result is non-empty and bounded by A's extent.
        assert!(d.n_faces() > 0);
        for p in &d.vertices {
            assert!(p.x >= -1e-6 && p.x <= 4.0 + 1e-6, "x={}", p.x);
        }
    }
}
