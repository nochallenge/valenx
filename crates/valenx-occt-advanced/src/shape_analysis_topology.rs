//! Phase 145 — `BRepCheck_Analyzer` full topology consistency check.
//!
//! ## What OCCT does
//!
//! `BRepCheck_Analyzer(shape, geom_controls)` runs an Eulerian-style
//! topology validation:
//!
//! 1. **Vertex consistency** — each vertex's coordinate must agree
//!    across all edges that reference it.
//! 2. **Edge consistency** — each edge must reference exactly two
//!    vertices (front/back) and at most one curve.
//! 3. **Face consistency** — each face must reference a surface plus
//!    a list of wires bounding it.
//! 4. **Shell consistency** — each shell must reference a non-empty
//!    list of faces forming a 2-manifold.
//! 5. **Solid consistency** — each solid must reference one outer
//!    shell plus optionally any number of inner ("void") shells.
//!
//! Reports the first failure as a `BRepCheck_Status` enum value;
//! the FreeCAD / Salome callers iterate over the analyzer to find
//! every defect.
//!
//! ## v1 status
//!
//! **Honest v1.** Walks the solid's BRep via truck-modeling
//! iterators, tallies counts at each topology level, and reports a
//! `TopologyReport`. The deep "each vertex coordinate agrees across
//! all referring edges" check is Phase 145.5 — depends on a vertex-to-
//! edge incidence index that truck doesn't expose directly.

use std::collections::HashSet;

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Summary returned by [`shape_analysis_topology`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TopologyReport {
    /// Total shells in the solid (a closed solid has 1 outer + N
    /// inner void shells; most solids have just 1).
    pub shells: usize,
    /// Total faces across all shells.
    pub faces: usize,
    /// Distinct edges (de-duplicated by edge ID — shared edges count
    /// once).
    pub edges: usize,
    /// Distinct vertices (de-duplicated by vertex ID).
    pub vertices: usize,
    /// Euler characteristic V − E + F. For a closed orientable
    /// 2-manifold, χ = 2 - 2g where g is the genus. χ = 2 means a
    /// sphere-equivalent (cube, sphere, convex polyhedron); χ = 0
    /// means a torus; χ < 0 means higher genus.
    pub euler_characteristic: i64,
}

/// Run topology consistency analysis on `solid`.
///
/// # Errors
///
/// - [`OcctAdvancedError::Backend`] for mesh-backed solids.
pub fn shape_analysis_topology(solid: &Solid) -> Result<TopologyReport, OcctAdvancedError> {
    let brep = match solid {
        Solid::Brep(b) => b,
        Solid::Mesh(_) => {
            return Err(OcctAdvancedError::Backend(
                "shape_analysis_topology: mesh-backed solid has no BRep topology".into(),
            ))
        }
    };

    // Walk shells / faces directly.
    let shells = brep.boundaries().len();
    let faces: usize = brep.boundaries().iter().map(|s| s.len()).sum();

    // De-duplicate edges + vertices by ID.
    let mut e_seen = HashSet::new();
    for e in brep.edge_iter() {
        e_seen.insert(e.id());
    }
    let mut v_seen = HashSet::new();
    for v in brep.vertex_iter() {
        v_seen.insert(v.id());
    }
    let edges = e_seen.len();
    let vertices = v_seen.len();

    let chi = vertices as i64 - edges as i64 + faces as i64;

    Ok(TopologyReport {
        shells,
        faces,
        edges,
        vertices,
        euler_characteristic: chi,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn cube_has_euler_two() {
        // Classic V - E + F: cube has 8 - 12 + 6 = 2.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let r = shape_analysis_topology(&cube).unwrap();
        assert_eq!(r.shells, 1);
        assert_eq!(r.faces, 6);
        assert_eq!(r.edges, 12);
        assert_eq!(r.vertices, 8);
        assert_eq!(r.euler_characteristic, 2);
    }

    #[test]
    fn mesh_backed_solid_rejected() {
        let mesh = valenx_mesh::Mesh::new("t");
        let s = Solid::from_mesh(mesh);
        let err = shape_analysis_topology(&s).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.backend");
    }
}
