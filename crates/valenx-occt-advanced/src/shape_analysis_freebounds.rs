//! Phase 139 — `ShapeAnalysis_FreeBounds` — extract free (non-shared)
//! boundary edges of a shape.
//!
//! ## What OCCT does
//!
//! `ShapeAnalysis_FreeBounds(shape, tolerance)` walks the shape's
//! topology counting how many faces share each edge:
//!
//! - **Shared edge** — an edge incident to exactly two faces (the
//!   normal case in a closed manifold solid).
//! - **Free edge** — an edge incident to only one face (boundary of
//!   an open shell or a sheet body).
//! - **Singular** — an edge incident to ≥3 faces (non-manifold
//!   geometry, typically a defect).
//!
//! Returns two compound shapes: the "free wires" (open wires built
//! from free edges in the right order) and the "closed wires" (closed
//! free-edge loops, which mean the shell has a hole). Critical input
//! to [`crate::shape_upgrade_close_open_wires()`] and to
//! [`crate::shape_analysis_orientedclosedsolid()`] for the closed-shell
//! check.
//!
//! ## v1 status
//!
//! **Honest v1.** Walks the solid's BRep topology via truck-modeling's
//! `boundaries()` + `edge_iter()` and counts face references per edge
//! ID. Mesh-backed solids return [`OcctAdvancedError::Backend`]
//! ("topology unavailable on mesh-backed solid").

use std::collections::HashMap;

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Summary returned by [`shape_analysis_freebounds`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FreeBounds {
    /// Number of distinct edges across the solid.
    pub total_edges: usize,
    /// Edges incident to exactly 2 faces (normal manifold case).
    pub shared_edges: usize,
    /// Edges incident to exactly 1 face (open boundary).
    pub free_edges: usize,
    /// Edges incident to ≥3 faces (non-manifold defect).
    pub singular_edges: usize,
}

/// Walk `solid`'s topology and report the free-boundary breakdown.
///
/// # Errors
///
/// - [`OcctAdvancedError::Backend`] for mesh-backed solids (no BRep
///   topology to walk).
pub fn shape_analysis_freebounds(solid: &Solid) -> Result<FreeBounds, OcctAdvancedError> {
    let brep = match solid {
        Solid::Brep(b) => b,
        Solid::Mesh(_) => {
            return Err(OcctAdvancedError::Backend(
                "shape_analysis_freebounds: mesh-backed solid has no BRep topology".into(),
            ))
        }
    };

    // Count face-incidence per edge ID. Walk every face's boundary
    // edges and tally how many faces each ID appears in.
    let mut counts: HashMap<_, usize> = HashMap::new();
    for shell in brep.boundaries() {
        for face in shell.iter() {
            // Each face's boundary is a list of wires; each wire is a
            // list of edges. Walk both layers to find every edge this
            // face touches.
            for wire in face.boundaries() {
                for edge in wire.iter() {
                    *counts.entry(edge.id()).or_insert(0) += 1;
                }
            }
        }
    }

    let mut shared = 0;
    let mut free = 0;
    let mut singular = 0;
    for &n in counts.values() {
        match n {
            0 => unreachable!("HashMap entry exists ⇒ count ≥ 1"),
            1 => free += 1,
            2 => shared += 1,
            _ => singular += 1,
        }
    }

    Ok(FreeBounds {
        total_edges: counts.len(),
        shared_edges: shared,
        free_edges: free,
        singular_edges: singular,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn cube_has_no_free_edges() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let fb = shape_analysis_freebounds(&cube).unwrap();
        // A closed cube: every edge is shared by exactly 2 faces.
        assert_eq!(fb.total_edges, 12, "cube has 12 distinct edges");
        assert_eq!(fb.shared_edges, 12, "all cube edges shared by 2 faces");
        assert_eq!(fb.free_edges, 0);
        assert_eq!(fb.singular_edges, 0);
    }

    #[test]
    fn mesh_backed_solid_rejected() {
        let mesh = valenx_mesh::Mesh::new("test");
        let s = Solid::from_mesh(mesh);
        let err = shape_analysis_freebounds(&s).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.backend");
    }

    #[test]
    fn freebounds_struct_equality() {
        // Sanity check that FreeBounds is Eq for downstream tests.
        let a = FreeBounds {
            total_edges: 12,
            shared_edges: 12,
            free_edges: 0,
            singular_edges: 0,
        };
        let b = a.clone();
        assert_eq!(a, b);
    }
}
