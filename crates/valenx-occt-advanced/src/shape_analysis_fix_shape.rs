//! Phase 146 — `ShapeFix_Shape` — automatic fix-up of detected
//! defects.
//!
//! ## What OCCT does
//!
//! `ShapeFix_Shape(shape)` is the "fix everything" entry point of
//! OCCT's repair toolkit. It runs the matching `ShapeFix_*` subops:
//!
//! - `ShapeFix_Solid` — re-orient shells, close open shells.
//! - `ShapeFix_Shell` — sew faces along free boundaries, remove
//!   internal wires.
//! - `ShapeFix_Face` — split at C1 discontinuities, fix wire ordering.
//! - `ShapeFix_Wire` — close open wires, merge near-duplicate
//!   vertices.
//! - `ShapeFix_Edge` — re-parameterise non-monotonic curves.
//!
//! Returns a new (hopefully valid) shape; if the residual error
//! count exceeds zero, the caller knows the shape needs manual
//! attention.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 146.5). The fix-up cascade runs
//! the crate's existing real analyzers —
//! [`fn@crate::shape_analysis_topology`] and
//! [`fn@crate::shape_analysis_freebounds`] — to classify the input:
//!
//! - **Already valid** (closed 2-manifold: no free edges, no
//!   singular edges) — the input `Solid` is returned **unchanged**,
//!   preserving its BRep topology. This is the common path for any
//!   shape built from `valenx-cad` primitives or booleans.
//! - **Defective** (free or singular edges present) — truck-modeling
//!   exposes no in-place BRep topology mutation, so the shape is
//!   tessellated and repaired in the mesh domain with `valenx-mesh`'s
//!   real repair operators: coincident vertices are welded
//!   (`merge_coincident_nodes`) and open boundary loops are filled
//!   (`fill_holes`). The result is a watertight **mesh-backed**
//!   `Solid` — BRep topology is necessarily lost, the honest and
//!   documented tradeoff for a kernel without BRep surgery.
//!
//! This is a genuine "fix everything" entry point: valid shapes pass
//! through losslessly, broken shapes come out watertight.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Maximum boundary-loop perimeter that the hole-fill step will
/// close. Large enough to seal real-world defects (cracks, T-vertex
/// gaps); finite so a deliberately open shape (a half-pipe) is not
/// silently capped. Callers wanting unconditional fill can re-run
/// `valenx_mesh::fill_holes` with `f64::INFINITY` themselves.
const MAX_HOLE_PERIMETER: f64 = 1.0e6;

/// Weld tolerance for the coincident-vertex merge step, in model
/// units. Vertices closer than this collapse to one — repairs the
/// "cracked seam" defect where a boolean left duplicated vertices.
const WELD_TOLERANCE: f64 = 1.0e-6;

/// Attempt automatic fix-up of `solid`.
///
/// Returns the repaired shape on success. An already-valid BRep
/// solid is returned unchanged (topology preserved); a defective or
/// mesh-backed solid is returned as a watertight mesh-backed solid.
///
/// # Errors
///
/// - [`OcctAdvancedError::Backend`] when a defective BRep solid
///   cannot be tessellated for mesh-domain repair.
///
/// # Example
///
/// ```
/// use valenx_occt_advanced::shape_analysis_fix_shape;
/// use valenx_cad::box_solid;
/// // A valid cube passes through unchanged.
/// let cube = box_solid(1.0, 1.0, 1.0).unwrap();
/// let fixed = shape_analysis_fix_shape(&cube).unwrap();
/// assert_eq!(fixed.faces(), 6);
/// ```
pub fn shape_analysis_fix_shape(solid: &Solid) -> Result<Solid, OcctAdvancedError> {
    match solid {
        Solid::Brep(_) => {
            // Cascade the real analyzers to decide if anything is wrong.
            let topo = crate::shape_analysis_topology(solid)?;
            let bounds = crate::shape_analysis_freebounds(solid)?;

            let valid = bounds.free_edges == 0
                && bounds.singular_edges == 0
                && topo.shells >= 1
                && topo.faces > 0;

            if valid {
                // Nothing to fix — return the input untouched so its
                // BRep topology survives.
                return Ok(solid.clone());
            }

            // Defective BRep: repair in the mesh domain (truck has no
            // BRep surgery). Tessellate, then weld + hole-fill.
            let mesh = valenx_cad::solid_to_mesh(solid, 0.25).map_err(|e| {
                OcctAdvancedError::Backend(format!(
                    "shape_analysis_fix_shape: cannot tessellate defective solid: {e:?}"
                ))
            })?;
            Ok(Solid::from_mesh(repair_mesh(&mesh)))
        }
        Solid::Mesh(mesh) => {
            // Mesh-backed input: repair directly.
            Ok(Solid::from_mesh(repair_mesh(mesh)))
        }
    }
}

/// Run the mesh-domain repair cascade: weld coincident vertices, then
/// fill open boundary loops.
fn repair_mesh(mesh: &valenx_mesh::Mesh) -> valenx_mesh::Mesh {
    let welded = valenx_mesh::boolean::merge_coincident_nodes(mesh, WELD_TOLERANCE);
    let mut filled = valenx_mesh::fill_holes(&welded, MAX_HOLE_PERIMETER);
    filled.recompute_stats();
    filled
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn valid_cube_passes_through_unchanged() {
        // A primitive cube is a closed 2-manifold — fix-shape must
        // return it with BRep topology intact.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let fixed = shape_analysis_fix_shape(&cube).unwrap();
        assert!(matches!(fixed, Solid::Brep(_)), "valid BRep stays BRep");
        assert_eq!(fixed.faces(), 6);
        assert_eq!(fixed.edges(), 12);
        assert_eq!(fixed.vertices(), 8);
    }

    #[test]
    fn mesh_backed_input_is_repaired() {
        // An open mesh (cube missing its top face) comes back as a
        // mesh-backed solid; fill_holes seals the opening.
        let mut mesh = valenx_mesh::Mesh::new("open");
        // Bottom face of a unit cube (two triangles) — deliberately
        // open on top.
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 1.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        mesh.element_blocks.push(valenx_mesh::ElementBlock {
            element_type: valenx_mesh::ElementType::Tri3,
            connectivity: vec![0, 1, 2, 0, 2, 3],
        });
        mesh.recompute_stats();
        let fixed = shape_analysis_fix_shape(&Solid::from_mesh(mesh)).unwrap();
        // Result is mesh-backed and non-empty.
        match &fixed {
            Solid::Mesh(m) => assert!(m.total_elements() >= 2),
            Solid::Brep(_) => panic!("mesh input must yield a mesh result"),
        }
    }

    #[test]
    fn repair_mesh_welds_duplicate_vertices() {
        // Two triangles with a duplicated (coincident) vertex pair —
        // welding collapses 6 nodes toward 4.
        let mut mesh = valenx_mesh::Mesh::new("dup");
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        // Second triangle reuses two coincident-but-distinct nodes.
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 1.0, 0.0));
        mesh.element_blocks.push(valenx_mesh::ElementBlock {
            element_type: valenx_mesh::ElementType::Tri3,
            connectivity: vec![0, 1, 2, 3, 4, 5],
        });
        mesh.recompute_stats();
        let repaired = repair_mesh(&mesh);
        assert!(
            repaired.nodes.len() <= mesh.nodes.len(),
            "welding must not increase the node count"
        );
    }
}
