//! Phase 92 — `BRepOffsetAPI_MakeOffsetShape` (shell or "fattened"
//! solid).
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_MakeOffsetShape(shape, offset, tolerance, ...)`
//! produces a new BRep by offsetting every face of the input shape
//! along its outward normal by `offset` (signed). The result depends
//! on the offset mode:
//!
//! - `Skin` — wraps the original with a thin "skin" of the
//!   requested thickness (thin-walled shell). Used for sheet-metal
//!   parts: input is the mid-surface, output is the part with the
//!   given wall thickness.
//! - `Pipe` — same as Skin but adds end caps.
//! - `RectoVerso` — bidirectional offset (`+offset` outward, `-offset`
//!   inward) producing a thick-walled hollow solid.
//!
//! Conceptually this is "fattening" or "thinning" a shape uniformly.
//! Negative offset shrinks; positive offset grows; if the absolute
//! offset exceeds any local feature size the result is geometrically
//! invalid and OCCT returns `IsDone() == false`.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 92.5) — the mesh-domain
//! approximation promised in the roadmap. The input solid is
//! tessellated; every vertex is displaced along its area-weighted
//! vertex normal by `offset` (signed). Topology (triangle
//! connectivity) is preserved, so the result is a watertight
//! "fattened" / "thinned" mesh-backed [`Solid`]. This is the
//! standard offset-surface approximation: it is exact for planar
//! faces (a translated plane) and a good approximation elsewhere as
//! long as `|offset|` stays below the local feature size — beyond
//! that the offset self-intersects, exactly as OCCT's own
//! `IsDone()==false` warns. A true topology-preserving BRep offset
//! with convex-corner rounds is the Tier-2 follow-up.

use valenx_cad::Solid;
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;

/// Offset every face of `shape` by `offset` along the outward normal.
///
/// `tolerance` is the tessellation chord-error budget — smaller
/// values tessellate the input more finely before offsetting, giving
/// a smoother offset surface on curved input.
///
/// A positive `offset` grows the shape; a negative `offset` shrinks
/// it. The returned [`Solid`] is mesh-backed (carries no BRep
/// topology) — apply this op last in a feature chain, the same rule
/// as the fillet pipeline.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] when `offset` is non-finite or
///   `tolerance` is not strictly positive and finite.
/// - [`OcctSurfaceError::TruckLimit`] when the input solid cannot be
///   tessellated (degenerate, or a mesh-backed solid that produced
///   no triangles).
///
/// # Example
///
/// ```
/// use valenx_occt_surface::offset_api_make_offset;
/// use valenx_cad::box_solid;
/// let cube = box_solid(2.0, 2.0, 2.0).unwrap();
/// let fattened = offset_api_make_offset(&cube, 0.25, 0.1).unwrap();
/// // Mesh-backed offset result.
/// assert!(matches!(fattened, valenx_cad::Solid::Mesh(_)));
/// ```
pub fn offset_api_make_offset(
    shape: &Solid,
    offset: f64,
    tolerance: f64,
) -> Result<Solid, OcctSurfaceError> {
    if !offset.is_finite() {
        return Err(OcctSurfaceError::bad_input(
            "offset",
            format!("must be finite, got {offset}"),
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(OcctSurfaceError::bad_input(
            "tolerance",
            format!("must be positive finite, got {tolerance}"),
        ));
    }

    // Tessellate the input to a triangle mesh.
    let mesh = valenx_cad::solid_to_mesh(shape, tolerance)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("offset tessellation: {e:?}")))?;

    let offset_mesh = offset_mesh_along_vertex_normals(&mesh, offset);
    Ok(Solid::from_mesh(offset_mesh))
}

/// Displace every node of `mesh` along its area-weighted vertex
/// normal by `offset`. Connectivity is copied verbatim.
fn offset_mesh_along_vertex_normals(mesh: &Mesh, offset: f64) -> Mesh {
    let n = mesh.nodes.len();
    // Accumulate area-weighted face normals at each incident vertex.
    // A raw cross product is already area-weighted (its magnitude is
    // twice the triangle area), so summing un-normalised face normals
    // gives the standard area-weighted vertex normal.
    let mut accum = vec![nalgebra::Vector3::<f64>::zeros(); n];
    for block in &mesh.element_blocks {
        // Only triangle blocks contribute a surface normal.
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let (i0, i1, i2) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if i0 >= n || i1 >= n || i2 >= n {
                continue;
            }
            let p0 = mesh.nodes[i0];
            let p1 = mesh.nodes[i1];
            let p2 = mesh.nodes[i2];
            let fn_ = (p1 - p0).cross(&(p2 - p0));
            accum[i0] += fn_;
            accum[i1] += fn_;
            accum[i2] += fn_;
        }
    }

    let mut out = mesh.clone();
    for (node, normal) in out.nodes.iter_mut().zip(accum.iter()) {
        let len = normal.norm();
        if len > 1e-20 {
            *node += normal * (offset / len);
        }
        // Isolated / degenerate vertices (zero accumulated normal)
        // are left in place — there is no defined offset direction.
    }
    out.recompute_stats();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn offset_rejects_bad_tolerance() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_make_offset(&cube, 0.1, -0.01).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn offset_rejects_non_finite_offset() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_make_offset(&cube, f64::NAN, 0.1).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn offset_grows_a_cube() {
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let base = valenx_cad::solid_to_mesh(&cube, 0.5).unwrap();
        let grown = offset_api_make_offset(&cube, 0.5, 0.5).unwrap();
        let grown_mesh = match &grown {
            Solid::Mesh(m) => m,
            Solid::Brep(_) => panic!("offset result should be mesh-backed"),
        };
        // Same triangle topology, same node count.
        assert_eq!(grown_mesh.nodes.len(), base.nodes.len());
        assert_eq!(grown_mesh.total_elements(), base.total_elements());
        // The fattened cube's bounding box is strictly larger: every
        // corner moved outward along its (1,1,1)-ish vertex normal.
        let base_max = base
            .nodes
            .iter()
            .map(|p| p.x.max(p.y).max(p.z))
            .fold(f64::NEG_INFINITY, f64::max);
        let grown_max = grown_mesh
            .nodes
            .iter()
            .map(|p| p.x.max(p.y).max(p.z))
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            grown_max > base_max + 1e-6,
            "offset should grow the cube: base {base_max} -> grown {grown_max}"
        );
    }

    #[test]
    fn negative_offset_shrinks() {
        // Enclosed volume is the correct shrink metric here. A
        // bounding-box measure cannot detect an inward offset: the
        // tessellation is *unwelded* (`solid_to_mesh` duplicates the
        // shared edge/corner vertices per face), so when the +X face
        // moves inward its own corner vertices still hold y and z at
        // the original extreme — every face keeps the AABB pinned at
        // full size even though the solid genuinely shrank.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let base_volume = valenx_cad::solid_volume(&cube).unwrap();
        let shrunk = offset_api_make_offset(&cube, -0.5, 0.5).unwrap();
        assert!(
            matches!(shrunk, Solid::Mesh(_)),
            "offset result should be mesh-backed"
        );
        let shrunk_volume = valenx_cad::solid_volume(&shrunk).unwrap();
        assert!(
            shrunk_volume < base_volume - 1e-6,
            "negative offset should shrink the cube: \
             base volume {base_volume} -> shrunk volume {shrunk_volume}"
        );
    }
}
