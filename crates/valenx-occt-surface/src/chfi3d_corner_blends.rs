//! Phase 74 — `ChFi3d` full corner blends at vertices where 3+
//! filleted edges meet.
//!
//! ## What OCCT does
//!
//! `ChFi3d_FilBuilder` is the OCCT topology-aware fillet builder.
//! Once the caller adds three or more concurrent edges that all carry
//! a fillet (`Add(radius, edge)`), the builder must construct a
//! **corner blend** — an N-sided patch that smoothly merges all
//! incoming fillet ribbons at the shared vertex. The math:
//!
//! 1. Each fillet ribbon is a `Geom_BSplineSurface` (tube wrapping the
//!    edge).
//! 2. At the vertex, the ribbons truncate against each other along
//!    their mutual intersection curves.
//! 3. The "hole" left between the ribbons is filled with a Coons-like
//!    N-sided patch (`GeomFill_ConstrainedFilling`) that interpolates
//!    each ribbon's boundary tangent.
//!
//! Without a corner blend the surface has G0 creases at the vertex,
//! which downstream STEP export will accept but renderers will show
//! as visible facet boundaries.
//!
//! ## v1 status
//!
//! Stub — depends on a true BRep fillet, which `valenx-cad::fillet`
//! itself defers (truck has no edge-fillet operator). Phase 74.5
//! ships once Phase 3.5 lands. Mesh-domain corner fillets are
//! available in `valenx-fillet` for tessellated workflows.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Apply a full corner blend at the vertex shared by `concurrent_edges`.
///
/// `radii_per_edge` carries the fillet radius for each edge in the
/// same order as `concurrent_edges`. Edge identity is given by the
/// edge index in the solid's edge_iter() order — see
/// `valenx-cad::Solid::edges`.
///
/// # Errors
///
/// Always [`OcctSurfaceError::NotYetImplemented`] in v1.
pub fn chfi3d_corner_blends(
    _solid: &Solid,
    concurrent_edges: &[usize],
    radii_per_edge: &[f64],
) -> Result<Solid, OcctSurfaceError> {
    if concurrent_edges.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "concurrent_edges",
            format!(
                "corner blends apply when 3+ filleted edges meet; got {}",
                concurrent_edges.len()
            ),
        ));
    }
    if concurrent_edges.len() != radii_per_edge.len() {
        return Err(OcctSurfaceError::bad_input(
            "radii_per_edge",
            format!(
                "must match concurrent_edges length ({} vs {})",
                radii_per_edge.len(),
                concurrent_edges.len()
            ),
        ));
    }
    for (i, &r) in radii_per_edge.iter().enumerate() {
        if !r.is_finite() || r <= 0.0 {
            return Err(OcctSurfaceError::bad_input(
                "radii_per_edge",
                format!("radii[{i}] must be positive finite, got {r}"),
            ));
        }
    }
    Err(OcctSurfaceError::not_yet("chfi3d_corner_blends"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn corner_blend_rejects_fewer_than_three_edges() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = chfi3d_corner_blends(&cube, &[0, 1], &[0.1, 0.1]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn corner_blend_rejects_length_mismatch() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = chfi3d_corner_blends(&cube, &[0, 1, 2], &[0.1, 0.1]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn corner_blend_rejects_bad_radius() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = chfi3d_corner_blends(&cube, &[0, 1, 2], &[0.1, -0.1, 0.1]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn corner_blend_is_stub_with_valid_inputs() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = chfi3d_corner_blends(&cube, &[0, 1, 2], &[0.1, 0.1, 0.1]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.not_yet_implemented");
    }
}
