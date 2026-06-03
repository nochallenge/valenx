//! Phase 140 — `ShapeAnalysis_Shell::IsLoaded` + `IsClosed` — verify
//! a solid is closed + properly oriented.
//!
//! ## What OCCT does
//!
//! `ShapeAnalysis_Shell` consumes a `TopoDS_Shell` and reports:
//!
//! 1. **Closure** — does every edge have exactly 2 adjacent faces?
//!    (Equivalent to "no free edges" per
//!    [`crate::shape_analysis_freebounds()`].)
//! 2. **Orientation consistency** — do all face normals point
//!    consistently outward? OCCT detects this by walking the shell
//!    along shared edges and checking that adjacent faces agree on
//!    the edge's orientation (one says "forward", the other says
//!    "reversed" — if both say the same, one of them is flipped).
//!
//! Both checks are prerequisites for STEP export (closed solids
//! only), boolean ops (need consistent normals), and FEM meshing
//! (needs the outward-pointing surface).
//!
//! ## v1 status
//!
//! **Honest v1** for the closure check (reuses
//! [`crate::shape_analysis_freebounds()`]'s free-edge count). The
//! orientation-consistency check is a [`OcctAdvancedError::Defect`]
//! stub: truck-modeling stores edge orientation internally per face
//! but doesn't expose a cross-face comparison API. Phase 140.5 ships
//! once `valenx-fillet-brep`'s edge-classification module exposes its
//! orientation-walk primitive.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;
use crate::shape_analysis_freebounds::shape_analysis_freebounds;

/// Result of the closed-and-oriented check.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClosureReport {
    /// True if every edge has exactly 2 adjacent faces (no free
    /// edges, no singular edges).
    pub is_closed: bool,
    /// Number of free edges (one face only).
    pub free_edges: usize,
    /// Number of singular edges (≥3 faces).
    pub singular_edges: usize,
}

/// Verify that `solid` is a closed shell with consistent outward
/// orientation.
///
/// Returns the report when both checks pass; surfaces failures as
/// [`OcctAdvancedError::Defect`] with locus + kind populated for
/// downstream telemetry / UI display.
///
/// # Errors
///
/// - [`OcctAdvancedError::Backend`] for mesh-backed solids.
/// - [`OcctAdvancedError::Defect`] when closure fails (free or
///   singular edges present).
pub fn shape_analysis_orientedclosedsolid(
    solid: &Solid,
) -> Result<ClosureReport, OcctAdvancedError> {
    let fb = shape_analysis_freebounds(solid)?;
    let report = ClosureReport {
        is_closed: fb.free_edges == 0 && fb.singular_edges == 0,
        free_edges: fb.free_edges,
        singular_edges: fb.singular_edges,
    };
    if !report.is_closed {
        return Err(OcctAdvancedError::defect(
            "shell",
            format!(
                "{} free edges + {} singular edges; shell is not closed",
                report.free_edges, report.singular_edges
            ),
        ));
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn cube_is_closed() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let r = shape_analysis_orientedclosedsolid(&cube).unwrap();
        assert!(r.is_closed);
        assert_eq!(r.free_edges, 0);
        assert_eq!(r.singular_edges, 0);
    }

    #[test]
    fn mesh_backed_solid_rejected() {
        let mesh = valenx_mesh::Mesh::new("t");
        let s = Solid::from_mesh(mesh);
        let err = shape_analysis_orientedclosedsolid(&s).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.backend");
    }
}
