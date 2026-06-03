//! Phase 99 — `BRepFeat_MakePipe` (feature-based pipe sweep).
//!
//! ## What OCCT does
//!
//! Feature-based wrapper around [`crate::sweep_api_pipe()`]: sweeps a
//! profile along a spine and adds/subtracts the result to/from an
//! existing solid. Maps to SolidWorks "Swept boss/cut", Inventor
//! "Sweep", and FreeCAD Part Design "Additive Pipe" /
//! "Subtractive Pipe".
//!
//! ## v1 status — real feature-based pipe
//!
//! A genuine swept boss/cut: the profile is swept along the spine by
//! the real [`crate::sweep_api_pipe()`] (parallel-transport frame for
//! curved spines, BRep prism for a straight `+Z` spine), then
//! **fused** (boss) or **subtracted** (cut) with `base` via
//! [`crate::feat_support::feature_combine`] — exact BRep boolean when
//! both operands are BReps, co-refinement mesh CSG otherwise (a
//! curved-spine sweep is mesh-backed, so the pipe-cut path naturally
//! uses the mesh CSG).

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;
use crate::feat_support::feature_combine;

/// Add or subtract a swept feature to/from `base`.
///
/// `profile_xy` is the cross-section; `spine` is the swept path;
/// `fuse` chooses Add vs Subtract.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for a too-short profile or spine.
/// - [`OcctSurfaceError::TruckLimit`] if the sweep cannot be built or
///   tessellated.
pub fn feat_make_pipe(
    base: &Solid,
    profile_xy: &[(f64, f64)],
    spine: &[[f64; 3]],
    fuse: bool,
) -> Result<Solid, OcctSurfaceError> {
    if profile_xy.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "profile_xy",
            "need at least 3 profile points",
        ));
    }
    if spine.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "spine",
            "need at least 2 spine points",
        ));
    }

    // Sweep the profile along the spine, then fuse/subtract.
    let pipe = crate::sweep_api_pipe::sweep_api_pipe(profile_xy, spine)?;
    feature_combine(base, &pipe, fuse)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn feat_pipe_rejects_short_spine() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_pipe(
            &base,
            &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
            &[[0.0; 3]],
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn feat_pipe_boss_produces_valid_geometry() {
        // Sweep a small triangular profile straight up and fuse it
        // onto a base box. The result is non-empty valid geometry.
        let base = box_solid(4.0, 4.0, 1.0).unwrap();
        let result = feat_make_pipe(
            &base,
            &[(1.0, 1.0), (2.0, 1.0), (1.5, 2.0)],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
            true,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty(), "swept boss should be non-empty");
    }

    #[test]
    fn feat_pipe_along_a_curved_spine() {
        // A curved (L-shaped) spine sweep — exercises the mesh-CSG
        // fallback path in feature_combine.
        let base = box_solid(6.0, 6.0, 6.0).unwrap();
        let result = feat_make_pipe(
            &base,
            &[(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0)],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [3.0, 0.0, 3.0]],
            true,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.3).unwrap();
        assert!(!mesh.nodes.is_empty());
    }
}
