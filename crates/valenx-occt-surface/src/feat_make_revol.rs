//! Phase 98 — `BRepFeat_MakeRevol` (feature-based revolution).
//!
//! ## What OCCT does
//!
//! Feature-based counterpart to [`crate::prim_api_revol()`]: revolves a
//! profile around an axis and adds/subtracts the result to/from an
//! existing solid in one call. Maps to SolidWorks "Revolved
//! boss/cut", Inventor "Revolve", and FreeCAD Part Design "Revolution"
//! / "Groove".
//!
//! ## v1 status — real feature-based revolution
//!
//! A genuine revolved boss/groove: the profile is revolved into a
//! solid of revolution by [`prim_api_revol`](crate::prim_api_revol())
//! (already placed in world space against the supplied axis), then
//! **fused** (boss)
//! or **subtracted** (groove) with `base` via
//! [`crate::feat_support::feature_combine`] — exact BRep boolean
//! when possible, co-refinement mesh CSG as the fallback.
//!
//! As with [`crate::feat_make_prism()`], the full `BRepFeat_MakeRevol`
//! topology graph with up-to-face depth resolution is not modelled —
//! this is the explicit-angle revolved feature, the common case.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;
use crate::feat_support::feature_combine;

/// Add or subtract a revolution feature to/from `base`.
///
/// `profile_xy` is the sketch profile; `axis_origin` and `axis_dir`
/// define the rotation axis; `angle_rad` is the sweep amount; `fuse`
/// chooses Add vs Subtract.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for a too-short profile, a zero
///   axis, or a non-finite / zero angle.
/// - [`OcctSurfaceError::TruckLimit`] if the revolution cannot be
///   built or tessellated.
pub fn feat_make_revol(
    base: &Solid,
    profile_xy: &[(f64, f64)],
    axis_origin: [f64; 3],
    axis_dir: [f64; 3],
    angle_rad: f64,
    fuse: bool,
) -> Result<Solid, OcctSurfaceError> {
    if profile_xy.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "profile_xy",
            "need at least 3 profile points",
        ));
    }
    let dir_len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
    if dir_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "axis_dir",
            "must be non-zero",
        ));
    }
    if !angle_rad.is_finite() || angle_rad.abs() < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "angle_rad",
            "must be non-zero finite",
        ));
    }

    // Build the solid of revolution (already world-placed), then fuse.
    let revol = crate::prim_api_revol::prim_api_revol(
        profile_xy,
        axis_origin,
        axis_dir,
        angle_rad,
    )?;
    feature_combine(base, &revol, fuse)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn feat_revol_rejects_zero_axis() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_revol(
            &base,
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
            [0.0; 3],
            [0.0; 3],
            1.0,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn feat_revol_boss_produces_valid_geometry() {
        // Revolve a small rectangle a full turn around the Y axis,
        // fusing the resulting ring onto a base box. The result must
        // be non-empty valid geometry.
        let base = box_solid(6.0, 2.0, 6.0).unwrap();
        let result = feat_make_revol(
            &base,
            &[(2.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0)],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            std::f64::consts::TAU,
            true,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty(), "revolved boss should be non-empty");
    }
}
