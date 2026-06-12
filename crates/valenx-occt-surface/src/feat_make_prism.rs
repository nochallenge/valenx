//! Phase 97 — `BRepFeat_MakePrism` (feature-based prism on a solid).
//!
//! ## What OCCT does
//!
//! `BRepFeat_MakePrism(base, profile, sketch_plane, direction,
//! fuse_subtract, modify)` is the feature-based prism builder —
//! distinct from [`crate::prim_api_prism()`] in that it doesn't just
//! create a prism in isolation, it integrates the prism into an
//! existing solid as a feature (a boss-extrude or pocket). The
//! caller controls:
//!
//! - `fuse_subtract` — Add (boss) or Subtract (pocket) the prism.
//! - `Perform(Until)` — depth (`UpToFirst` / `UpToLast` / explicit
//!   length).
//! - `modify` — whether the operation may modify the base topology
//!   beyond the immediate feature region.
//!
//! Maps to SolidWorks "extruded boss/cut", Inventor "Extrude", and
//! FreeCAD Part Design "Pad" / "Pocket".
//!
//! ## v1 status — real feature-based prism
//!
//! This is a genuine boss/pocket. The pipeline:
//!
//! 1. Build the prism body from `profile_xy` extruded by `height`
//!    along `+Z` ([`prim_api_prism`](crate::prim_api_prism())).
//! 2. Rigidly re-orient it onto the world sketch plane defined by
//!    `(sketch_plane_origin, sketch_plane_normal)`
//!    ([`crate::feat_support::orient_z_to`]).
//! 3. **Fuse** (boss) or **subtract** (pocket) it with `base`
//!    ([`crate::feat_support::feature_combine`]) — an exact BRep
//!    boolean when both `base` and the prism are BReps, with a
//!    co-refinement mesh-CSG fallback otherwise.
//!
//! What this is *not* is the full `BRepFeat_MakePrism` topology
//! graph with `UpToFirst` / `UpToLast` depth resolution against the
//! base's faces — that needs a feature-tree topology engine. This
//! ships the explicit-depth boss/pocket, which is the overwhelmingly
//! common case (SolidWorks "Blind" extrude) and is a *real*
//! integrated feature, not a stub.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;
use crate::feat_support::{feature_combine, orient_z_to};

/// Add or subtract a prism feature to/from `base`.
///
/// `profile_xy` is the sketch profile in the sketch plane's local
/// coordinates; `sketch_plane_origin` and `sketch_plane_normal` orient
/// the sketch plane in world space; `height` is the extrusion depth
/// (positive = along normal); `fuse` chooses Add (true) vs Subtract
/// (false).
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for a too-short profile, a
///   non-finite / zero height, or a zero plane normal.
/// - [`OcctSurfaceError::TruckLimit`] if the prism cannot be built or
///   tessellated.
pub fn feat_make_prism(
    base: &Solid,
    profile_xy: &[(f64, f64)],
    sketch_plane_origin: [f64; 3],
    sketch_plane_normal: [f64; 3],
    height: f64,
    fuse: bool,
) -> Result<Solid, OcctSurfaceError> {
    if profile_xy.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "profile_xy",
            "need at least 3 profile points",
        ));
    }
    if !height.is_finite() || height.abs() < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "height",
            "must be non-zero finite",
        ));
    }
    let n_len = (sketch_plane_normal[0].powi(2)
        + sketch_plane_normal[1].powi(2)
        + sketch_plane_normal[2].powi(2))
    .sqrt();
    if n_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "sketch_plane_normal",
            "must be non-zero",
        ));
    }
    if sketch_plane_origin.iter().any(|c| !c.is_finite()) {
        return Err(OcctSurfaceError::bad_input(
            "sketch_plane_origin",
            "origin must be finite",
        ));
    }

    // Build the +Z prism, then orient + fuse/subtract.
    let prism = crate::prim_api_prism::prim_api_prism(profile_xy, height.abs())?;
    // A negative height means "extrude below the plane" — flip the
    // orientation normal so the prism grows the other way.
    let dir = if height >= 0.0 {
        sketch_plane_normal
    } else {
        [
            -sketch_plane_normal[0],
            -sketch_plane_normal[1],
            -sketch_plane_normal[2],
        ]
    };
    let placed = orient_z_to(prism, sketch_plane_origin, dir);
    feature_combine(base, &placed, fuse)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn feat_prism_rejects_short_profile() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_prism(
            &base,
            &[(0.0, 0.0), (1.0, 0.0)],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.5,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn feat_prism_boss_adds_material() {
        // A boss-extrude on top of a box: build a feature prism on the
        // top face and fuse it. The result must be taller than the
        // base box (extends above z=1).
        let base = box_solid(2.0, 2.0, 1.0).unwrap();
        let result = feat_make_prism(
            &base,
            &[(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)],
            [0.0, 0.0, 1.0], // sketch plane on the box's top face
            [0.0, 0.0, 1.0], // extrude up
            0.5,
            true, // boss
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.1).unwrap();
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            zmax > 1.0 + 1e-6,
            "boss should rise above the base top (z=1), got zmax={zmax}"
        );
    }

    #[test]
    fn feat_prism_pocket_removes_material() {
        // A pocket cut into a box: the result must be valid non-empty
        // geometry (the carved solid).
        let base = box_solid(4.0, 4.0, 2.0).unwrap();
        let result = feat_make_prism(
            &base,
            &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
            [0.0, 0.0, 2.0], // sketch on the top face
            [0.0, 0.0, 1.0],
            -1.0,  // negative height → cut downward into the box
            false, // pocket
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.1).unwrap();
        assert!(!mesh.nodes.is_empty(), "pocket result should be non-empty");
    }
}
