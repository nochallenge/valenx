//! Phase 137 — `BRepFeat_MakeDPrism` subtractive / pocket variant.
//!
//! ## What OCCT does
//!
//! Counterpart to Phase 136's [`crate::feat_make_dgreater_pad()`] — the
//! "DSubtract" flag flips the prism from additive to subtractive, so
//! the feature carves into `base` rather than growing it. Same
//! topology-aware integration: the pocket can extend past the sketch
//! profile's planar footprint, following the base's curvature so the
//! resulting hole wraps around fillets / curved walls. Used for
//! through-holes on cylindrical shafts, milling pockets on curved
//! aerospace skins, drainage channels on cast parts.
//!
//! ## v1 status — real subtractive pocket
//!
//! A genuine subtractive pocket. The pipeline:
//!
//! 1. Build a prism body of the requested `depth` (which may exceed
//!    the profile's planar footprint — the DSubtract regime).
//! 2. Orient it onto the world sketch plane, extruding *into* the
//!    solid (the prism grows opposite the sketch normal).
//! 3. **Subtract** it from `base` via the real Phase-97
//!    [`valenx_occt_surface::feat_support::feature_combine`] — an
//!    exact BRep boolean when both operands are BReps, with a
//!    co-refinement mesh-CSG fallback otherwise.
//!
//! This delegates to the real Phase-97
//! [`feat_make_prism`](fn@valenx_occt_surface::feat_make_prism) with `fuse = false` and a
//! negative height (so the cutter pokes *into* the base from the
//! sketch plane). As with Phase 136 the `BRepFeat_MakeDPrism`
//! footprint-redirection topology engine (the pocket walls bending
//! along the base's curved faces) stays Tier-3; the boolean-subtracted
//! pocket is real, carved geometry.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Apply a subtractive "DSubtract" pocket to `base`.
///
/// `profile_xy` is the closed sketch profile in the sketch-plane XY;
/// `sketch_plane_origin` / `sketch_plane_normal` orient the plane.
/// `depth` is the pocket depth — positive; the cutter pokes *into* the
/// base (opposite the sketch normal). It may exceed the profile's
/// planar footprint (the DSubtract regime).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for a too-short profile, a
///   non-positive / non-finite depth, or a zero plane normal.
/// - [`OcctAdvancedError::Backend`] if the pocket cannot be built.
pub fn feat_make_dsubtract_pocket(
    base: &Solid,
    profile_xy: &[(f64, f64)],
    sketch_plane_origin: [f64; 3],
    sketch_plane_normal: [f64; 3],
    depth: f64,
) -> Result<Solid, OcctAdvancedError> {
    if profile_xy.len() < 3 {
        return Err(OcctAdvancedError::bad_input(
            "profile_xy",
            format!(
                "need ≥3 points for a closed profile; got {}",
                profile_xy.len()
            ),
        ));
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "depth",
            "must be positive finite (subtractive depth)",
        ));
    }
    let n_norm = (sketch_plane_normal[0].powi(2)
        + sketch_plane_normal[1].powi(2)
        + sketch_plane_normal[2].powi(2))
    .sqrt();
    if n_norm < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "sketch_plane_normal",
            "must be non-zero",
        ));
    }

    // Subtractive: a negative height makes `feat_make_prism` extrude
    // the cutter *into* the solid (opposite the sketch normal);
    // `fuse = false` removes it.
    valenx_occt_surface::feat_make_prism::feat_make_prism(
        base,
        profile_xy,
        sketch_plane_origin,
        sketch_plane_normal,
        -depth, // negative → cutter pokes into the base
        false,  // subtract — subtractive
    )
    .map_err(|e| OcctAdvancedError::Backend(format!("feat_make_dsubtract_pocket: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    fn square() -> Vec<(f64, f64)> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    }

    #[test]
    fn rejects_short_profile() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_dsubtract_pocket(
            &s,
            &[(0.0, 0.0), (1.0, 0.0)],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.5,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_negative_depth() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err =
            feat_make_dsubtract_pocket(&s, &square(), [0.0; 3], [0.0, 0.0, 1.0], -0.5).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_depth() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err =
            feat_make_dsubtract_pocket(&s, &square(), [0.0; 3], [0.0, 0.0, 1.0], 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_normal() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_dsubtract_pocket(&s, &square(), [0.0; 3], [0.0; 3], 0.5).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn subtractive_pocket_carves_the_base() {
        // A pocket cut into a box from its top face — the result is
        // valid non-empty carved geometry.
        let base = box_solid(6.0, 6.0, 3.0).unwrap();
        let result = feat_make_dsubtract_pocket(
            &base,
            &[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)],
            [0.0, 0.0, 3.0], // sketch plane on the box top
            [0.0, 0.0, 1.0],
            1.5, // pocket depth
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty(), "pocket result should be non-empty");
    }
}
