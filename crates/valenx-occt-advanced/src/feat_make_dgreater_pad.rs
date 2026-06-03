//! Phase 136 — `BRepFeat_MakeDPrism` "additive" / pad variant
//! (feature-based pad larger than the sketch profile).
//!
//! ## What OCCT does
//!
//! The "DGreater" suffix in OCCT-land flags the additive `MakeDPrism`
//! mode where the prism *grows* the base solid (boss-extrude) — vs the
//! subtractive `DSubtract` (Phase 137) which carves a pocket. The
//! distinguishing feature: the prism height can exceed the sketch
//! profile's planar footprint so the resulting feature wraps around
//! the base's curved boundaries. Used for adding a circular boss to a
//! cylindrical shaft where the boss must fully encircle the shaft, or
//! for adding ribs that wrap around fillets.
//!
//! ## v1 status — real additive pad
//!
//! A genuine additive boss-extrude. The pipeline:
//!
//! 1. Build the prism body from `profile_xy` extruded by `length`
//!    along `+Z` ([`prim_api_prism`](fn@valenx_occt_surface::prim_api_prism)).
//! 2. Orient it onto the world sketch plane.
//! 3. **Fuse** it onto `base` via the real Phase-97
//!    [`valenx_occt_surface::feat_support::feature_combine`] — an
//!    exact BRep boolean when both operands are BReps, with a
//!    co-refinement mesh-CSG fallback otherwise.
//!
//! This delegates to the real Phase-97
//! [`feat_make_prism`](fn@valenx_occt_surface::feat_make_prism) with `fuse = true`. The
//! `length`-may-exceed-the-footprint contract is honoured: a tall pad
//! built on a curved base genuinely fuses through the base's
//! curvature via the boolean. What this does *not* model is the
//! `BRepFeat_MakeDPrism` topology graph that *redirects* the prism's
//! side walls along the base's curved faces (so the boss's footprint
//! literally bends with the shaft) — that needs the BRep feature
//! topology engine and stays Tier-3. The boolean-fused oversized pad
//! is the common case and is real, carved geometry.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Apply an additive "DGreater" pad to `base`.
///
/// `profile_xy` is the closed sketch profile in the sketch-plane XY;
/// `sketch_plane_origin` / `sketch_plane_normal` orient the plane.
/// `length` is the pad height — it may exceed the profile's planar
/// extent (the DGreater regime).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for a too-short profile, a
///   non-positive / non-finite length, or a zero plane normal.
/// - [`OcctAdvancedError::Backend`] if the pad cannot be built.
pub fn feat_make_dgreater_pad(
    base: &Solid,
    profile_xy: &[(f64, f64)],
    sketch_plane_origin: [f64; 3],
    sketch_plane_normal: [f64; 3],
    length: f64,
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
    if !length.is_finite() || length <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "length",
            "must be positive finite (pad grows outward)",
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

    // Additive: fuse a prism of the (possibly oversized) length.
    valenx_occt_surface::feat_make_prism::feat_make_prism(
        base,
        profile_xy,
        sketch_plane_origin,
        sketch_plane_normal,
        length,
        true, // fuse — additive
    )
    .map_err(|e| OcctAdvancedError::Backend(format!("feat_make_dgreater_pad: {e}")))
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
        let err = feat_make_dgreater_pad(
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
    fn rejects_zero_length() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err =
            feat_make_dgreater_pad(&s, &square(), [0.0; 3], [0.0, 0.0, 1.0], 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_negative_length() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err =
            feat_make_dgreater_pad(&s, &square(), [0.0; 3], [0.0, 0.0, 1.0], -0.5).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_normal() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_dgreater_pad(&s, &square(), [0.0; 3], [0.0; 3], 0.5).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn additive_pad_grows_the_base() {
        // A tall pad (length 3 > the 1×1 footprint — the DGreater
        // regime) fused onto a box's top face must rise above it.
        let base = box_solid(4.0, 4.0, 1.0).unwrap();
        let result = feat_make_dgreater_pad(
            &base,
            &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
            [0.0, 0.0, 1.0], // sketch plane on the box top
            [0.0, 0.0, 1.0],
            3.0, // length far exceeds the 2×2 footprint
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.2).unwrap();
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            zmax > 3.5,
            "the oversized pad should rise well above the base, got zmax={zmax}"
        );
    }
}
