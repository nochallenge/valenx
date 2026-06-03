//! Phase 138 — feature-based loft with guide rails (parallel to
//! Phase 13's LoftWithRails feature-tree node).
//!
//! ## What OCCT does
//!
//! Combines [`crate::offset_api_thru_sections_with_guides()`] with
//! `BRepFeat_MakeRevol`-style topology integration: takes a list of
//! cross-section sketches *plus* a list of guide rails, builds the
//! loft surface, closes it into a solid, and fuses or subtracts the
//! result into `base`. Vs the bare ThruSections+guides
//! ([`crate::offset_api_thru_sections_with_guides()`]) which returns a
//! standalone solid, this op participates in the host solid's feature
//! tree (sharing edges where the loft butts up against existing faces).
//!
//! Maps to FreeCAD Part Design "Additive Loft" with rails enabled +
//! Pad-style integration, SolidWorks "Lofted Boss/Base" with guide
//! curves applied to a base feature.
//!
//! ## v1 status — real feature-based guide loft
//!
//! A genuine boss/pocket guide-loft. The pipeline:
//!
//! 1. Build the guide-warped loft *solid* via the same
//!    [`crate::guide_loft`] machinery Phase 131 uses (capped, so the
//!    loft is a closed body).
//! 2. **Fuse** (boss) or **subtract** (pocket) it with `base` via the
//!    real Phase-97 [`valenx_occt_surface::feat_support::feature_combine`]
//!    — exact BRep boolean when both operands are BReps, co-refinement
//!    mesh CSG otherwise (the loft is mesh-backed, so the mesh-CSG
//!    path is the one taken).
//!
//! Honest scope: inherits [`crate::guide_loft`]'s mesh-domain limit.
//! The loft has no BRep faces and the guide warp is a rigid-translate
//! plus uniform-radial-scale, not OCCT's per-control-point skinning.
//! The full `BRepFeat` topology graph (loft sharing parametric edges
//! with the base's faces) stays Tier-3.

use valenx_cad::Solid;
use valenx_occt_surface::feat_support::feature_combine;

use crate::error::OcctAdvancedError;
use crate::guide_loft::guide_loft_solid;
use crate::offset_api_thru_sections_with_guides::validate;

/// 3-component point used for both `sections` (closed cross-section
/// polylines) and `rails` (open guide polylines).
pub type Pt3 = [f64; 3];

/// Apply a loft-with-rails feature to `base`.
///
/// `sections` — at least 2 closed cross-section polylines.
/// `rails` — at least 1 guide-rail polyline that intersects every
/// section.
/// `fuse_subtract` — true for additive loft (boss), false for
/// subtractive (pocket).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for malformed inputs.
/// - [`OcctAdvancedError::Backend`] if the feature combine fails even
///   the mesh-CSG fallback.
pub fn feat_make_loft_with_rails(
    base: &Solid,
    sections: &[Vec<Pt3>],
    rails: &[Vec<Pt3>],
    fuse_subtract: bool,
) -> Result<Solid, OcctAdvancedError> {
    validate(sections, rails)?;
    // The loft must be a closed body to fuse/subtract — cap the ends.
    let loft = guide_loft_solid(sections, rails, true);
    feature_combine(base, &loft, fuse_subtract)
        .map_err(|e| OcctAdvancedError::Backend(format!("loft-with-rails combine: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    fn square(z: f64, half: f64) -> Vec<Pt3> {
        vec![
            [-half, -half, z],
            [half, -half, z],
            [half, half, z],
            [-half, half, z],
        ]
    }

    fn make_valid_inputs() -> (Vec<Vec<Pt3>>, Vec<Vec<Pt3>>) {
        let secs = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]],
        ];
        let rails = vec![vec![[0.5, 0.5, 0.0], [0.5, 0.5, 1.0]]];
        (secs, rails)
    }

    #[test]
    fn rejects_one_section() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let (_, rails) = make_valid_inputs();
        let secs = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]]];
        let err = feat_make_loft_with_rails(&base, &secs, &rails, true).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_rails() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let (secs, _) = make_valid_inputs();
        let err = feat_make_loft_with_rails(&base, &secs, &[], true).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_degenerate_section() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let (_, rails) = make_valid_inputs();
        let secs = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]],
        ];
        let err = feat_make_loft_with_rails(&base, &secs, &rails, true).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn boss_loft_produces_valid_geometry() {
        // Fuse a guide-lofted body onto a base box; the result must
        // be non-empty carved geometry.
        let base = box_solid(8.0, 8.0, 1.0).unwrap();
        let secs = vec![square(0.0, 1.5), square(4.0, 1.5)];
        let rails = vec![vec![[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]]];
        let result = feat_make_loft_with_rails(&base, &secs, &rails, true).unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.3).unwrap();
        assert!(!mesh.nodes.is_empty(), "boss-loft result should be non-empty");
        // The loft rises above the 1-tall base.
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmax > 1.0 + 1e-6, "loft boss should rise above the base");
    }

    #[test]
    fn pocket_loft_carves_the_base() {
        let base = box_solid(8.0, 8.0, 6.0).unwrap();
        let secs = vec![square(-1.0, 1.0), square(7.0, 1.0)];
        let rails = vec![vec![[0.0, 0.0, -1.0], [0.0, 0.0, 7.0]]];
        let result = feat_make_loft_with_rails(&base, &secs, &rails, false).unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.3).unwrap();
        assert!(!mesh.nodes.is_empty(), "pocket-loft result should be non-empty");
    }
}
