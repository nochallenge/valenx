//! Phase 131 — `BRepOffsetAPI_ThruSections` with explicit guide
//! curves controlling intermediate cross-sections.
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_ThruSections(is_solid, is_ruled, presision)` builds
//! a loft surface through an ordered sequence of cross-section wires.
//! The guide-curve variant takes additional `Geom_Curve` rails that
//! the swept surface must pass through tangentially at the
//! intermediate `v` parameters between sections — this is how a
//! turbine-blade leading edge is built (sections at hub / mid /
//! shroud, guide curves along the LE / TE / suction-side / pressure-
//! side). Without guides, the loft picks the smoothest path between
//! sections and intermediate cross-sections may bulge in unwanted
//! directions.
//!
//! Maps to FreeCAD Part Design "Additive Loft" with "Rails" option,
//! SolidWorks "Lofted Boss/Base" with "Guide Curves", Inventor "Loft"
//! with "Rails".
//!
//! Block 1's `valenx_occt_surface::sweep_api_thru_sections()` does
//! the bare `is_solid`/`is_ruled` form; this phase adds the guides.
//!
//! ## v1 status — real guide-warped loft
//!
//! This is a genuine guide-constrained loft. Built on the Phase-90
//! loft, then each smooth intermediate ring is warped so it follows
//! the guide curves at the matching arc-length station — see
//! [`crate::guide_loft`] for the full algorithm and its honest scope.
//! The cross-section rings that define the loft (the input sections)
//! are left un-warped so the surface still interpolates them exactly.
//!
//! Honest scope (also documented on [`crate::guide_loft`]): the
//! result is **mesh-domain** (no BRep faces) and the guide warp is a
//! rigid-translate + uniform-radial-scale, not the per-control-point
//! tangential `Geom_BSplineSurface` skinning OCCT runs. Exact surface
//! tangency at each rail crossing stays a Tier-3 follow-up gated on
//! the BRep substrate.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;
use crate::guide_loft::guide_loft_solid;

/// 3-component point used as both `sections` (closed polylines per
/// cross-section) and `guides` (open polylines per rail).
pub type Pt3 = [f64; 3];

/// Build a loft solid through `sections` constrained to pass through
/// `guides` at intermediate parameters.
///
/// Each entry in `sections` is one closed cross-section wire as a
/// polyline of 3D points; each entry in `guides` is one open rail
/// polyline. Real OCCT requires `guides.len() >= 1` (a single rail
/// is the minimum useful case — pinch a soap-bubble loft to a single
/// curve), and each guide must intersect every section.
///
/// `is_solid` adds triangulated planar end caps so the result is a
/// closed (watertight) mesh-backed [`Solid`].
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for malformed inputs (fewer than
///   2 sections, a section with < 3 points, a guide with < 2 points,
///   or a non-finite coordinate).
pub fn offset_api_thru_sections_with_guides(
    sections: &[Vec<Pt3>],
    guides: &[Vec<Pt3>],
) -> Result<Solid, OcctAdvancedError> {
    offset_api_thru_sections_with_guides_ex(sections, guides, true)
}

/// Guide-warped loft with explicit `is_solid` control.
///
/// `is_solid = true` caps the ends (closed solid); `is_solid = false`
/// returns the open lofted shell.
///
/// # Errors
///
/// [`OcctAdvancedError::BadInput`] — see
/// [`offset_api_thru_sections_with_guides`].
pub fn offset_api_thru_sections_with_guides_ex(
    sections: &[Vec<Pt3>],
    guides: &[Vec<Pt3>],
    is_solid: bool,
) -> Result<Solid, OcctAdvancedError> {
    validate(sections, guides)?;
    Ok(guide_loft_solid(sections, guides, is_solid))
}

/// Shared input validation for the guide-loft entry points (and
/// [`feat_make_loft_with_rails`](fn@crate::feat_make_loft_with_rails),
/// which has the same shape of inputs).
pub(crate) fn validate(
    sections: &[Vec<Pt3>],
    guides: &[Vec<Pt3>],
) -> Result<(), OcctAdvancedError> {
    if sections.len() < 2 {
        return Err(OcctAdvancedError::bad_input(
            "sections",
            format!(
                "need at least 2 cross-sections to loft; got {}",
                sections.len()
            ),
        ));
    }
    if guides.is_empty() {
        return Err(OcctAdvancedError::bad_input(
            "guides",
            "guide-curve variant requires at least 1 rail; \
             for guide-free loft use valenx_occt_surface::sweep_api_thru_sections",
        ));
    }
    for (i, sec) in sections.iter().enumerate() {
        if sec.len() < 3 {
            return Err(OcctAdvancedError::bad_input(
                "sections",
                format!("sections[{i}] must have ≥3 points; got {}", sec.len()),
            ));
        }
        for p in sec {
            if p.iter().any(|c| !c.is_finite()) {
                return Err(OcctAdvancedError::bad_input(
                    "sections",
                    format!("sections[{i}] has a non-finite coordinate"),
                ));
            }
        }
    }
    for (i, g) in guides.iter().enumerate() {
        if g.len() < 2 {
            return Err(OcctAdvancedError::bad_input(
                "guides",
                format!("guides[{i}] must have ≥2 points; got {}", g.len()),
            ));
        }
        for p in g {
            if p.iter().any(|c| !c.is_finite()) {
                return Err(OcctAdvancedError::bad_input(
                    "guides",
                    format!("guides[{i}] has a non-finite coordinate"),
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(z: f64, half: f64) -> Vec<Pt3> {
        vec![
            [-half, -half, z],
            [half, -half, z],
            [half, half, z],
            [-half, half, z],
        ]
    }

    #[test]
    fn rejects_one_section() {
        let secs = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]]];
        let guides = vec![vec![[0.5, 0.5, 0.0], [0.5, 0.5, 1.0]]];
        let err = offset_api_thru_sections_with_guides(&secs, &guides).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_guides() {
        let secs = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]],
        ];
        let err = offset_api_thru_sections_with_guides(&secs, &[]).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_degenerate_section() {
        let secs = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]],
        ];
        let guides = vec![vec![[0.5, 0.5, 0.0], [0.5, 0.5, 1.0]]];
        let err = offset_api_thru_sections_with_guides(&secs, &guides).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_degenerate_guide() {
        let secs = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]],
        ];
        let guides = vec![vec![[0.5, 0.5, 0.0]]];
        let err = offset_api_thru_sections_with_guides(&secs, &guides).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_non_finite_coordinate() {
        let secs = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, f64::NAN], [0.5, 1.0, 0.0]],
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]],
        ];
        let guides = vec![vec![[0.5, 0.5, 0.0], [0.5, 0.5, 1.0]]];
        let err = offset_api_thru_sections_with_guides(&secs, &guides).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn lofts_two_squares_with_a_guide() {
        // A real guide-warped loft of two squares — the result
        // tessellates to non-empty geometry spanning both heights.
        let secs = vec![square(0.0, 1.0), square(4.0, 1.0)];
        let guides = vec![vec![[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]]];
        let solid = offset_api_thru_sections_with_guides(&secs, &guides).unwrap();
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty());
        let zmin = mesh.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmin.abs() < 1e-6, "loft should reach z=0, got {zmin}");
        assert!((zmax - 4.0).abs() < 1e-6, "loft should reach z=4, got {zmax}");
    }

    #[test]
    fn guide_curvature_bends_the_loft() {
        // A guide that bows out to +X pulls the loft's mid section
        // toward +X — the mesh's max X exceeds the un-warped square's
        // half-width of 1.0.
        let secs = vec![square(0.0, 1.0), square(6.0, 1.0)];
        let bowed = vec![vec![
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 3.0],
            [0.0, 0.0, 6.0],
        ]];
        let solid = offset_api_thru_sections_with_guides(&secs, &bowed).unwrap();
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.2).unwrap();
        let xmax = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            xmax > 1.5,
            "the bowed guide should push the loft past x=1.5, got xmax={xmax}"
        );
    }

    #[test]
    fn open_shell_variant_skips_caps() {
        let secs = vec![square(0.0, 1.0), square(2.0, 1.0)];
        let guides = vec![vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]]];
        let solid_capped = offset_api_thru_sections_with_guides_ex(&secs, &guides, true).unwrap();
        let shell = offset_api_thru_sections_with_guides_ex(&secs, &guides, false).unwrap();
        let cap_tris = valenx_cad::solid_to_mesh(&solid_capped, 0.2)
            .unwrap()
            .total_elements();
        let shell_tris = valenx_cad::solid_to_mesh(&shell, 0.2)
            .unwrap()
            .total_elements();
        assert!(cap_tris > shell_tris, "capped solid should have more triangles");
    }
}
