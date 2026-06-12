//! Phase 133 — `BRepFeat_MakePrism` driven by a sketch (rather than
//! an explicit profile wire).
//!
//! ## What OCCT does
//!
//! `BRepFeat_MakePrism(base, sketch, sketch_plane, direction, fuse,
//! modify)` takes a 2D **sketch** as input — a closed wire on a
//! parametric plane — and pads/pockets it into `base`. Block 1's
//! `valenx_occt_surface::feat_make_prism()` takes a flat profile
//! polyline; this phase consumes a sketch entity directly, which is
//! how FreeCAD Part Design's "Pad" / "Pocket" tools and SolidWorks'
//! "Extruded Boss/Cut" actually feed the kernel. The win is
//! re-parameterisation: edit the sketch dimensions and the feature
//! tree re-evaluates the pad with the new boundary.
//!
//! ## v1 status — real sketch-driven prism
//!
//! A genuine sketch-driven pad/pocket. The pipeline:
//!
//! 1. Extract the sketch's closed wire as a profile polyline via
//!    [`valenx_sketch::extrude::extract_profile_lines`] — the same
//!    accessor the Phase-13 Pad/Pocket feature ops use.
//! 2. Delegate to the real Phase-97
//!    [`feat_make_prism`](fn@valenx_occt_surface::feat_make_prism), which builds the
//!    prism, orients it onto the sketch plane, and fuses (pad) or
//!    subtracts (pocket) it into `base`.
//!
//! What this is *not* is the full `BRepFeat_MakePrism` topology graph
//! with `UpToFirst` / `UpToLast` depth resolution against the base's
//! faces — that needs a feature-tree topology engine. This ships the
//! explicit-depth sketch-driven boss/pocket (SolidWorks "Blind"
//! extrude), which is the overwhelmingly common case.

use valenx_cad::Solid;
use valenx_sketch::Sketch;

use crate::error::OcctAdvancedError;

/// Profile-extraction tolerance — two wire endpoints within this
/// distance are treated as connected. Matches the value the Phase-13
/// feature ops pass.
pub(crate) const SKETCH_PROFILE_TOL: f64 = 1e-6;

/// Apply a sketch-driven prism (pad / pocket) to `base`.
///
/// `sketch` is a `valenx-sketch` 2D sketch — its closed line-loop is
/// extracted as the profile. `sketch_plane_origin` /
/// `sketch_plane_normal` orient the sketch in world space; `length` is
/// the extrusion depth (positive = along normal); `fuse_subtract`
/// chooses pad (true) vs pocket (false).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for a zero / non-finite length,
///   a zero plane normal, or a sketch whose profile is not a closed
///   loop of ≥3 line segments.
/// - [`OcctAdvancedError::Backend`] if the prism cannot be built.
pub fn feat_make_prism_with_sketch(
    base: &Solid,
    sketch: &Sketch,
    sketch_plane_origin: [f64; 3],
    sketch_plane_normal: [f64; 3],
    length: f64,
    fuse_subtract: bool,
) -> Result<Solid, OcctAdvancedError> {
    if !length.is_finite() || length.abs() < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "length",
            "must be non-zero finite",
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
    let profile = sketch_profile(sketch)?;
    valenx_occt_surface::feat_make_prism::feat_make_prism(
        base,
        &profile,
        sketch_plane_origin,
        sketch_plane_normal,
        length,
        fuse_subtract,
    )
    .map_err(|e| OcctAdvancedError::Backend(format!("feat_make_prism_with_sketch: {e}")))
}

/// Extract a sketch's closed profile as an `(x, y)` polyline, mapping
/// the sketcher's error model onto [`OcctAdvancedError::BadInput`].
///
/// Shared by [`feat_make_prism_with_sketch`] and
/// [`crate::feat_make_revol_with_sketch`].
pub(crate) fn sketch_profile(sketch: &Sketch) -> Result<Vec<(f64, f64)>, OcctAdvancedError> {
    let profile = valenx_sketch::extrude::extract_profile_lines(sketch, SKETCH_PROFILE_TOL)
        .map_err(|e| {
            OcctAdvancedError::bad_input(
                "sketch",
                format!("could not extract a closed profile from the sketch: {e}"),
            )
        })?;
    if profile.len() < 3 {
        return Err(OcctAdvancedError::bad_input(
            "sketch",
            format!(
                "sketch profile needs ≥3 line segments forming a closed loop; got {}",
                profile.len()
            ),
        ));
    }
    Ok(profile)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    /// A `side`×`side` square sketch with its lower-left corner at the
    /// origin.
    fn square_sketch(side: f64) -> Sketch {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(side, 0.0);
        let c = s.add_point(side, side);
        let d = s.add_point(0.0, side);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    fn rejects_zero_length() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_prism_with_sketch(
            &base,
            &square_sketch(0.5),
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.0,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_normal() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err =
            feat_make_prism_with_sketch(&base, &square_sketch(0.5), [0.0; 3], [0.0; 3], 0.5, true)
                .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_empty_sketch() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_prism_with_sketch(
            &base,
            &Sketch::new(),
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.5,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn sketch_driven_boss_adds_material() {
        // Pad a square sketch onto a box's top face — the result must
        // rise above the base.
        let base = box_solid(4.0, 4.0, 1.0).unwrap();
        let result = feat_make_prism_with_sketch(
            &base,
            &square_sketch(2.0),
            [0.0, 0.0, 1.0], // sketch plane on the box top
            [0.0, 0.0, 1.0],
            0.75,
            true,
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
            "sketch-driven boss should rise above z=1"
        );
    }

    #[test]
    fn sketch_driven_pocket_removes_material() {
        let base = box_solid(6.0, 6.0, 2.0).unwrap();
        let result = feat_make_prism_with_sketch(
            &base,
            &square_sketch(2.0),
            [1.0, 1.0, 2.0], // sketch plane on the box top
            [0.0, 0.0, 1.0],
            -1.0, // negative → cut downward
            false,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.1).unwrap();
        assert!(!mesh.nodes.is_empty(), "pocket result should be non-empty");
    }
}
