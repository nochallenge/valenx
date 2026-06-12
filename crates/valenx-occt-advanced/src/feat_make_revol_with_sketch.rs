//! Phase 134 — `BRepFeat_MakeRevol` driven by a sketch.
//!
//! ## What OCCT does
//!
//! Sketch-driven analog of [`crate::feat_make_prism_with_sketch()`] —
//! takes a 2D sketch on a parametric plane, sweeps it around an axis
//! by `angle_rad` (typically 2π for a full revolution), and fuses or
//! subtracts the result into `base`. Maps to FreeCAD Part Design
//! "Revolution" / "Groove", SolidWorks "Revolved Boss/Base" /
//! "Revolved Cut".
//!
//! ## v1 status — real sketch-driven revolution
//!
//! A genuine sketch-driven revolved boss/groove. The pipeline:
//!
//! 1. Extract the sketch's closed wire as a profile polyline via
//!    [`valenx_sketch::extrude::extract_profile_lines`].
//! 2. Delegate to the real Phase-98
//!    [`feat_make_revol`](fn@valenx_occt_surface::feat_make_revol), which revolves the
//!    profile around the axis and fuses (boss) or subtracts (groove)
//!    it into `base`.
//!
//! Same scope note as [`crate::feat_make_prism_with_sketch()`]: this
//! is the explicit-angle revolved feature (the common case), not the
//! full `BRepFeat` topology graph with up-to-face depth resolution.
//!
//! Note the sketch profile lives in its own XY plane; for a valid
//! solid of revolution it must lie entirely on one side of the
//! revolution axis (the same constraint OCCT imposes).

use valenx_cad::Solid;
use valenx_sketch::Sketch;

use crate::error::OcctAdvancedError;
use crate::feat_make_prism_with_sketch::sketch_profile;

/// Apply a sketch-driven revolution (boss / groove) to `base`.
///
/// `sketch` is a `valenx-sketch` 2D sketch — its closed line-loop is
/// the profile. `axis_origin` / `axis_direction` define the
/// revolution axis in world space; `angle_rad` is the sweep angle
/// (positive = right-hand rule); `fuse_subtract` chooses boss-revolve
/// (true) vs cut-revolve (false).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for a zero / non-finite angle, a
///   zero axis, or a sketch whose profile is not a closed loop of ≥3
///   line segments.
/// - [`OcctAdvancedError::Backend`] if the revolution cannot be built.
pub fn feat_make_revol_with_sketch(
    base: &Solid,
    sketch: &Sketch,
    axis_origin: [f64; 3],
    axis_direction: [f64; 3],
    angle_rad: f64,
    fuse_subtract: bool,
) -> Result<Solid, OcctAdvancedError> {
    if !angle_rad.is_finite() || angle_rad.abs() < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "angle_rad",
            "must be non-zero finite",
        ));
    }
    let a_norm =
        (axis_direction[0].powi(2) + axis_direction[1].powi(2) + axis_direction[2].powi(2)).sqrt();
    if a_norm < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "axis_direction",
            "must be non-zero",
        ));
    }
    let profile = sketch_profile(sketch)?;
    valenx_occt_surface::feat_make_revol::feat_make_revol(
        base,
        &profile,
        axis_origin,
        axis_direction,
        angle_rad,
        fuse_subtract,
    )
    .map_err(|e| OcctAdvancedError::Backend(format!("feat_make_revol_with_sketch: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;
    use valenx_cad::box_solid;

    /// A unit-square sketch offset to `x ∈ [x0, x0+1]`, `y ∈ [0, 1]`
    /// so it sits entirely off a Y-axis revolution axis.
    fn offset_square_sketch(x0: f64) -> Sketch {
        let mut s = Sketch::new();
        let a = s.add_point(x0, 0.0);
        let b = s.add_point(x0 + 1.0, 0.0);
        let c = s.add_point(x0 + 1.0, 1.0);
        let d = s.add_point(x0, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    fn rejects_zero_angle() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_revol_with_sketch(
            &base,
            &offset_square_sketch(2.0),
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.0,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_axis() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_revol_with_sketch(
            &base,
            &offset_square_sketch(2.0),
            [0.0; 3],
            [0.0; 3],
            TAU,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_empty_sketch() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_revol_with_sketch(
            &base,
            &Sketch::new(),
            [0.0; 3],
            [0.0, 1.0, 0.0],
            TAU,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn sketch_driven_revolved_boss_produces_geometry() {
        // Revolve an off-axis square sketch a full turn around Y,
        // fusing the resulting ring onto a base box.
        let base = box_solid(8.0, 2.0, 8.0)
            .unwrap()
            .translated(-4.0, 0.0, -4.0)
            .unwrap();
        let result = feat_make_revol_with_sketch(
            &base,
            &offset_square_sketch(2.0),
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            TAU,
            true,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty(), "revolved boss should be non-empty");
    }
}
