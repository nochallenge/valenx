//! Phase 76 — `GeomFill_SectionLaw`: sweep with a section that
//! evolves along the path.
//!
//! ## What OCCT does
//!
//! `GeomFill_SectionLaw` is the abstract base for "section evolution
//! laws" used by `GeomFill_Sweep`. The classical use case is a sweep
//! where the cross-section morphs from one shape at the start of the
//! path to a different shape at the end — for instance a rocket-
//! engine nozzle whose throat is circular at the inlet and elliptical
//! at the outlet. Three concrete implementations ship with OCCT:
//!
//! - `GeomFill_UniformSection` — constant cross-section (degenerate
//!   case; equivalent to `BRepPrimAPI_MakePrism`).
//! - `GeomFill_EvolvedSection` — linear interpolation between two
//!   `Geom_Curve`s.
//! - `GeomFill_NSections` — N-way blend through M intermediate
//!   sections, with the per-section parameter values supplied by the
//!   caller.
//!
//! Combined with a path curve, the result is a `Geom_BSplineSurface`
//! that represents the swept body. This is the geometric engine
//! behind [`crate::pipe_shell()`] and OCCT's `BRepOffsetAPI_MakePipeShell`.
//!
//! ## v1 status — real implementation (N-section skinning)
//!
//! This module implements the `GeomFill_NSections` case directly: a
//! tensor-product B-spline surface skinned through the supplied
//! cross-sections.
//!
//! 1. Each section [`NurbsCurve`] is **sampled** at a common count
//!    `SECTION_SAMPLES` of points across its parameter range. The
//!    samples become one row of a structured `(N_sections ×
//!    SECTION_SAMPLES)` data grid.
//! 2. The grid is fitted with [`valenx_surface::fit::nurbs_surface_through_grid`]
//!    — a v-direction interpolation through the section samples and a
//!    u-direction interpolation across the sections. The surface's
//!    v parameter runs *along* the path (section evolution); u runs
//!    *around* the cross-section.
//! 3. `section_params` order the sections; their values are validated
//!    strictly increasing so the skinning direction is well-defined.
//!    (The numeric values themselves seed the v-direction
//!    parameterisation conceptually — the fit uses a uniform
//!    parameterisation, a documented simplification.)
//!
//! The `path` curve is accepted for API parity with OCCT's
//! `GeomFill_Sweep`; in the pure section-law (skinning) form the
//! surface is fully determined by the sections, so `path` only
//! affects the result through the *placement* of the sections, which
//! the caller bakes into the section geometry. A spine-driven sweep
//! (where the path also bends the section frames) is
//! [`crate::sweep_api_pipe_shell()`].
//!
//! ### Limits
//!
//! - Sections must already sit in 3D where the caller wants them
//!   (this is skinning, not transport — the frame-transport variant
//!   is `sweep_api_pipe_shell`).
//! - At least **2** sections are required: a single section defines
//!   no evolution and cannot skin a surface.

use valenx_surface::{NurbsCurve, NurbsSurface};

use crate::error::OcctSurfaceError;

/// Number of points each section curve is sampled at to form a row of
/// the skinning grid. Higher = more faithful section shape, denser
/// surface.
pub const SECTION_SAMPLES: usize = 16;

/// Sweep `sections` along `path`, interpolating the cross-section
/// shape at each path parameter (the `GeomFill_NSections` skinning
/// law).
///
/// `section_params` carries the t-values along `path` at which each
/// section is sampled — must be the same length as `sections` and
/// strictly increasing.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] for malformed inputs (empty /
/// single section, length mismatch, non-increasing params, or a
/// section / fit that fails to build).
pub fn geom_fill_section_law(
    sections: &[NurbsCurve],
    section_params: &[f64],
    _path: &NurbsCurve,
) -> Result<NurbsSurface, OcctSurfaceError> {
    if sections.is_empty() {
        return Err(OcctSurfaceError::bad_input(
            "sections",
            "need at least one cross-section",
        ));
    }
    if sections.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "sections",
            "skinning needs at least 2 cross-sections to define an \
             evolution; for a constant section use a prism sweep",
        ));
    }
    if sections.len() != section_params.len() {
        return Err(OcctSurfaceError::bad_input(
            "section_params",
            format!(
                "must match sections length ({} vs {})",
                section_params.len(),
                sections.len()
            ),
        ));
    }
    // Strictly-increasing check.
    for w in section_params.windows(2) {
        if w[0] >= w[1] {
            return Err(OcctSurfaceError::bad_input(
                "section_params",
                "must be strictly increasing",
            ));
        }
    }

    // Sample every section at SECTION_SAMPLES points → one grid row.
    let mut grid: Vec<Vec<nalgebra::Vector3<f64>>> = Vec::with_capacity(sections.len());
    for section in sections {
        let (t0, t1) = section.parameter_range();
        let mut row = Vec::with_capacity(SECTION_SAMPLES);
        for k in 0..SECTION_SAMPLES {
            let t = if SECTION_SAMPLES == 1 {
                t0
            } else {
                t0 + (t1 - t0) * (k as f64 / (SECTION_SAMPLES - 1) as f64)
            };
            row.push(section.evaluate(t));
        }
        grid.push(row);
    }

    // Choose surface degrees + control-point counts that the grid can
    // support. v runs along the cross-section samples; u across the
    // sections.
    let degree_u = pick_degree(sections.len());
    let degree_v = pick_degree(SECTION_SAMPLES);
    // Fit through every sample (one CP per data point) for a faithful
    // interpolating skin.
    let n_cps_u = sections.len();
    let n_cps_v = SECTION_SAMPLES;

    let fit = valenx_surface::fit::nurbs_surface_through_grid(
        &grid, degree_u, degree_v, n_cps_u, n_cps_v,
    )
    .map_err(|e| OcctSurfaceError::bad_input("sections", format!("skinning fit failed: {e}")))?;
    Ok(fit.surface)
}

/// Pick a B-spline degree appropriate for `n` control points: cubic
/// when there's room, otherwise the highest degree the count allows
/// (degree ≤ n - 1), clamped to ≥ 1.
fn pick_degree(n: usize) -> usize {
    3.min(n.saturating_sub(1)).max(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn line(a: Vector3<f64>, b: Vector3<f64>) -> NurbsCurve {
        NurbsCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![a, b], vec![1.0, 1.0]).unwrap()
    }

    /// A square cross-section as a degree-1 closed-ish polyline curve
    /// at height `z`, scaled by `scale`.
    fn square_section(z: f64, scale: f64) -> NurbsCurve {
        let s = scale;
        let pts = vec![
            Vector3::new(-s, -s, z),
            Vector3::new(s, -s, z),
            Vector3::new(s, s, z),
            Vector3::new(-s, s, z),
            Vector3::new(-s, -s, z),
        ];
        let n = pts.len();
        // Open-uniform degree-1 knot vector.
        let mut knots = vec![0.0; n + 2];
        for (i, k) in knots.iter_mut().enumerate() {
            *k = (i.saturating_sub(1)).min(n - 1) as f64;
        }
        let weights = vec![1.0; n];
        NurbsCurve::new(1, knots, pts, weights).unwrap()
    }

    #[test]
    fn section_law_rejects_empty_sections() {
        let path = line(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0));
        let err = geom_fill_section_law(&[], &[], &path).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn section_law_rejects_single_section() {
        // A single section defines no evolution — must be rejected.
        let s = square_section(0.0, 1.0);
        let path = line(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0));
        let err = geom_fill_section_law(std::slice::from_ref(&s), &[0.0], &path).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn section_law_rejects_non_increasing_params() {
        let s1 = square_section(0.0, 1.0);
        let s2 = square_section(1.0, 1.0);
        let path = line(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0));
        let err = geom_fill_section_law(&[s1, s2], &[0.5, 0.3], &path).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn section_law_rejects_length_mismatch() {
        let s1 = square_section(0.0, 1.0);
        let s2 = square_section(1.0, 1.0);
        let path = line(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0));
        let err = geom_fill_section_law(&[s1, s2], &[0.0], &path).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn section_law_skins_two_sections_into_a_surface() {
        // A square at z=0 and a bigger square at z=2 — skinning gives
        // a real frustum-shell surface.
        let s1 = square_section(0.0, 1.0);
        let s2 = square_section(2.0, 2.0);
        let path = line(Vector3::zeros(), Vector3::new(0.0, 0.0, 2.0));
        let surf = geom_fill_section_law(&[s1, s2], &[0.0, 1.0], &path).unwrap();
        // The surface spans both sections in z.
        let bottom = surf.evaluate(0.0, 0.5);
        let top = surf.evaluate(1.0, 0.5);
        assert!(bottom.z.abs() < 0.3, "bottom z ≈ 0, got {}", bottom.z);
        assert!((top.z - 2.0).abs() < 0.3, "top z ≈ 2, got {}", top.z);
        // The top section is wider than the bottom — a point on the
        // top edge is further from the z-axis than the bottom.
        let bottom_r = (bottom.x * bottom.x + bottom.y * bottom.y).sqrt();
        let top_r = (top.x * top.x + top.y * top.y).sqrt();
        assert!(
            top_r > bottom_r,
            "top section should be wider: top_r={top_r}, bottom_r={bottom_r}"
        );
    }

    #[test]
    fn section_law_handles_three_sections() {
        let s1 = square_section(0.0, 1.0);
        let s2 = square_section(1.0, 1.5);
        let s3 = square_section(2.0, 0.5);
        let path = line(Vector3::zeros(), Vector3::new(0.0, 0.0, 2.0));
        let surf = geom_fill_section_law(&[s1, s2, s3], &[0.0, 0.5, 1.0], &path).unwrap();
        // Three sections → a valid surface; the middle bulges out.
        assert!(surf.nu() >= 3);
        let mid = surf.evaluate(0.5, 0.5);
        assert!(mid.z > 0.5 && mid.z < 1.5, "mid z in band, got {}", mid.z);
    }
}
