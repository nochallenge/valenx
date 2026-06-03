//! Phase 135 — `BRepFeat_MakePipe` with path-following constraint.
//!
//! ## What OCCT does
//!
//! `BRepFeat_MakePipe(base, profile, spine, fuse_subtract)` sweeps
//! `profile` along `spine` and integrates the result into `base` as a
//! tube-style feature. The path-following constraint means the profile
//! is transported along the spine with a specific frame law
//! (Frenet / corrected-frenet / discrete-tangent) so it stays normal
//! to the spine — vs an unconstrained sweep that may twist arbitrarily.
//! Used for HVAC ducts, cable routing through a chassis, water-cooling
//! channels through a heat exchanger.
//!
//! Block 1's `valenx_occt_surface::feat_make_pipe()` is the bare
//! "sweep this profile along this path, fuse the result" variant;
//! this phase adds the explicit frame-law control.
//!
//! ## v1 status — real path-constrained pipe feature
//!
//! A genuine path-following swept boss/cut. The pipeline:
//!
//! 1. Sweep the profile along the spine with the explicit
//!    [`FrameLaw`] — the section is locked normal to the spine via a
//!    Bishop / rotation-minimising parallel transport (the Phase-89
//!    [`sweep_api_pipe_shell`](fn@valenx_occt_surface::sweep_api_pipe_shell) machinery, used
//!    here without an auxiliary spine so the only constraint is the
//!    path itself).
//! 2. **Fuse** (boss) or **subtract** (cut) the swept body with
//!    `base` via the real Phase-97
//!    [`valenx_occt_surface::feat_support::feature_combine`].
//!
//! ### How the frame laws map
//!
//! Valenx's sweep transport is a **rotation-minimising frame** — the
//! Bishop frame, which is exactly OCCT's *corrected Frenet*: the
//! section never twists about the tangent beyond what the path
//! curvature forces.
//!
//! - [`FrameLaw::CorrectedFrenet`] is this RMF directly — no roll
//!   accumulation, smooth across inflection points.
//! - [`FrameLaw::DiscreteTangent`] uses the same per-segment minimal
//!   rotation (straight transport between sample points), which is
//!   what the RMF already does on a polyline spine — so it is the RMF
//!   as well.
//! - [`FrameLaw::Frenet`] in OCCT is the *raw* Frenet frame (uses the
//!   curvature normal, flips at inflection points). Valenx has no raw
//!   Frenet sweep; the path-constrained pipe falls back to the RMF
//!   and reports that honestly via the doc here — the RMF is a
//!   strictly better section frame, so the substitution never
//!   produces *worse* geometry, only a non-twisting one.
//!
//! Honest scope: the swept body is **mesh-domain** (the sweep carries
//! no BRep faces, the same rule as every Valenx sweep), so the
//! feature combine takes the co-refinement mesh-CSG path. Raw Frenet
//! and the full `BRepFeat` topology graph stay Tier-3.

use valenx_cad::Solid;
use valenx_occt_surface::feat_support::feature_combine;
use valenx_occt_surface::sweep_api_pipe_shell::sweep_api_pipe_shell;

use crate::error::OcctAdvancedError;

/// Frame transport law for the profile as it moves along the spine.
///
/// Maps 1:1 onto OCCT's `BRepBuilderAPI_TransitionMode`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum FrameLaw {
    /// Frenet frame — profile aligned to spine tangent + normal +
    /// binormal. Fast but flips on inflection points.
    Frenet,
    /// "Corrected" Frenet — smoothed across inflection points by
    /// integrating the binormal rotation. Default in most CAD systems.
    CorrectedFrenet,
    /// Discrete tangent — profile rotated only at sample points,
    /// straight transport in between. Cheapest, suitable for polyline
    /// spines.
    DiscreteTangent,
}

impl FrameLaw {
    /// Whether this law is delivered exactly by Valenx's
    /// rotation-minimising frame. `CorrectedFrenet` and
    /// `DiscreteTangent` are; raw `Frenet` is not (Valenx substitutes
    /// the RMF — see the module docs).
    pub fn is_exact_in_v1(self) -> bool {
        !matches!(self, FrameLaw::Frenet)
    }
}

/// Apply a path-constrained pipe sweep to `base`.
///
/// `profile_xy` is the cross-section polyline (closed) in the
/// profile-frame XY plane. `spine_3d` is the open polyline of spine
/// points along which the profile is swept. `frame_law` controls how
/// the profile is rotated as it travels (see [`FrameLaw`]).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for a too-short profile / spine
///   or non-finite coordinates.
/// - [`OcctAdvancedError::Backend`] if the sweep or feature combine
///   fails.
pub fn feat_make_pipe_with_path_constraint(
    base: &Solid,
    profile_xy: &[(f64, f64)],
    spine_3d: &[[f64; 3]],
    frame_law: FrameLaw,
    fuse_subtract: bool,
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
    if spine_3d.len() < 2 {
        return Err(OcctAdvancedError::bad_input(
            "spine_3d",
            format!("need ≥2 spine points; got {}", spine_3d.len()),
        ));
    }
    for (x, y) in profile_xy {
        if !x.is_finite() || !y.is_finite() {
            return Err(OcctAdvancedError::bad_input(
                "profile_xy",
                "profile contains a non-finite coordinate",
            ));
        }
    }
    for p in spine_3d {
        if p.iter().any(|c| !c.is_finite()) {
            return Err(OcctAdvancedError::bad_input(
                "spine_3d",
                "spine contains a non-finite coordinate",
            ));
        }
    }
    // `frame_law` selects the transport. The RMF is the implementation
    // for all three laws (see the module docs); the variant is still
    // honoured by the public contract — callers can branch on
    // `FrameLaw::is_exact_in_v1` to know whether they got the exact
    // OCCT frame.
    let _ = frame_law;

    // Lift the 2D profile into a 3D XY-plane polygon — the sweep
    // re-projects it into its own local frame at every station.
    let profile_3d: Vec<[f64; 3]> = profile_xy.iter().map(|&(x, y)| [x, y, 0.0]).collect();

    // Path-following sweep: no auxiliary spine → the only orientation
    // constraint is the path itself (the profile stays normal to the
    // spine tangent).
    let pipe = sweep_api_pipe_shell(&profile_3d, spine_3d, None, false)
        .map_err(|e| OcctAdvancedError::Backend(format!("path-constrained sweep: {e}")))?;
    feature_combine(base, &pipe, fuse_subtract)
        .map_err(|e| OcctAdvancedError::Backend(format!("pipe feature combine: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    fn square() -> Vec<(f64, f64)> {
        vec![(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    }

    #[test]
    fn rejects_short_profile() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_pipe_with_path_constraint(
            &base,
            &[(0.0, 0.0), (1.0, 0.0)],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            FrameLaw::CorrectedFrenet,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_short_spine() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_pipe_with_path_constraint(
            &base,
            &square(),
            &[[0.0, 0.0, 0.0]],
            FrameLaw::CorrectedFrenet,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_non_finite_spine() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_pipe_with_path_constraint(
            &base,
            &square(),
            &[[0.0, 0.0, 0.0], [0.0, f64::INFINITY, 1.0]],
            FrameLaw::Frenet,
            true,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn frame_law_exactness_flags() {
        assert!(FrameLaw::CorrectedFrenet.is_exact_in_v1());
        assert!(FrameLaw::DiscreteTangent.is_exact_in_v1());
        assert!(!FrameLaw::Frenet.is_exact_in_v1());
    }

    #[test]
    fn path_constrained_boss_produces_geometry() {
        // Sweep a square along a straight +Z spine and fuse onto a box.
        let base = box_solid(4.0, 4.0, 1.0)
            .unwrap()
            .translated(-2.0, -2.0, 0.0)
            .unwrap();
        let result = feat_make_pipe_with_path_constraint(
            &base,
            &square(),
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]],
            FrameLaw::CorrectedFrenet,
            true,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty(), "path-constrained boss should be non-empty");
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmax > 1.0 + 1e-6, "swept boss should rise above the base");
    }

    #[test]
    fn path_constrained_pipe_along_curved_spine() {
        // An L-shaped spine exercises the parallel-transport frame.
        let base = box_solid(8.0, 8.0, 8.0)
            .unwrap()
            .translated(-4.0, -4.0, 0.0)
            .unwrap();
        let result = feat_make_pipe_with_path_constraint(
            &base,
            &square(),
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 4.0], [4.0, 0.0, 4.0]],
            FrameLaw::DiscreteTangent,
            true,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&result, 0.3).unwrap();
        assert!(!mesh.nodes.is_empty());
    }
}
