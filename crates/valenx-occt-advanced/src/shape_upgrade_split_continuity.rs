//! Phase 150 — `ShapeUpgrade_FaceDivide` (continuity-driven) — split
//! a face at C1 discontinuities of its surface.
//!
//! ## What OCCT does
//!
//! `ShapeUpgrade_FaceDivide(face)` with the `splitter` set to
//! `ShapeUpgrade_SplitContinuity` walks the face's surface parameter
//! grid, evaluates the continuity tag at every sample (`C0` /
//! `C1` / `C2` / `G1` / …) and inserts new edges along the iso-
//! parametric lines where the continuity drops below the target.
//! Each region between cuts becomes a separate face whose surface is
//! at least the target continuity throughout.
//!
//! Required when exporting B-spline surfaces with internal C0 seams
//! (e.g. from sculpting tools that don't enforce knot multiplicity)
//! to STEP — many receiving systems silently corrupt the surface if
//! a single face crosses a C0 seam.
//!
//! ## v1 status
//!
//! Stub — needs surface continuity evaluation (compute second
//! derivatives across knot crossings, test for jump discontinuities)
//! plus face subdivision (the topology mutation truck doesn't
//! expose). Phase 150.5 lands with valenx-surface's continuity probe
//! + the face-divide topology op.

use valenx_cad::Solid;

use crate::error::OcctAdvancedError;

/// Required continuity level when subdividing faces.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum TargetContinuity {
    /// G0 (position-continuous) — splits at gaps only.
    G0,
    /// C1 (first-derivative-continuous) — the OCCT default.
    C1,
    /// C2 (second-derivative-continuous) — strict, needed for
    /// downstream curvature analysis.
    C2,
}

/// Walk `solid` and split every face whose surface has internal
/// continuity below `target`.
///
/// # Errors
///
/// Always [`OcctAdvancedError::NotYetImplemented`] in v1.
pub fn shape_upgrade_split_continuity(
    _solid: &Solid,
    _target: TargetContinuity,
) -> Result<Solid, OcctAdvancedError> {
    Err(OcctAdvancedError::not_yet("shape_upgrade_split_continuity"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn stub_with_cube_input() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_split_continuity(&cube, TargetContinuity::C1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
    }
}
