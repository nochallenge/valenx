//! Fillet evaluator — round sharp convex edges of a target feature's
//! solid.
//!
//! # Two-tier dispatch (Phase 14 / 14.5)
//!
//! 1. **BRep path.** If the target is a [`Solid::Brep`], try
//!    `valenx_fillet_brep::fillet::fillet_solid_edges`. As of Phase
//!    14.5 this is a **real BRep fillet** — it builds the cutter +
//!    fillet-bar prisms and runs the `truck_shapeops` booleans. When
//!    it succeeds the result is a true [`Solid::Brep`] with a genuine
//!    circular-arc fillet face that round-trips through STEP and
//!    survives downstream booleans.
//! 2. **Mesh fallback.** When the BRep path returns a soft error
//!    (non-planar adjacency, non-convex edge, mesh-backed input), or
//!    the `truck_shapeops` boolean cannot resolve the flush cutter
//!    faces (`TruckOp`), or the target is already [`Solid::Mesh`],
//!    the evaluator falls through to the Phase 3 mesh-domain
//!    pipeline: tessellate → weld → mesh-fillet → wrap. The
//!    mesh-domain pipeline produces the same fillet shape.
//!
//! Genuinely-bad parameters from the BRep path (`RadiusTooLarge`,
//! `BadParameter`) bypass the fallback and surface to the user
//! verbatim — a malformed radius shouldn't silently swallow into
//! mesh-domain with a different radius.
//!
//! # Edge selection
//!
//! - When [`FilletParams::edge_indices`] is `Some`, those indices are
//!   interpreted against `valenx_fillet_brep::bridge::unique_edges`
//!   for BRep targets. Mesh-backed targets ignore the indices (mesh-
//!   domain has no concept of BRep edge IDs).
//! - When `edge_indices` is `None`, BRep targets auto-select by
//!   dihedral angle (via `valenx_fillet_brep::bridge::
//!   edges_for_threshold`); mesh targets use the Phase 3 mesh-domain
//!   angle classifier.
//!
//! # Scope of the BRep path
//!
//! The BRep fillet handles convex straight edges between two planar
//! faces, constant radius (see `valenx-fillet-brep`); when 3 such
//! filleted edges meet at an **orthogonal convex corner** (the corner
//! of a box) `fillet_solid_edges` also blends that corner with a real
//! rolling-ball spherical corner. General N-edge / non-orthogonal /
//! concave corners, curved adjacent faces, and concave edges remain
//! mesh-domain or deferred. Because the cutter trims flush with the
//! adjacent faces, some geometry trips the boolean kernel's
//! coincident-face handling; that case falls through to mesh-domain
//! rather than failing the feature (a corner blend that trips it
//! falls back to the independent per-edge fillets).

use std::collections::HashMap;

use valenx_cad::{solid_to_mesh, Solid};
use valenx_fillet::apply_fillet;
use valenx_fillet_brep::error::FilletBrepError;

use crate::feature::{FeatureId, FilletParams};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Tessellation tolerance for the BRep → mesh step. Matches the
/// default the mesh-toolbox uses for viewport tessellation
/// ([`valenx_cad::DEFAULT_TESS_TOLERANCE`] = 0.5).
pub const TESS_TOLERANCE: f64 = valenx_cad::DEFAULT_TESS_TOLERANCE;

/// Spatial tolerance for welding tessellation output vertices before
/// fillet/chamfer. See Phase 3 docs.
pub const MERGE_TOLERANCE: f64 = 1e-4;

/// Evaluate a Fillet: try the BRep path, fall through to mesh-domain
/// on soft errors.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &FilletParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    // ---- parameter validation ----
    if !p.radius.is_finite() || p.radius <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "radius",
            reason: format!("must be > 0 and finite, got {}", p.radius),
        });
    }
    if !p.threshold_deg.is_finite() || p.threshold_deg < 0.0 {
        return Err(FeatureError::BadParameter {
            name: "threshold_deg",
            reason: format!("must be >= 0 and finite, got {}", p.threshold_deg),
        });
    }

    // ---- target lookup ----
    let target_result = prior
        .get(&p.target)
        .ok_or_else(|| FeatureError::BadParameter {
            name: "target",
            reason: format!(
                "feature {} has not been evaluated before this fillet (forward / self reference?)",
                p.target.0
            ),
        })?;
    let original = match target_result {
        FeatureResult::Solid(s) => s,
        FeatureResult::Suppressed => {
            return Err(FeatureError::BadParameter {
                name: "target",
                reason: format!(
                    "feature {} is suppressed; fillet needs a live target",
                    p.target.0
                ),
            });
        }
    };

    // ---- BRep path (try first) ----
    if let Solid::Brep(brep) = original {
        let threshold_rad = p.threshold_deg.to_radians();
        let edges = match &p.edge_indices {
            Some(idxs) => valenx_fillet_brep::bridge::edges_for_indices(brep, idxs),
            None => valenx_fillet_brep::bridge::edges_for_threshold(brep, threshold_rad),
        };
        // Only attempt BRep if we actually have edges to fillet —
        // an empty list is a soft signal to fall through to mesh
        // (the mesh-domain pipeline has its own classifier).
        if !edges.is_empty() {
            match valenx_fillet_brep::fillet::fillet_solid_edges(brep, &edges, p.radius) {
                Ok(brep_solid) => {
                    return Ok(Solid::from_truck(brep_solid));
                }
                Err(FilletBrepError::NotPlanarFaces)
                | Err(FilletBrepError::NonConvexEdge)
                | Err(FilletBrepError::MeshBackedSolid)
                | Err(FilletBrepError::TruckSubstitutionUnavailable) => {
                    // Soft errors — log and fall through to the
                    // mesh-domain pipeline, which has its own
                    // classifier and produces the same fillet shape.
                    tracing::debug!(
                        target: "valenx_feature_tree::ops::fillet",
                        "BRep fillet fell through to mesh-domain (unsupported geometry)"
                    );
                }
                Err(FilletBrepError::RadiusTooLarge {
                    radius,
                    min_edge_length,
                }) => {
                    return Err(FeatureError::BadParameter {
                        name: "radius",
                        reason: format!(
                            "BRep fillet rejected radius {radius}: too large for edge of length {min_edge_length}"
                        ),
                    });
                }
                Err(FilletBrepError::BadParameter { name, reason }) => {
                    return Err(FeatureError::BadParameter { name, reason });
                }
                Err(FilletBrepError::TruckOp(msg)) => {
                    // The real BRep fillet (Phase 14.5) builds the
                    // cutter + fillet-bar prisms and runs the
                    // `truck_shapeops` booleans. A `TruckOp` here
                    // means the boolean kernel could not resolve the
                    // flush cutter faces for this geometry — the BRep
                    // fillet genuinely failed. Rather than aborting
                    // the feature, fall through to the mesh-domain
                    // pipeline, which builds the same fillet shape by
                    // strip insertion.
                    tracing::debug!(
                        target: "valenx_feature_tree::ops::fillet",
                        msg = %msg,
                        "BRep fillet boolean failed; falling through to mesh-domain"
                    );
                }
            }
        }
    }

    // ---- Mesh-domain fallback (Phase 3 path) ----
    let raw = solid_to_mesh(original, TESS_TOLERANCE)?;
    let mesh = valenx_mesh::boolean::merge_coincident_nodes(&raw, MERGE_TOLERANCE);
    let threshold_rad = p.threshold_deg.to_radians();
    let filleted = apply_fillet(&mesh, p.radius, threshold_rad)?;
    Ok(Solid::from_mesh(filleted))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams, SketchRef};
    use crate::replay::replay;

    /// Build a unit-square sketch.
    fn square_sketch() -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(1.0, 1.0);
        let d = s.add_point(0.0, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    fn pad_then_fillet_produces_a_solid() {
        // Pad a unit square 1 unit up → unit cube. Fillet with radius
        // 0.1 at 45°. The dispatcher tries the real BRep fillet first;
        // it either succeeds (a true Brep result) or — if the
        // coincident-face boolean trips the kernel — falls through to
        // the mesh-domain pipeline. Either way the result is a valid
        // filleted solid with more triangles than the bare cube.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: pad_id,
                radius: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Fillet 1mm",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        // Accept either backend (real BRep result, or mesh-domain
        // fall-through): both produce a valid filleted solid.
        let mesh = solid_to_mesh(&solid, 0.5).expect("round-trip");
        assert!(
            mesh.total_elements() > 12,
            "filleted unit cube should have more triangles than the original 12, got {}",
            mesh.total_elements()
        );
    }

    #[test]
    fn fillet_targeting_suppressed_feature_returns_bad_parameter() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Suppressed Pad",
        );
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: pad_id,
                radius: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Dangling Fillet",
        );
        tree.set_suppressed(pad_id, true).unwrap();
        let err = replay(&tree).unwrap_err();
        match err {
            FeatureError::BadParameter { name, reason } => {
                assert_eq!(name, "target");
                assert!(reason.contains("suppressed"));
            }
            other => panic!("expected BadParameter(target), got {other:?}"),
        }
    }

    #[test]
    fn fillet_targeting_forward_reference_returns_bad_parameter() {
        let mut tree = FeatureTree::new();
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: FeatureId(99),
                radius: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Forward Fillet",
        );
        let s = tree.add_sketch(square_sketch());
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Pad After Fillet",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter { name: "target", .. }
        ));
    }

    #[test]
    fn fillet_rejects_zero_radius() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Pad",
        );
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: pad_id,
                radius: 0.0,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Bad Fillet",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter { name: "radius", .. }
        ));
    }

    #[test]
    fn fillet_rejects_negative_threshold() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Pad",
        );
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: pad_id,
                radius: 0.1,
                threshold_deg: -10.0,
                edge_indices: None,
            }),
            "Bad Fillet",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "threshold_deg",
                ..
            }
        ));
    }

    #[test]
    fn fillet_with_unknown_target_returns_bad_parameter() {
        let mut tree = FeatureTree::new();
        let _ = tree.add_sketch(valenx_sketch::Sketch::new());
        let _ = SketchRef(0);
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: FeatureId(0),
                radius: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Standalone Fillet",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter { name: "target", .. }
        ));
    }

    #[test]
    fn pad_then_fillet_with_explicit_edge_indices_replays() {
        // edge_indices=Some(...) targets the BRep path explicitly.
        // Under v1 it still falls through (substitution unavailable);
        // the test just confirms the dispatch wiring doesn't blow
        // up when indices are provided. Phase 14.5 upgrades this to
        // a topology assertion on the resulting BRep.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: pad_id,
                radius: 0.1,
                threshold_deg: 45.0,
                edge_indices: Some(vec![0, 1, 2, 3]),
            }),
            "Fillet 4 edges",
        );
        let solid = replay(&tree).unwrap().unwrap();
        // Should always be representable as a mesh (BRep or fallback).
        let mesh = solid_to_mesh(&solid, 0.5).unwrap();
        assert!(mesh.total_elements() > 0);
    }
}
