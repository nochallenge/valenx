//! Chamfer evaluator — replace sharp convex edges of a target
//! feature's solid with flat bevels of constant width.
//!
//! Two-tier dispatch mirrors [`crate::ops::fillet`]: BRep path first
//! (via `valenx_fillet_brep::chamfer`), fall through to Phase 3
//! mesh-domain on soft errors. v1 always falls through pending the
//! Phase 14.5 truck-substitution work.

use std::collections::HashMap;

use valenx_cad::{solid_to_mesh, Solid};
use valenx_fillet::apply_chamfer;
use valenx_fillet_brep::error::FilletBrepError;

use crate::feature::{ChamferParams, FeatureId};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Same tessellation tolerance as the Fillet pipeline.
pub use crate::ops::fillet::{MERGE_TOLERANCE, TESS_TOLERANCE};

/// Evaluate a Chamfer: try BRep path, fall through to mesh-domain on
/// soft errors.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &ChamferParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    // ---- parameter validation ----
    if !p.distance.is_finite() || p.distance <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "distance",
            reason: format!("must be > 0 and finite, got {}", p.distance),
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
                "feature {} has not been evaluated before this chamfer (forward / self reference?)",
                p.target.0
            ),
        })?;
    let original = match target_result {
        FeatureResult::Solid(s) => s,
        FeatureResult::Suppressed => {
            return Err(FeatureError::BadParameter {
                name: "target",
                reason: format!(
                    "feature {} is suppressed; chamfer needs a live target",
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
        if !edges.is_empty() {
            match valenx_fillet_brep::chamfer::chamfer_solid_edges(brep, &edges, p.distance) {
                Ok(brep_solid) => {
                    return Ok(Solid::from_truck(brep_solid));
                }
                Err(FilletBrepError::NotPlanarFaces)
                | Err(FilletBrepError::NonConvexEdge)
                | Err(FilletBrepError::MeshBackedSolid)
                | Err(FilletBrepError::TruckSubstitutionUnavailable) => {
                    tracing::debug!(
                        target: "valenx_feature_tree::ops::chamfer",
                        "BRep chamfer fell through to mesh-domain (Phase 14 soft error)"
                    );
                }
                Err(FilletBrepError::RadiusTooLarge {
                    radius,
                    min_edge_length,
                }) => {
                    return Err(FeatureError::BadParameter {
                        name: "distance",
                        reason: format!(
                            "BRep chamfer rejected distance {radius}: too large for edge of length {min_edge_length}"
                        ),
                    });
                }
                Err(FilletBrepError::BadParameter { name, reason }) => {
                    return Err(FeatureError::BadParameter { name, reason });
                }
                Err(FilletBrepError::TruckOp(msg)) => {
                    // The BRep chamfer builds a cutter prism flush with
                    // the adjacent faces and runs the `truck_shapeops`
                    // booleans; a `TruckOp` means the boolean kernel
                    // could not resolve the coincident cutter faces for
                    // this geometry. Rather than aborting the feature,
                    // fall through to the mesh-domain pipeline — it
                    // builds the same flat-bevel shape by strip
                    // insertion. This matches `ops::fillet`'s
                    // robustness-minded handling of the same error.
                    tracing::debug!(
                        target: "valenx_feature_tree::ops::chamfer",
                        msg = %msg,
                        "BRep chamfer boolean failed; falling through to mesh-domain"
                    );
                }
            }
        }
    }

    // ---- Mesh-domain fallback (Phase 3 path) ----
    let raw = solid_to_mesh(original, TESS_TOLERANCE)?;
    let mesh = valenx_mesh::boolean::merge_coincident_nodes(&raw, MERGE_TOLERANCE);
    let threshold_rad = p.threshold_deg.to_radians();
    let chamfered = apply_chamfer(&mesh, p.distance, threshold_rad)?;
    Ok(Solid::from_mesh(chamfered))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams};
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
    fn pad_then_chamfer_produces_a_solid() {
        // Same shape as the Fillet test: under v1 BRep falls through
        // to mesh-domain; once Phase 14.5 ships the result is BRep.
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
            Feature::Chamfer(ChamferParams {
                target: pad_id,
                distance: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Chamfer 1mm",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        let mesh = solid_to_mesh(&solid, 0.5).expect("round-trip");
        assert!(
            mesh.total_elements() > 12,
            "chamfered unit cube should have more triangles than the original 12, got {}",
            mesh.total_elements()
        );
    }

    #[test]
    fn chamfer_rejects_zero_distance() {
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
            Feature::Chamfer(ChamferParams {
                target: pad_id,
                distance: 0.0,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Bad Chamfer",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "distance",
                ..
            }
        ));
    }

    #[test]
    fn chamfer_targeting_suppressed_feature_returns_bad_parameter() {
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
            Feature::Chamfer(ChamferParams {
                target: pad_id,
                distance: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Dangling Chamfer",
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
}
