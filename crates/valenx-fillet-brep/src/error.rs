//! Error taxonomy for the BRep fillet/chamfer pipeline.
//!
//! Same shape as `valenx_fillet::FilletError` (`code()` +
//! [`ErrorCategory`]) so the feature-tree's `FeatureError` can route
//! both layers uniformly. Two of the variants —
//! [`FilletBrepError::NotPlanarFaces`] and
//! [`FilletBrepError::NonConvexEdge`] — are *soft* errors: the
//! upstream BRep is valid, the BRep fillet just can't handle it in
//! v1. Callers should treat these as a signal to fall through to the
//! mesh-domain pipeline rather than surfacing the error to the user.

use thiserror::Error;

/// Errors raised by `valenx-fillet-brep`'s fillet / chamfer functions.
#[derive(Debug, Error)]
pub enum FilletBrepError {
    /// The two faces adjacent to the target edge are not both planar.
    /// v1 only fillets between two planar faces; non-planar adjacency
    /// is a Phase 14.5+ problem (offset-surface intersection).
    ///
    /// Soft error: callers should fall through to mesh-domain fillet.
    #[error("BRep fillet v1 requires both adjacent faces to be planar")]
    NotPlanarFaces,

    /// The target edge is concave (the two adjacent faces meet with
    /// the convex side pointing inward). v1 only handles convex
    /// edges — concave-edge fillet ("fillet a re-entrant corner")
    /// requires a sphere blend that's Phase 14.5+.
    ///
    /// Soft error: callers should fall through to mesh-domain fillet
    /// (which also handles only convex edges but flags them as
    /// such in its own classification step).
    #[error("BRep fillet v1 requires a convex edge")]
    NonConvexEdge,

    /// Requested radius is too large for the edge — the offset
    /// tangent point would land beyond the edge's other endpoint,
    /// producing a self-intersecting fillet surface.
    #[error("fillet radius {radius} too large for edge of length {min_edge_length}")]
    RadiusTooLarge {
        /// Radius the caller asked for.
        radius: f64,
        /// Length of the shortest edge that violated the bound.
        min_edge_length: f64,
    },

    /// Caller passed a [`valenx_cad::Solid::Mesh`] (no BRep topology
    /// to operate on). The BRep fillet pipeline cannot fix a
    /// mesh-backed solid — fall through to the mesh-domain pipeline
    /// (which is what produced the mesh in the first place).
    ///
    /// Soft error: callers should fall through to mesh-domain.
    #[error("BRep fillet requires a BRep solid; input was mesh-backed")]
    MeshBackedSolid,

    /// truck-modeling refused a construction step (degenerate sweep,
    /// non-manifold result, NURBS evaluation failure, etc.). String-
    /// wrapped because `truck_modeling::errors::Error` mixes
    /// I/O variants we never produce; carrying the full type would
    /// leak STEP/IGES symbols into a crate that doesn't depend on
    /// those.
    #[error("truck op failed: {0}")]
    TruckOp(String),

    /// truck-modeling 0.6 does not expose a generic face-trim-and-
    /// substitute primitive — the BRep fillet surgery would have to
    /// reconstruct the entire solid topology by hand for each edge.
    /// v1 surfaces this honestly as a soft error so the dispatcher
    /// falls through to the mesh-domain pipeline, where the same
    /// fillet shape is achievable via per-triangle strip insertion.
    ///
    /// Full BRep substitution waits for a truck release that ships
    /// `truck_shapeops::trim_face` (or equivalent); see Phase 14.5+.
    ///
    /// Soft error: callers should fall through to mesh-domain.
    #[error(
        "BRep face-trim-and-substitute is not exposed by truck-modeling 0.6 \
         (Phase 14.5+); v1 falls through to mesh-domain"
    )]
    TruckSubstitutionUnavailable,

    /// User-supplied parameter was out of range or nonsense.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"radius"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },
}

/// Coarse category an LLM / UI can branch on to decide who to escalate
/// the error to.
///
/// Mirrors `valenx_fillet::ErrorCategory` and
/// `valenx_feature_tree::ErrorCategory` so the three layers can share
/// routing logic.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Input geometry can't be handled by v1 — try a different op
    /// (typically mesh-domain fillet) or fix the model.
    Input,
    /// User-tunable knob out of range (fix the radius).
    Config,
    /// Transient failure inside truck or another runtime dependency.
    Runtime,
    /// Bug in valenx-fillet-brep or a truck op behaving unexpectedly.
    Internal,
}

impl FilletBrepError {
    /// Stable kebab-cased identifier; never changes across versions.
    /// New variants get new codes.
    pub fn code(&self) -> &'static str {
        match self {
            FilletBrepError::NotPlanarFaces => "fillet_brep.not_planar_faces",
            FilletBrepError::NonConvexEdge => "fillet_brep.non_convex_edge",
            FilletBrepError::RadiusTooLarge { .. } => "fillet_brep.radius_too_large",
            FilletBrepError::MeshBackedSolid => "fillet_brep.mesh_backed_solid",
            FilletBrepError::TruckOp(_) => "fillet_brep.truck_op",
            FilletBrepError::TruckSubstitutionUnavailable => {
                "fillet_brep.truck_substitution_unavailable"
            }
            FilletBrepError::BadParameter { .. } => "fillet_brep.bad_parameter",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FilletBrepError::NotPlanarFaces
            | FilletBrepError::NonConvexEdge
            | FilletBrepError::MeshBackedSolid
            | FilletBrepError::TruckSubstitutionUnavailable => ErrorCategory::Input,
            FilletBrepError::RadiusTooLarge { .. } | FilletBrepError::BadParameter { .. } => {
                ErrorCategory::Config
            }
            FilletBrepError::TruckOp(_) => ErrorCategory::Internal,
        }
    }

    /// True if this is a *soft* error — the caller should fall through
    /// to the mesh-domain pipeline rather than surfacing to the user.
    /// Hard errors (RadiusTooLarge, BadParameter, TruckOp) are
    /// genuine problems and should be reported up.
    pub fn is_soft(&self) -> bool {
        matches!(
            self,
            FilletBrepError::NotPlanarFaces
                | FilletBrepError::NonConvexEdge
                | FilletBrepError::MeshBackedSolid
                | FilletBrepError::TruckSubstitutionUnavailable
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(FilletBrepError, &'static str, ErrorCategory)> = vec![
            (
                FilletBrepError::NotPlanarFaces,
                "fillet_brep.not_planar_faces",
                ErrorCategory::Input,
            ),
            (
                FilletBrepError::NonConvexEdge,
                "fillet_brep.non_convex_edge",
                ErrorCategory::Input,
            ),
            (
                FilletBrepError::RadiusTooLarge {
                    radius: 1.0,
                    min_edge_length: 0.5,
                },
                "fillet_brep.radius_too_large",
                ErrorCategory::Config,
            ),
            (
                FilletBrepError::MeshBackedSolid,
                "fillet_brep.mesh_backed_solid",
                ErrorCategory::Input,
            ),
            (
                FilletBrepError::TruckOp("sweep failed".into()),
                "fillet_brep.truck_op",
                ErrorCategory::Internal,
            ),
            (
                FilletBrepError::TruckSubstitutionUnavailable,
                "fillet_brep.truck_substitution_unavailable",
                ErrorCategory::Input,
            ),
            (
                FilletBrepError::BadParameter {
                    name: "radius",
                    reason: "must be > 0".into(),
                },
                "fillet_brep.bad_parameter",
                ErrorCategory::Config,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "code mismatch for {err:?}");
            assert_eq!(
                err.category(),
                expected_cat,
                "category mismatch for {err:?}"
            );
        }
    }

    #[test]
    fn soft_vs_hard_classification() {
        // Soft errors: caller should fall through to mesh-domain.
        assert!(FilletBrepError::NotPlanarFaces.is_soft());
        assert!(FilletBrepError::NonConvexEdge.is_soft());
        assert!(FilletBrepError::MeshBackedSolid.is_soft());
        assert!(FilletBrepError::TruckSubstitutionUnavailable.is_soft());
        // Hard errors: surface to user.
        assert!(!FilletBrepError::RadiusTooLarge {
            radius: 1.0,
            min_edge_length: 0.5
        }
        .is_soft());
        assert!(!FilletBrepError::TruckOp("x".into()).is_soft());
        assert!(!FilletBrepError::BadParameter {
            name: "x",
            reason: "y".into(),
        }
        .is_soft());
    }

    #[test]
    fn display_includes_useful_context() {
        // RadiusTooLarge should mention both numbers so users can fix
        // the input without reading the source.
        let err = FilletBrepError::RadiusTooLarge {
            radius: 2.5,
            min_edge_length: 1.0,
        };
        let msg = err.to_string();
        assert!(msg.contains("2.5"));
        assert!(msg.contains('1'));
    }
}
