//! Pad evaluator — extrudes a sketch profile into a solid.
//!
//! The Pad operation is the bread-and-butter of part design: a closed 2D
//! profile drawn on the XY working plane gets swept along ±Z by `depth`
//! units to produce a prism. Internally we lean on Phase 1's
//! [`valenx_sketch::Sketch::extrude`] (which itself wraps truck-modeling's
//! `tsweep`) so the brep-construction pipeline is shared end-to-end with
//! the standalone Sketcher panel.

use valenx_cad::Solid;
use valenx_sketch::extrude::extract_profile_lines;
use valenx_spreadsheet::Spreadsheet;

use crate::feature::PadParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Evaluate a Pad: look up the sketch, validate the profile + depth,
/// and extrude.
///
/// `direction_positive == true` extrudes along +Z; `false` flips the
/// sweep vector by negating the depth before handing it to truck.
///
/// `ss` is the spreadsheet context used to resolve any
/// [`valenx_feature_tree::feature::Value::Expression`] depth value. If
/// the call site has no live spreadsheet, pass `&Spreadsheet::new()`
/// — literal depths still pass through.
pub(crate) fn evaluate(
    tree: &FeatureTree,
    p: &PadParams,
    ss: &Spreadsheet,
) -> Result<Solid, FeatureError> {
    let depth = p
        .depth
        .resolve(ss)
        .map_err(|e| FeatureError::BadParameter {
            name: "depth",
            reason: format!("could not resolve depth expression: {e}"),
        })?;
    if !depth.is_finite() || depth.abs() < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "depth",
            reason: format!("must be nonzero and finite, got {depth}"),
        });
    }
    let sketch = tree.get_sketch(p.sketch)?;
    let waypoints = extract_profile_lines(sketch, 1e-6)?;
    if waypoints.len() < 3 {
        return Err(FeatureError::EmptyProfile);
    }
    let signed_depth = if p.direction_positive { depth } else { -depth };
    // Delegate to the sketch's own extrude (reuses the truck-modeling
    // pipeline that Phase 1 already validated).
    let solid = sketch.extrude(signed_depth)?;
    Ok(solid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::SketchRef;

    /// Build a sketch containing a closed triangular loop:
    /// (0,0) → (1,0) → (0.5,1) → (0,0).
    fn triangle_sketch() -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(0.5, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, a).unwrap();
        s
    }

    #[test]
    fn pad_triangle_profile_produces_solid() {
        // Standalone Pad evaluate(): wrap the triangle in a FeatureTree
        // and confirm the extrude returns a valid Solid that tessellates
        // to a non-empty mesh.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(triangle_sketch());
        let params = PadParams {
            sketch: s,
            depth: 1.5.into(),
            direction_positive: true,
        };
        let solid = evaluate(&tree, &params, &Spreadsheet::new()).expect("pad evaluation succeeds");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.5).unwrap();
        assert!(!mesh.nodes.is_empty(), "tessellation produced nodes");
        assert!(mesh.total_elements() > 0, "tessellation produced triangles");
    }

    #[test]
    fn pad_unknown_sketch_errors() {
        // SketchRef pointing past the end of the sketch table — should
        // surface as UnknownSketch, not panic.
        let tree = FeatureTree::new();
        let params = PadParams {
            sketch: SketchRef(99),
            depth: 1.0.into(),
            direction_positive: true,
        };
        let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
        assert_eq!(err.code(), "feature_tree.unknown_sketch");
    }

    #[test]
    fn pad_with_negative_depth_produces_solid_below_plane() {
        // direction_positive = false should flip the sweep vector to -Z.
        // The triangle profile lives at z = 0, so the resulting prism's
        // bounding box should sit at -depth..=0 along Z (not 0..=+depth).
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(triangle_sketch());
        let params = PadParams {
            sketch: s,
            depth: 2.0.into(),
            direction_positive: false,
        };
        let solid = evaluate(&tree, &params, &Spreadsheet::new()).expect("pad evaluation succeeds");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.5).unwrap();
        assert!(!mesh.nodes.is_empty(), "tessellation produced nodes");
        let min_z = mesh.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let max_z = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        // Allow a small slop for tessellation; the prism should span
        // approximately -2.0..=0.0.
        assert!(
            min_z <= -2.0 + 1e-6,
            "min z should reach near -2.0, got {min_z}"
        );
        assert!(
            (-1e-6..=1e-6).contains(&max_z),
            "max z should be approximately 0.0, got {max_z}"
        );
    }

    #[test]
    fn pad_with_empty_sketch_returns_empty_profile() {
        // Brand-new sketch has no entities, so extract_profile_lines
        // returns an empty waypoint list; we must reject with the
        // dedicated EmptyProfile variant (not a generic CAD error).
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(valenx_sketch::Sketch::new());
        let params = PadParams {
            sketch: s,
            depth: 1.0.into(),
            direction_positive: true,
        };
        let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
        assert!(
            matches!(err, FeatureError::EmptyProfile),
            "expected EmptyProfile, got {err:?}"
        );
        assert_eq!(err.code(), "feature_tree.empty_profile");
    }

    #[test]
    fn pad_with_zero_depth_returns_bad_parameter() {
        // depth = 0 fails the early guard before we ever touch the
        // sketch, so the profile being valid is irrelevant. We also
        // check that the parameter name in the diagnostic is "depth".
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(triangle_sketch());
        let params = PadParams {
            sketch: s,
            depth: 0.0.into(),
            direction_positive: true,
        };
        let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
        match err {
            FeatureError::BadParameter { name, .. } => {
                assert_eq!(name, "depth");
            }
            other => panic!("expected BadParameter, got {other:?}"),
        }
    }

    #[test]
    fn pad_with_non_finite_depth_returns_bad_parameter() {
        // NaN / infinity must also be rejected by the same guard;
        // this is a belt-and-braces check that the finite-check fires.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(triangle_sketch());
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let params = PadParams {
                sketch: s,
                depth: bad.into(),
                direction_positive: true,
            };
            let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
            assert!(
                matches!(err, FeatureError::BadParameter { name: "depth", .. }),
                "expected BadParameter(depth) for {bad}, got {err:?}"
            );
        }
    }

    #[test]
    fn pad_via_replay_matches_direct_sketch_extrude() {
        // Sanity: routing the Pad through `FeatureTree::replay` must
        // produce a solid with the same topology as calling
        // `Sketch::extrude` directly on the same sketch + depth.
        // The two paths share the underlying truck pipeline, so this
        // is really checking that the replay shell isn't dropping or
        // mutating anything on the way through.
        use crate::feature::Feature;
        use crate::replay::replay;

        let sketch = triangle_sketch();
        let depth = 2.5;
        let direct = sketch.extrude(depth).expect("direct extrude succeeds");

        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(sketch);
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: depth.into(),
                direction_positive: true,
            }),
            "Pad",
        );
        let via_replay = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");

        // Same topology counts (faces / edges / vertices) is the
        // cheapest invariant we can check across two independent brep
        // builds; if any of these diverge the extrude paths have
        // drifted apart and the test should fail loudly.
        assert_eq!(direct.faces(), via_replay.faces(), "face count diverged");
        assert_eq!(direct.edges(), via_replay.edges(), "edge count diverged");
        assert_eq!(
            direct.vertices(),
            via_replay.vertices(),
            "vertex count diverged"
        );

        // And a coarse tessellation should yield the same number of
        // triangles — flexible enough that float-noise from two
        // separate truck builds doesn't matter.
        let mesh_a = valenx_cad::solid_to_mesh(&direct, 0.5).unwrap();
        let mesh_b = valenx_cad::solid_to_mesh(&via_replay, 0.5).unwrap();
        assert_eq!(
            mesh_a.total_elements(),
            mesh_b.total_elements(),
            "tessellation triangle count diverged"
        );
    }
}
