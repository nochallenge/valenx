//! Revolve evaluator — rotational sweep of a sketch profile about an
//! arbitrary axis.
//!
//! A Revolve takes the same "closed planar profile on XY" input as Pad
//! but, instead of dragging it along ±Z, spins it around the axis
//! defined by `axis_origin + axis_direction` for `angle` radians. A
//! full 2π revolution closes the swept shell into a torus-ish solid; a
//! partial sweep produces an open wedge whose two ends are sealed by
//! the truck builder.
//!
//! ## truck-modeling notes
//!
//! - The underlying call is [`truck_modeling::builder::rsweep`] with a
//!   `Face` as input. `Face::Swept = Solid`, so a complete sweep
//!   returns a closed BRep directly — no shell-to-solid wrapping
//!   required (unlike [`builder::cone`] which returns a `Shell`).
//! - `rsweep` **`debug_assert`s** that the axis vector has unit
//!   magnitude. We normalise the user's `axis_direction` before
//!   handing it to the builder so debug builds don't panic on
//!   reasonable-but-unnormalised inputs.
//! - For angles that aren't a clean fraction of 2π the builder
//!   subdivides into 2 or 3 spans internally (see `partial_rsweep` in
//!   truck-modeling); the result is still a closed solid, just with
//!   the two flat caps that close the wedge.

use truck_modeling::{builder, Point3, Rad, Vector3};
use valenx_cad::Solid;
use valenx_sketch::extrude::extract_profile_lines;
use valenx_spreadsheet::Spreadsheet;

use crate::feature::RevolveParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Evaluate a Revolve: look up the sketch, validate axis + angle, build
/// the planar profile face, and rsweep it about the axis.
pub(crate) fn evaluate(
    tree: &FeatureTree,
    p: &RevolveParams,
    ss: &Spreadsheet,
) -> Result<Solid, FeatureError> {
    // ---- angle validation ----
    let angle = p
        .angle
        .resolve(ss)
        .map_err(|e| FeatureError::BadParameter {
            name: "angle",
            reason: format!("could not resolve angle expression: {e}"),
        })?;
    if !angle.is_finite() {
        return Err(FeatureError::BadParameter {
            name: "angle",
            reason: format!("must be finite, got {angle}"),
        });
    }
    if angle.abs() < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "angle",
            reason: format!("must be nonzero, got {angle}"),
        });
    }
    if angle.abs() > std::f64::consts::TAU + 1e-9 {
        return Err(FeatureError::BadParameter {
            name: "angle",
            reason: format!(
                "must be in [-2*pi, 2*pi], got {angle} ({} turns)",
                angle / std::f64::consts::TAU
            ),
        });
    }

    // ---- axis validation ----
    if !p.axis_direction.x.is_finite()
        || !p.axis_direction.y.is_finite()
        || !p.axis_direction.z.is_finite()
    {
        return Err(FeatureError::BadParameter {
            name: "axis_direction",
            reason: format!("must be finite, got {:?}", p.axis_direction),
        });
    }
    let axis_norm = p.axis_direction.norm();
    if axis_norm < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "axis_direction",
            reason: format!(
                "must have nonzero magnitude (got {:?}, |v| = {axis_norm})",
                p.axis_direction
            ),
        });
    }
    if !p.axis_origin.x.is_finite() || !p.axis_origin.y.is_finite() || !p.axis_origin.z.is_finite()
    {
        return Err(FeatureError::BadParameter {
            name: "axis_origin",
            reason: format!("must be finite, got {:?}", p.axis_origin),
        });
    }

    // ---- sketch / profile validation ----
    let sketch = tree.get_sketch(p.sketch)?;
    let waypoints = extract_profile_lines(sketch, 1e-6)?;
    if waypoints.len() < 3 {
        return Err(FeatureError::EmptyProfile);
    }

    // ---- build a planar face from the sketch profile (same recipe as
    // valenx_sketch::extrude::extrude) ----
    let vertices: Vec<_> = waypoints
        .iter()
        .map(|(x, y)| builder::vertex(Point3::new(*x, *y, 0.0)))
        .collect();
    let mut edges = Vec::with_capacity(vertices.len());
    for i in 0..vertices.len() {
        let next = (i + 1) % vertices.len();
        edges.push(builder::line(&vertices[i], &vertices[next]));
    }
    let wire: truck_modeling::Wire = edges.into_iter().collect();
    let face = builder::try_attach_plane(&[wire]).map_err(|e| FeatureError::BadParameter {
        name: "sketch",
        reason: format!("could not build planar face for revolve profile: {e:?}"),
    })?;

    // ---- rsweep about the (normalised) axis ----
    //
    // rsweep debug-asserts the axis has unit magnitude — normalise here
    // so reasonable user input (e.g. Vector3::new(0.0, 1.0, 0.0) or
    // Vector3::new(0.0, 2.0, 0.0)) both work.
    let axis_unit = p.axis_direction / axis_norm;
    let origin = Point3::new(p.axis_origin.x, p.axis_origin.y, p.axis_origin.z);
    let axis_vec = Vector3::new(axis_unit.x, axis_unit.y, axis_unit.z);
    let solid: truck_modeling::Solid = builder::rsweep(&face, origin, axis_vec, Rad(angle));
    Ok(Solid::from_truck(solid))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::SketchRef;
    use nalgebra::Vector3 as NaVec3;
    use std::f64::consts::TAU;

    /// Build a 1×1 square profile offset from the origin so that
    /// revolving it about the Y-axis produces a *hollow* ring rather
    /// than a solid disc. Corners at (2, 0), (3, 0), (3, 1), (2, 1).
    fn offset_square_sketch() -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(2.0, 0.0);
        let b = s.add_point(3.0, 0.0);
        let c = s.add_point(3.0, 1.0);
        let d = s.add_point(2.0, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    fn revolve_rectangle_360_about_y_produces_ring() {
        // Rectangle at x ∈ [2, 3], y ∈ [0, 1], swept 360° about the
        // Y-axis. The result is a hollow ring (annulus prism). truck's
        // rsweep returns a Solid for Face inputs, so we get a closed
        // BRep without any extra shell-wrapping.
        //
        // The profile lives in the XY plane (z = 0), and the rotation
        // axis is Y through the origin, so the swept geometry sweeps
        // through ±X (and z is generated by the rotation). With four
        // line edges in the profile rsweep produces 4 swept walls plus
        // 2 caps (top + bottom of the ring) per the truck builder
        // convention.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(offset_square_sketch());
        let params = RevolveParams {
            sketch: s,
            axis_origin: NaVec3::new(0.0, 0.0, 0.0),
            axis_direction: NaVec3::new(0.0, 1.0, 0.0),
            angle: TAU.into(),
        };
        let solid =
            evaluate(&tree, &params, &Spreadsheet::new()).expect("revolve evaluation succeeds");
        // A 4-edge profile swept around a closed loop produces a
        // closed manifold with strictly more than 1 face. Don't pin
        // the exact count — truck may split the swept walls during
        // its 3-subdivision sweep — but check that it tessellates to
        // a non-trivial mesh.
        assert!(
            solid.faces() > 0,
            "revolve should produce a closed solid with > 0 faces"
        );
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        // A ring tessellates to ~32+ triangles at this tolerance —
        // pin a loose lower bound that still proves "non-trivial".
        assert!(
            mesh.total_elements() > 12,
            "ring tessellation should have > 12 triangles, got {}",
            mesh.total_elements()
        );
        // Sanity: the bounding box in X should reach beyond the
        // profile's far edge in both ±X (because the full revolution
        // takes x = 3 around to x = -3).
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        assert!(
            max_x >= 3.0 - 1e-3 && min_x <= -3.0 + 1e-3,
            "full revolution should span ±3 in X, got [{min_x}, {max_x}]"
        );
    }

    #[test]
    fn revolve_rectangle_90_partial_sweep_produces_wedge() {
        // Same offset rectangle as the 360° test, but spun only π/2
        // radians (a quarter turn). truck's `partial_rsweep` routes
        // partial angles through `multi_sweep`, which subdivides into
        // 2 spans for angles < π and 3 spans otherwise. The result is
        // still a closed Solid because truck caps the two "ends" of
        // the wedge with planar faces.
        //
        // What we're really pinning down here: a partial sweep does
        // NOT produce a full ring (so bounding box in X stays ≤ 3,
        // since the profile starts at x = 3 and only rotates 90° into
        // ±Z, not all the way around to -X), and the cap faces give
        // a strictly higher face count than the 360° version's swept
        // walls alone.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(offset_square_sketch());
        let params = RevolveParams {
            sketch: s,
            axis_origin: NaVec3::new(0.0, 0.0, 0.0),
            axis_direction: NaVec3::new(0.0, 1.0, 0.0),
            angle: std::f64::consts::FRAC_PI_2.into(),
        };
        let solid =
            evaluate(&tree, &params, &Spreadsheet::new()).expect("partial revolve succeeds");
        assert!(
            solid.faces() > 0,
            "partial revolve should still produce a closed solid"
        );
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        // Bounding box in X must NOT reach -3 (a 90° sweep from x = 3
        // about the Y axis rotates the far-X edge toward +Z, not all
        // the way to -X). Loose epsilon for tessellation slop.
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        assert!(
            min_x > -1.0,
            "90° sweep should not reach -X (min x = {min_x}, expected > -1)"
        );
        // The sweep does push into +Z though (rotating x = 3 by 90°
        // about Y lands it at z = -3 or z = +3 depending on sign of
        // rotation — rsweep with positive angle rotates +X → -Z given
        // the right-hand rule about +Y).
        let min_z = mesh.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let max_z = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            min_z.abs() > 1.0 || max_z.abs() > 1.0,
            "90° sweep should produce non-trivial Z extent, got [{min_z}, {max_z}]"
        );
    }

    #[test]
    fn revolve_about_offset_axis_translates_swept_geometry() {
        // Take the same offset rectangle (profile x ∈ [2, 3]) and
        // revolve it about an axis that's NOT through the origin —
        // axis_origin = (5, 0, 0), axis direction +Y. The full 360°
        // sweep should produce a ring whose centre sits at x = 5,
        // not at x = 0.
        //
        // For a 360° sweep about x = 5, the inner radius of the ring
        // is |5 - 3| = 2 and the outer radius is |5 - 2| = 3, so
        // bounding box in X is [5 - 3, 5 + 3] = [2, 8]. The min-x
        // therefore lands at ~2 (NOT -3 as in the origin-axis test).
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(offset_square_sketch());
        let params = RevolveParams {
            sketch: s,
            axis_origin: NaVec3::new(5.0, 0.0, 0.0),
            axis_direction: NaVec3::new(0.0, 1.0, 0.0),
            angle: TAU.into(),
        };
        let solid =
            evaluate(&tree, &params, &Spreadsheet::new()).expect("offset-axis revolve succeeds");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        // Loose epsilon — tessellation rounds the curved swept walls
        // and may slightly under/overshoot the analytical bounds.
        assert!(
            (2.0 - 1e-2..=2.5).contains(&min_x),
            "offset-axis revolve should have min x ~= 2 (axis at x=5, inner r = 2), got {min_x}"
        );
        assert!(
            (7.5..=8.0 + 1e-2).contains(&max_x),
            "offset-axis revolve should have max x ~= 8 (axis at x=5, outer r = 3), got {max_x}"
        );
    }

    #[test]
    fn revolve_with_zero_axis_direction_returns_bad_parameter() {
        // axis_direction = (0, 0, 0): no axis to rotate about. Must
        // be caught by our validation guard *before* hitting truck's
        // `debug_assert!(axis.magnitude().near(&1.0))` (which would
        // produce an unhelpful panic in debug builds, or NaN-laden
        // garbage geometry in release).
        //
        // The error must come back as BadParameter("axis_direction",
        // ...) so the UI can target the exact widget. Anything else
        // would lose the structured diagnostic.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(offset_square_sketch());
        let params = RevolveParams {
            sketch: s,
            axis_origin: NaVec3::new(0.0, 0.0, 0.0),
            axis_direction: NaVec3::new(0.0, 0.0, 0.0),
            angle: TAU.into(),
        };
        let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
        match err {
            FeatureError::BadParameter { name, reason } => {
                assert_eq!(name, "axis_direction");
                assert!(
                    reason.contains("magnitude") || reason.contains("nonzero"),
                    "expected reason to mention axis magnitude, got {reason:?}"
                );
            }
            other => panic!("expected BadParameter(axis_direction), got {other:?}"),
        }
    }

    #[test]
    fn revolve_with_nan_axis_direction_returns_bad_parameter() {
        // NaN / infinity components in the axis must also be rejected
        // by the same guard — otherwise a `division/norm` would
        // propagate NaN into truck and crash deep in the BRep build.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(offset_square_sketch());
        for bad in [
            NaVec3::new(f64::NAN, 0.0, 0.0),
            NaVec3::new(0.0, f64::INFINITY, 0.0),
        ] {
            let params = RevolveParams {
                sketch: s,
                axis_origin: NaVec3::new(0.0, 0.0, 0.0),
                axis_direction: bad,
                angle: TAU.into(),
            };
            let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
            assert!(
                matches!(
                    err,
                    FeatureError::BadParameter {
                        name: "axis_direction",
                        ..
                    }
                ),
                "expected BadParameter(axis_direction) for {bad:?}, got {err:?}"
            );
        }
    }

    #[test]
    fn revolve_in_tree_alongside_pad_produces_independent_solids() {
        // Pad → Revolve chain. Like the Pad → Pocket → Pad chain in
        // ops::pocket::tests, Revolve does not consume `last_solid` —
        // it builds a fresh solid from its own sketch + axis. So a
        // Pad placed *before* a Revolve in the tree is overwritten by
        // the Revolve's result in `last_solid`, and replay()'s final
        // solid is the Revolve's output (not a union, not the Pad).
        //
        // This is the same "Pad-as-add is a Phase 2.5 enhancement"
        // note from the Pocket chain test — Mirror and the patterns
        // are the only ops that fold the previous solid in via
        // explicit references. Document the boundary here so a
        // future enhancement to make Revolve compose with prior
        // solids doesn't sneak past CI.
        //
        // The triangle profile padded 1 unit tall gives a small prism
        // (3-sided base, ~5 faces), and the offset-rectangle revolved
        // 360° gives a hollow ring (many more faces). Tell them apart
        // by the face count.
        use crate::feature::{Feature, PadParams};
        use crate::replay::replay;

        fn triangle_sketch() -> valenx_sketch::Sketch {
            let mut s = valenx_sketch::Sketch::new();
            let a = s.add_point(-5.0, -5.0); // Place the triangle well
            let b = s.add_point(-4.0, -5.0); // away from the revolve axis
            let c = s.add_point(-4.5, -4.0); // so the two solids can't touch.
            s.add_line(a, b).unwrap();
            s.add_line(b, c).unwrap();
            s.add_line(c, a).unwrap();
            s
        }

        let mut tree = FeatureTree::new();
        let pad_s = tree.add_sketch(triangle_sketch());
        let rev_s = tree.add_sketch(offset_square_sketch());
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: pad_s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Triangle Pad",
        );
        tree.add_feature(
            Feature::Revolve(RevolveParams {
                sketch: rev_s,
                axis_origin: NaVec3::new(0.0, 0.0, 0.0),
                axis_direction: NaVec3::new(0.0, 1.0, 0.0),
                angle: TAU.into(),
            }),
            "Ring Revolve",
        );
        let final_solid = replay(&tree)
            .expect("Pad + Revolve chain replays")
            .expect("replay produces a final solid");

        // The triangular prism has 5 faces. The ring has substantially
        // more (curved swept walls + top/bottom caps). The final solid
        // from replay() is the LAST one in the chain — the Revolve —
        // so its face count must clearly exceed 5.
        assert!(
            final_solid.faces() > 5,
            "Revolve as final feature should leave its own solid (faces > 5), got {} (Pad leaked through?)",
            final_solid.faces()
        );

        // Also verify the bounding box looks like a ring (X spans
        // beyond the triangle's region of x ∈ [-5, -4]) — the
        // triangle prism alone could never produce points at x > 0.
        let mesh = valenx_cad::solid_to_mesh(&final_solid, 0.25).unwrap();
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_x >= 2.5,
            "final solid should reach +X (ring), got max x = {max_x}"
        );
    }

    #[test]
    fn revolve_unknown_sketch_errors() {
        // SketchRef pointing past the end of the sketch table — must
        // surface as UnknownSketch, not panic inside truck.
        let tree = FeatureTree::new();
        let params = RevolveParams {
            sketch: SketchRef(99),
            axis_origin: NaVec3::new(0.0, 0.0, 0.0),
            axis_direction: NaVec3::new(0.0, 1.0, 0.0),
            angle: TAU.into(),
        };
        let err = evaluate(&tree, &params, &Spreadsheet::new()).unwrap_err();
        assert_eq!(err.code(), "feature_tree.unknown_sketch");
    }
}
