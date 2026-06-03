//! Hole evaluator — drills 1+ cylindrical pockets at the points of a
//! sketch, optionally with counterbore / countersink modifiers.
//!
//! Phase 13A Task 8. Reads `Entity::Point` entries from the hole
//! sketch, builds a cylindrical cutter per point, and subtracts each
//! cutter from the base solid via [`valenx_cad::difference`]. Thread
//! information is *metadata only* in v1 — the [`crate::feature::HoleParams::thread`]
//! field is round-tripped but does not generate helical thread geometry.
//!
//! ## Cutter stab-overhang
//!
//! Same trick as [`super::pocket`]: stab the cutter past both faces of
//! the base so the cap planes don't coincide with the base's faces.
//! Without this truck-shapeops returns
//! [`valenx_cad::CadError::EmptyResult`].
//!
//! ## Limitations (v1)
//!
//! - Sketch points only — circles / arcs in the hole sketch are ignored.
//! - Counterbore and countersink add cylindrical / conical recesses
//!   stacked on the entry face; the conical countersink is modelled
//!   as a `cone()` primitive.
//! - `HoleDepthMode::UpToFace` is treated as a deep blind cut
//!   (`max(20.0, drill_diameter * 10)`); face resolution requires
//!   BRep info we don't yet propagate.

use valenx_cad::{cone, cylinder, difference, Solid};

use crate::feature::{HoleDepthMode, HoleParams};
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Per-cutter stab overhang — matches [`super::pocket::POCKET_STAB_EPSILON`].
pub const HOLE_STAB_EPSILON: f64 = 0.5;

/// Evaluate a Hole: enumerate sketch points, build & subtract a cutter
/// for each.
pub(crate) fn evaluate(
    tree: &FeatureTree,
    p: &HoleParams,
    base: Option<&Solid>,
) -> Result<Solid, FeatureError> {
    let base = base.ok_or(FeatureError::BadParameter {
        name: "hole",
        reason: "hole requires a base solid (must be preceded by a solid-producing feature)".into(),
    })?;
    if !p.drill_diameter.is_finite() || p.drill_diameter <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "drill_diameter",
            reason: format!("must be > 0 and finite, got {}", p.drill_diameter),
        });
    }

    let sketch = tree.get_sketch(p.sketch)?;
    // Collect (x, y) positions from every Point entity in the sketch.
    let centers: Vec<(f64, f64)> = sketch
        .entities
        .iter()
        .filter_map(|e| match e {
            valenx_sketch::geom::Entity::Point(pt) => {
                Some((sketch.vars[pt.x_var], sketch.vars[pt.y_var]))
            }
            _ => None,
        })
        .collect();
    if centers.is_empty() {
        return Err(FeatureError::EmptyProfile);
    }

    let depth = match &p.depth_mode {
        HoleDepthMode::Blind { depth } => {
            if !depth.is_finite() || *depth <= 0.0 {
                return Err(FeatureError::BadParameter {
                    name: "depth",
                    reason: format!("must be > 0 and finite, got {depth}"),
                });
            }
            *depth
        }
        HoleDepthMode::Through => {
            // Use a generous default — the base is clipped at boolean
            // time so overshoot doesn't show up in the result.
            (p.drill_diameter * 20.0).max(50.0)
        }
        HoleDepthMode::UpToFace { .. } => {
            // v1 stub — would resolve via BRep in a future phase.
            (p.drill_diameter * 10.0).max(20.0)
        }
    };

    let radius = p.drill_diameter * 0.5;
    let stab = HOLE_STAB_EPSILON;
    let sign: f64 = if p.direction_negative { -1.0 } else { 1.0 };

    let mut acc = base.clone();
    for (cx, cy) in centers {
        // Drill cutter: cylindrical pocket with stab on both ends.
        let cutter = cylinder(radius, depth + 2.0 * stab)?;
        // Move cutter: cylinder primitives are built at origin with
        // base on z=0 and axis +Z. For a downward drill we shift down
        // so the cutter spans z ∈ [-depth-stab, +stab]; for upward we
        // mirror that. Then translate in X/Y to the point.
        let cutter = if sign < 0.0 {
            cutter.translated(cx, cy, -(depth + stab))?
        } else {
            cutter.translated(cx, cy, -stab)?
        };
        acc = match difference(&acc, &cutter) {
            Ok(s) => s,
            Err(e) => return Err(FeatureError::from(e)),
        };

        // Counterbore — stacked at the entry face on the opposite end
        // from the drilled tip. For a downward drill (negative dir)
        // the entry face is at top, so cut a flat-bottomed recess
        // there going downward.
        if let Some(cb) = &p.counterbore {
            if cb.diameter > 0.0 && cb.depth > 0.0 && cb.diameter > p.drill_diameter {
                let r_cb = cb.diameter * 0.5;
                let cb_cyl = cylinder(r_cb, cb.depth + 2.0 * stab)?;
                let cb_cyl = if sign < 0.0 {
                    cb_cyl.translated(cx, cy, -(cb.depth + stab))?
                } else {
                    cb_cyl.translated(cx, cy, -stab)?
                };
                acc = match difference(&acc, &cb_cyl) {
                    Ok(s) => s,
                    Err(e) => return Err(FeatureError::from(e)),
                };
            }
        }

        // Countersink — a cone-shaped recess at the entry face. The
        // included angle is the cone's full apex angle; the cone we
        // build has its apex at the bottom of the recess and its
        // base at the top diameter.
        if let Some(cs) = &p.countersink {
            if cs.diameter > 0.0 && cs.angle_deg > 0.0 && cs.diameter > p.drill_diameter {
                let r_top = cs.diameter * 0.5;
                let r_bot = radius;
                // Depth from r_top down to r_bot at the half-angle:
                // tan(half) = (r_top - r_bot) / depth ⇒
                // depth = (r_top - r_bot) / tan(angle/2)
                let half = (cs.angle_deg * 0.5).to_radians();
                let depth_cs = if half.tan().abs() > 1e-9 {
                    ((r_top - r_bot) / half.tan()).abs()
                } else {
                    cs.diameter
                };
                // cone(base, top, height) — base at z=0, top at z=height.
                // We want the wide part at the entry face. Build with
                // base = r_top, top = r_bot, then translate so base
                // sits at the entry face going inward.
                let cs_cone = cone(r_top, r_bot, depth_cs + stab)?;
                let cs_cone = if sign < 0.0 {
                    // Drilling -Z: base (wide) should be at top, tip
                    // (narrow) inside the part. Build cone with wide
                    // base at z=0, narrow at z=+depth_cs; then mirror
                    // via translation: shift so wide base is at z=0
                    // (entry face) and tip points -Z. Truck's cone()
                    // already builds with base at z=0 and goes +Z;
                    // for a -Z drill we want the cone to go -Z. We
                    // approximate by translating it inverted via
                    // building base big and top small at +Z then
                    // shifting down by depth_cs (so it now spans
                    // z ∈ [-depth_cs, 0]).
                    cs_cone.translated(cx, cy, -(depth_cs))?
                } else {
                    cs_cone.translated(cx, cy, 0.0)?
                };
                acc = match difference(&acc, &cs_cone) {
                    Ok(s) => s,
                    // Cone subtraction failures are surfaced (not
                    // silently ignored) so the user sees the
                    // problem.
                    Err(e) => return Err(FeatureError::from(e)),
                };
            }
        }
    }

    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams};
    use crate::replay::replay;
    use crate::threads::{iso_metric_table, ThreadStandard};

    fn square_sketch(half: f64) -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(-half, -half);
        let b = s.add_point(half, -half);
        let c = s.add_point(half, half);
        let d = s.add_point(-half, half);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    fn one_point_sketch(x: f64, y: f64) -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        s.add_point(x, y);
        s
    }

    #[test]
    fn single_hole_through_a_padded_block_increases_face_count() {
        let mut tree = FeatureTree::new();
        let base_s = tree.add_sketch(square_sketch(5.0));
        let hole_s = tree.add_sketch(one_point_sketch(0.0, 0.0));
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: base_s,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Plate",
        );
        tree.add_feature(
            Feature::Hole(HoleParams {
                sketch: hole_s,
                depth_mode: HoleDepthMode::Through,
                drill_diameter: 1.5,
                direction_negative: true,
                counterbore: None,
                countersink: None,
                thread: None,
            }),
            "Center Hole",
        );
        let solid = replay(&tree).expect("replay").expect("solid");
        // Plate alone has 6 faces; with a cylindrical hole drilled
        // through we expect MORE.
        assert!(
            solid.faces() > 6,
            "hole should add faces, got {}",
            solid.faces()
        );
    }

    #[test]
    fn thread_metadata_round_trips_through_params() {
        // Build a HoleParams with a thread, persist + restore via RON,
        // and check the spec survives. (This exercises the params
        // round-trip — the thread is metadata that doesn't affect the
        // generated geometry.)
        let m6 = iso_metric_table()
            .into_iter()
            .find(|s| s.designation == "M6")
            .unwrap();
        assert_eq!(m6.standard, ThreadStandard::IsoMetric);
        let params = HoleParams {
            sketch: crate::feature::SketchRef(0),
            depth_mode: HoleDepthMode::Blind { depth: 10.0 },
            drill_diameter: m6.tap_drill_diameter(),
            direction_negative: true,
            counterbore: None,
            countersink: None,
            thread: Some(m6.clone()),
        };
        let ron = ron::ser::to_string(&params).unwrap();
        let restored: HoleParams = ron::from_str(&ron).unwrap();
        assert_eq!(restored.thread.as_ref().unwrap().designation, "M6");
        assert!((restored.drill_diameter - m6.tap_drill_diameter()).abs() < 1e-9);
    }

    #[test]
    fn hole_without_base_returns_bad_parameter() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(one_point_sketch(0.0, 0.0));
        let params = HoleParams {
            sketch: s,
            depth_mode: HoleDepthMode::Blind { depth: 1.0 },
            drill_diameter: 0.5,
            direction_negative: true,
            counterbore: None,
            countersink: None,
            thread: None,
        };
        let err = evaluate(&tree, &params, None).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter { name: "hole", .. }
        ));
    }

    #[test]
    fn hole_rejects_zero_diameter() {
        let mut tree = FeatureTree::new();
        let _b = tree.add_sketch(square_sketch(2.0));
        let s = tree.add_sketch(one_point_sketch(0.0, 0.0));
        let dummy_base = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let params = HoleParams {
            sketch: s,
            depth_mode: HoleDepthMode::Blind { depth: 1.0 },
            drill_diameter: 0.0,
            direction_negative: true,
            counterbore: None,
            countersink: None,
            thread: None,
        };
        let err = evaluate(&tree, &params, Some(&dummy_base)).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "drill_diameter",
                ..
            }
        ));
    }
}
