//! Parametric feature tree — Fusion-style feature history.
//!
//! A [`FeatureTimeline`] is an ordered list of [`Step`]s. Each step takes a
//! modelling [`Feature`] (box / cylinder / extrude), **places** it at a
//! parametric position, and **combines** it with the running body via an
//! [`Op`] (start a new body, join, cut, or intersect). Every dimension and
//! every placement coordinate is a *parameter expression*, resolved against a
//! [`ParameterTable`].
//!
//! [`FeatureTimeline::rebuild`] folds the steps into a single solid — the same
//! base-feature-plus-modifications model Fusion and SolidWorks expose. Edit a
//! named parameter and rebuild, and the whole tree re-drives: a hole moves, a
//! boss grows, the model updates. `rebuild` also returns a snapshot of the
//! running body after each step, so a timeline UI can scrub history.
//!
//! The geometry kernel is `valenx-cad`: primitives (`box_solid` / `cylinder` /
//! `prism`) plus the robust boolean wrappers (`union` / `difference` /
//! `intersection`), which never panic and never return a silently-invalid
//! solid — a failed boolean surfaces here as [`TimelineError::Cad`].

use serde::{Deserialize, Serialize};
use valenx_cad::{
    box_solid, cone, cylinder, difference, intersection, prism, sphere, torus, union, CadError,
    Solid,
};

use crate::parameters::{ParamError, ParameterTable};

/// One modelling primitive. Dimensions are parameter-expression strings
/// (e.g. `"width"`, `"base * 2"`, `"40"`) resolved against a [`ParameterTable`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Feature {
    /// An axis-aligned box `dx × dy × dz`, corner at the step's placement.
    Box {
        /// Width expression (x).
        dx: String,
        /// Depth expression (y).
        dy: String,
        /// Height expression (z).
        dz: String,
    },
    /// A cylinder of the given radius and height, base centred on the step's
    /// placement, axis along +Z.
    Cylinder {
        /// Radius expression.
        radius: String,
        /// Height expression.
        height: String,
    },
    /// A prism: a fixed 2-D profile extruded along +Z to the given height.
    Extrude {
        /// Closed 2-D profile, as `(x, y)` points.
        profile: Vec<(f64, f64)>,
        /// Extrusion-height expression.
        height: String,
    },
    /// A sphere of the given radius, centred on the step's placement.
    Sphere {
        /// Radius expression.
        radius: String,
    },
    /// A (truncated) cone — base radius, top radius (`0` = pointed apex), and
    /// height. Base disk on the placement, axis along +Z.
    Cone {
        /// Base-radius expression.
        base_radius: String,
        /// Top-radius expression (`0` for a pointed cone).
        top_radius: String,
        /// Height expression.
        height: String,
    },
    /// A torus — `major_radius` (centre-circle) and `minor_radius` (tube).
    /// Major axis along +Z, centred on the placement.
    Torus {
        /// Major-radius (centre-circle) expression.
        major_radius: String,
        /// Minor-radius (tube) expression.
        minor_radius: String,
    },
}

/// How a step's feature combines with the body built so far.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Op {
    /// Start a fresh body from this feature, discarding any running body.
    /// The first step of a timeline must be `New`.
    New,
    /// Weld this feature onto the running body (boolean union).
    Join,
    /// Subtract this feature from the running body (boolean difference).
    Cut,
    /// Keep only the overlap of this feature and the running body
    /// (boolean intersection).
    Intersect,
}

/// The no-rotation default (`0°` about each axis), used by [`Step::placed`]
/// and as the serde default so documents written before rotation existed
/// still deserialise.
fn zero_rotation() -> [String; 3] {
    ["0".to_string(), "0".to_string(), "0".to_string()]
}

/// One timeline step: a placed feature plus how it combines with the
/// running body.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Step {
    /// How this feature combines with the body built by earlier steps.
    pub op: Op,
    /// The feature to build at this step.
    pub feature: Feature,
    /// Placement `(x, y, z)` as parameter expressions — where the feature's
    /// origin lands before the combine. `["0", "0", "0"]` for no offset.
    pub at: [String; 3],
    /// Rotation about the feature's own origin, as parameter expressions in
    /// **degrees** about the X, Y, Z axes (applied in that order, before the
    /// placement offset). `["0", "0", "0"]` for no rotation. Defaulted on load
    /// so documents written before rotation existed still parse.
    #[serde(default = "zero_rotation")]
    pub rotate_deg: [String; 3],
}

impl Step {
    /// A step that places `feature` at `(x, y, z)` (parameter expressions)
    /// and combines it with the running body via `op`.
    pub fn placed(
        op: Op,
        feature: Feature,
        x: impl Into<String>,
        y: impl Into<String>,
        z: impl Into<String>,
    ) -> Self {
        Self {
            op,
            feature,
            at: [x.into(), y.into(), z.into()],
            rotate_deg: zero_rotation(),
        }
    }

    /// A step at the origin — placement `(0, 0, 0)`.
    pub fn at_origin(op: Op, feature: Feature) -> Self {
        Self::placed(op, feature, "0", "0", "0")
    }
}

/// The result of rebuilding a timeline: every body in the model, plus a
/// snapshot of the full body set after each step.
///
/// A timeline can hold more than one body: each [`Op::New`] finalises the
/// current body and starts a fresh one (instead of replacing it), so a tree
/// like `New … New …` ends with two separate bodies. `snapshots[k]` is the set
/// of all bodies present after step `k` (`snapshots.len() == steps.len()`, and
/// the last entry equals `bodies`).
pub struct RebuiltModel {
    /// Every body in the final model (one entry per `New` group).
    pub bodies: Vec<Solid>,
    /// The full set of bodies after each step, in order.
    pub snapshots: Vec<Vec<Solid>>,
}

/// An ordered list of parametric steps forming a feature tree.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct FeatureTimeline {
    /// The steps, in build order.
    pub steps: Vec<Step>,
}

/// A failure while rebuilding the timeline.
#[derive(Debug)]
pub enum TimelineError {
    /// A dimension or placement expression failed to resolve.
    Param(ParamError),
    /// A geometry-kernel operation failed (degenerate dimension, an empty or
    /// disjoint boolean result, a mesh-backed boolean operand, …).
    Cad(CadError),
    /// A `Join` / `Cut` / `Intersect` step ran with no body to act on — the
    /// first step of a timeline must be [`Op::New`].
    NoBaseBody,
    /// `rebuild` was called on a timeline with no steps.
    EmptyTimeline,
}

impl std::fmt::Display for TimelineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TimelineError::Param(e) => write!(f, "parameter: {e}"),
            TimelineError::Cad(e) => write!(f, "geometry: {e}"),
            TimelineError::NoBaseBody => {
                write!(f, "the first step must start a new body (Op::New)")
            }
            TimelineError::EmptyTimeline => write!(f, "the timeline has no steps"),
        }
    }
}

impl std::error::Error for TimelineError {}

impl FeatureTimeline {
    /// An empty timeline.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a step.
    pub fn push(&mut self, step: Step) {
        self.steps.push(step);
    }

    /// Number of steps.
    pub fn len(&self) -> usize {
        self.steps.len()
    }

    /// Whether the timeline has no steps.
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Rebuild the whole tree against `params`, folding every step into a
    /// single body. Editing a parameter and rebuilding re-drives the model.
    pub fn rebuild(&self, params: &ParameterTable) -> Result<RebuiltModel, TimelineError> {
        let resolve = |expr: &str| params.compute(expr).map_err(TimelineError::Param);
        // Finalised bodies (from earlier `New`s) plus the body currently being
        // built. `New` finalises `active` and starts a fresh body.
        let mut finished: Vec<Solid> = Vec::new();
        let mut active: Option<Solid> = None;
        let mut snapshots: Vec<Vec<Solid>> = Vec::with_capacity(self.steps.len());

        for step in &self.steps {
            // 1. Build the feature's primitive solid.
            let solid = match &step.feature {
                Feature::Box { dx, dy, dz } => box_solid(resolve(dx)?, resolve(dy)?, resolve(dz)?)
                    .map_err(TimelineError::Cad)?,
                Feature::Cylinder { radius, height } => {
                    cylinder(resolve(radius)?, resolve(height)?).map_err(TimelineError::Cad)?
                }
                Feature::Extrude { profile, height } => {
                    prism(profile, resolve(height)?).map_err(TimelineError::Cad)?
                }
                Feature::Sphere { radius } => {
                    sphere(resolve(radius)?).map_err(TimelineError::Cad)?
                }
                Feature::Cone {
                    base_radius,
                    top_radius,
                    height,
                } => cone(
                    resolve(base_radius)?,
                    resolve(top_radius)?,
                    resolve(height)?,
                )
                .map_err(TimelineError::Cad)?,
                Feature::Torus {
                    major_radius,
                    minor_radius,
                } => torus(resolve(major_radius)?, resolve(minor_radius)?)
                    .map_err(TimelineError::Cad)?,
            };

            // 2a. Orient it: rotate about the feature's own origin (degrees →
            // radians), applying X, then Y, then Z. Identity rotations are
            // skipped so an un-rotated feature isn't needlessly re-transformed.
            let rot = [
                ((1.0, 0.0, 0.0), resolve(&step.rotate_deg[0])?.to_radians()),
                ((0.0, 1.0, 0.0), resolve(&step.rotate_deg[1])?.to_radians()),
                ((0.0, 0.0, 1.0), resolve(&step.rotate_deg[2])?.to_radians()),
            ];
            let mut oriented = solid;
            for (axis, ang) in rot {
                if ang != 0.0 {
                    oriented = oriented
                        .rotated((0.0, 0.0, 0.0), axis, ang)
                        .map_err(TimelineError::Cad)?;
                }
            }
            // 2b. Place it (parametric offset).
            let x = resolve(&step.at[0])?;
            let y = resolve(&step.at[1])?;
            let z = resolve(&step.at[2])?;
            let placed = oriented.translated(x, y, z).map_err(TimelineError::Cad)?;

            // 3. Combine with the active body. `New` finalises the active body
            // (keeping it as a separate body) and starts a fresh one.
            match (active.take(), step.op) {
                (None, Op::New) => active = Some(placed),
                (None, _) => return Err(TimelineError::NoBaseBody),
                (Some(prev), Op::New) => {
                    finished.push(prev);
                    active = Some(placed);
                }
                (Some(r), Op::Join) => {
                    active = Some(union(&r, &placed).map_err(TimelineError::Cad)?)
                }
                (Some(r), Op::Cut) => {
                    active = Some(difference(&r, &placed).map_err(TimelineError::Cad)?)
                }
                (Some(r), Op::Intersect) => {
                    active = Some(intersection(&r, &placed).map_err(TimelineError::Cad)?)
                }
            }

            // Snapshot the full set of bodies present after this step.
            let mut snap = finished.clone();
            if let Some(a) = &active {
                snap.push(a.clone());
            }
            snapshots.push(snap);
        }

        if let Some(a) = active {
            finished.push(a);
        }
        if finished.is_empty() {
            return Err(TimelineError::EmptyTimeline);
        }
        Ok(RebuiltModel {
            bodies: finished,
            snapshots,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> ParameterTable {
        let mut p = ParameterTable::new();
        p.set("size", "1");
        p.set("hole_r", "0.25");
        p.set("hole_h", "2");
        p
    }

    fn unit_box() -> Feature {
        Feature::Box {
            dx: "size".into(),
            dy: "size".into(),
            dz: "size".into(),
        }
    }

    /// The cut/join tests reuse `valenx-cad`'s proven "punched cube" boolean
    /// geometry: a 1×1×1 box at the origin and a cylinder (r=0.25, h=2) whose
    /// base sits at (0.5, 0.5, -0.5), so it straddles the box through its top
    /// and bottom faces — a known-good shapeops input.
    fn through_pin() -> Feature {
        Feature::Cylinder {
            radius: "hole_r".into(),
            height: "hole_h".into(),
        }
    }

    #[test]
    fn base_box_seeds_the_model() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(Op::New, unit_box()));
        let m = tl.rebuild(&params()).expect("rebuild");
        assert_eq!(m.bodies[0].faces(), 6, "a box has six faces");
        assert_eq!(m.snapshots.len(), 1, "one snapshot per step");
    }

    #[test]
    fn cut_punches_a_hole() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(Op::New, unit_box()));
        tl.push(Step::placed(Op::Cut, through_pin(), "0.5", "0.5", "-0.5"));
        let m = tl.rebuild(&params()).expect("punch");
        // A punched cube has the cube's 6 outer faces plus the inner hole
        // walls — strictly more than 6.
        assert!(
            m.bodies[0].faces() > 6,
            "punched body should have hole walls ({} vs 6)",
            m.bodies[0].faces()
        );
        assert_eq!(m.snapshots.len(), 2);
    }

    #[test]
    fn join_welds_a_protrusion() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(Op::New, unit_box()));
        tl.push(Step::placed(Op::Join, through_pin(), "0.5", "0.5", "-0.5"));
        let m = tl.rebuild(&params()).expect("weld");
        assert!(m.bodies[0].faces() > 0, "welded body is a real solid");
    }

    #[test]
    fn cut_without_a_base_body_errors() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(Op::Cut, unit_box()));
        assert!(matches!(
            tl.rebuild(&params()),
            Err(TimelineError::NoBaseBody)
        ));
    }

    #[test]
    fn empty_timeline_errors() {
        let tl = FeatureTimeline::new();
        assert!(matches!(
            tl.rebuild(&params()),
            Err(TimelineError::EmptyTimeline)
        ));
    }

    #[test]
    fn undefined_parameter_errors() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(
            Op::New,
            Feature::Cylinder {
                radius: "missing".into(),
                height: "1".into(),
            },
        ));
        assert!(matches!(
            tl.rebuild(&ParameterTable::new()),
            Err(TimelineError::Param(_))
        ));
    }

    #[test]
    fn editing_a_parameter_redrives_the_build() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(Op::New, unit_box()));
        let mut p = ParameterTable::new();
        p.set("size", "1");
        assert!(tl.rebuild(&p).is_ok(), "a valid size builds");
        // Edit the parameter to a broken expression — re-resolved on rebuild,
        // so the edit surfaces as an error (proving the value is re-read).
        p.set("size", "1 +");
        assert!(
            matches!(tl.rebuild(&p), Err(TimelineError::Param(_))),
            "the edited parameter is re-resolved each rebuild"
        );
    }

    #[test]
    fn sphere_cone_torus_build_as_standalone_bodies() {
        let mut p = ParameterTable::new();
        p.set("r", "0.5");
        let cases = [
            Feature::Sphere { radius: "r".into() },
            Feature::Cone {
                base_radius: "r".into(),
                top_radius: "0".into(),
                height: "1".into(),
            },
            Feature::Torus {
                major_radius: "1".into(),
                minor_radius: "r".into(),
            },
        ];
        for feat in cases {
            let mut tl = FeatureTimeline::new();
            tl.push(Step::at_origin(Op::New, feat));
            let m = tl.rebuild(&p).expect("primitive builds");
            assert!(m.bodies[0].faces() > 0, "primitive should have faces");
        }
    }

    #[test]
    fn rotation_is_applied_and_validated() {
        let mut p = ParameterTable::new();
        p.set("ang", "45");
        // A rotated box is still a box (6 faces) — the rotation path runs cleanly.
        let mut step = Step::at_origin(
            Op::New,
            Feature::Box {
                dx: "1".into(),
                dy: "1".into(),
                dz: "1".into(),
            },
        );
        step.rotate_deg = ["ang".into(), "0".into(), "0".into()];
        let mut tl = FeatureTimeline::new();
        tl.push(step);
        let m = tl.rebuild(&p).expect("rotated box builds");
        assert_eq!(m.bodies[0].faces(), 6, "a rotated box is still a box");

        // A broken rotation expression surfaces as a Param error, not a panic.
        let mut bad = Step::at_origin(
            Op::New,
            Feature::Box {
                dx: "1".into(),
                dy: "1".into(),
                dz: "1".into(),
            },
        );
        bad.rotate_deg = ["1 +".into(), "0".into(), "0".into()];
        let mut tl2 = FeatureTimeline::new();
        tl2.push(bad);
        assert!(matches!(
            tl2.rebuild(&ParameterTable::new()),
            Err(TimelineError::Param(_))
        ));
    }

    #[test]
    fn extrude_builds_from_a_profile() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(
            Op::New,
            Feature::Extrude {
                profile: vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)],
                height: "1".into(),
            },
        ));
        let m = tl.rebuild(&ParameterTable::new()).expect("extrude builds");
        assert!(
            m.bodies[0].faces() >= 5,
            "an extruded quad is a box-like solid"
        );
    }

    #[test]
    fn new_keeps_a_separate_body() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(
            Op::New,
            Feature::Box {
                dx: "1".into(),
                dy: "1".into(),
                dz: "1".into(),
            },
        ));
        // A second New starts a separate body instead of replacing the first.
        tl.push(Step::placed(
            Op::New,
            Feature::Box {
                dx: "1".into(),
                dy: "1".into(),
                dz: "1".into(),
            },
            "3",
            "0",
            "0",
        ));
        let m = tl
            .rebuild(&ParameterTable::new())
            .expect("two bodies build");
        assert_eq!(m.bodies.len(), 2, "a second New keeps a separate body");
        assert_eq!(m.snapshots[0].len(), 1, "one body after the first step");
        assert_eq!(
            m.snapshots.last().unwrap().len(),
            2,
            "both bodies in the final snapshot"
        );
    }
}
