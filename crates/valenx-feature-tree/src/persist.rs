//! `.valenx` project persistence — a RON envelope wrapping the
//! [`FeatureTree`] with format-version metadata.
//!
//! Modeled after `valenx_sketch::persist::SketchFile`. Bumping
//! [`ValenxProject::VERSION`] flags an on-disk schema change.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::tree::FeatureTree;
use crate::FeatureError;

/// On-disk envelope wrapping a feature tree with format-version
/// metadata. The tree itself owns its sketches, so a round-trip
/// through [`Self::write_to`] / [`Self::read_from`] is self-contained.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ValenxProject {
    /// Format version — bumped when the on-disk schema changes.
    pub version: u32,
    /// The feature tree payload (sketches + features).
    pub feature_tree: FeatureTree,
}

impl ValenxProject {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a tree as a project envelope at the current version.
    pub fn from_tree(tree: &FeatureTree) -> Self {
        Self {
            version: Self::VERSION,
            feature_tree: tree.clone(),
        }
    }

    /// Serialize to a pretty-printed RON string.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::Ron`] if RON serialization fails. In
    /// practice this means a non-`Serialize`-clean value snuck into a
    /// feature parameter — investigate the offending variant.
    pub fn to_ron(&self) -> Result<String, FeatureError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| FeatureError::Ron(e.to_string()))
    }

    /// Parse a project envelope from a RON string.
    ///
    /// R33 H1: after structural deserialization every sketch in the
    /// loaded tree is run through [`FeatureTree::validate`], so a
    /// hand-edited / corrupt / version-skewed `.valenx` carrying an
    /// out-of-range sketch variable handle is rejected here rather than
    /// panicking "index out of bounds" during replay on the UI thread
    /// (which has no `catch_unwind`). Covers the
    /// `ValenxProject::read_from` load path used by the Part Design
    /// panel's "Open Valenx project".
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::Ron`] if the input is not valid RON or
    /// does not match the expected shape, or
    /// [`FeatureError::SketchError`] if a loaded sketch fails handle
    /// validation.
    pub fn from_ron(s: &str) -> Result<Self, FeatureError> {
        let project: Self = ron::from_str(s).map_err(|e| FeatureError::Ron(e.to_string()))?;
        project.feature_tree.validate()?;
        Ok(project)
    }

    /// Write the envelope to `path` as RON.
    ///
    /// Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::Ron`] for serialization failures or
    /// [`FeatureError::Io`] for filesystem failures.
    pub fn write_to(&self, path: &Path) -> Result<(), FeatureError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Read an envelope from `path`.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::Io`] when the file cannot be read or
    /// [`FeatureError::Ron`] when parsing fails.
    pub fn read_from(path: &Path) -> Result<Self, FeatureError> {
        // R29 D: canonical valenx_core::io_caps::read_capped_to_string at
        // MAX_DOC_FILE_BYTES (16 MiB), replacing the private dupe.
        let s = valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?;
        Self::from_ron(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, FeatureId, MirrorParams, PadParams, PocketParams};
    use nalgebra::Vector3;

    #[test]
    fn round_trips_empty_tree() {
        let tree = FeatureTree::new();
        let project = ValenxProject::from_tree(&tree);
        assert_eq!(project.version, ValenxProject::VERSION);
        let ron = project.to_ron().expect("serialize");
        assert!(ron.contains("version: 1"));
        let parsed = ValenxProject::from_ron(&ron).expect("parse");
        assert_eq!(parsed.version, ValenxProject::VERSION);
        assert_eq!(parsed.feature_tree.sketches.len(), 0);
        assert_eq!(parsed.feature_tree.features.len(), 0);
    }

    #[test]
    fn round_trips_tree_with_one_pad() {
        let mut tree = FeatureTree::new();
        // Sketch with a single point so the persisted payload is not
        // trivially zero-sized.
        let mut sk = valenx_sketch::Sketch::new();
        sk.add_point(1.0, 2.0);
        let s_ref = tree.add_sketch(sk);
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s_ref,
                depth: 7.5.into(),
                direction_positive: true,
            }),
            "Bottom Pad",
        );

        let ron = ValenxProject::from_tree(&tree).to_ron().expect("serialize");
        let parsed = ValenxProject::from_ron(&ron).expect("parse");

        assert_eq!(parsed.version, ValenxProject::VERSION);
        assert_eq!(parsed.feature_tree.sketches.len(), 1);
        assert_eq!(parsed.feature_tree.features.len(), 1);
        assert_eq!(parsed.feature_tree.features[0].name, "Bottom Pad");
        assert!(!parsed.feature_tree.features[0].suppressed);
        match &parsed.feature_tree.features[0].feature {
            Feature::Pad(p) => {
                assert_eq!(p.depth.literal(), Some(7.5));
                assert!(p.direction_positive);
                assert_eq!(p.sketch.0, 0);
            }
            other => panic!("expected Pad, got {other:?}"),
        }
        assert_eq!(parsed.feature_tree.sketches[0].vars, vec![1.0, 2.0]);
    }

    #[test]
    fn round_trips_phase_13_variants() {
        // Phase 13G Task 47 — every new Feature variant survives a RON
        // round-trip with field values intact. Build a tree that uses
        // each of the 10 new variants at least once, persist to a
        // string, parse it back, and spot-check a couple of distinctive
        // params per variant.
        use crate::feature::{
            BoolKind, BooleanHistoryParams, CounterboreParams, CountersinkParams, DraftAngleParams,
            Feature, HelixParams, HoleDepthMode, HoleParams, LoftParams, MultiTransformParams,
            PipeParams, ShellParams, ShellSide, SweepParams, ThicknessParams, TransformOp,
        };
        use crate::threads::{iso_metric_table, ThreadStandard};

        let mut tree = FeatureTree::new();
        // Two sketches: one with points (for Hole), one square profile.
        let mut sk_points = valenx_sketch::Sketch::new();
        sk_points.add_point(1.0, 2.0);
        sk_points.add_point(-1.0, -2.0);
        let s_pts = tree.add_sketch(sk_points);
        let mut sk_square = valenx_sketch::Sketch::new();
        let a = sk_square.add_point(0.0, 0.0);
        let b = sk_square.add_point(1.0, 0.0);
        let c = sk_square.add_point(1.0, 1.0);
        let d = sk_square.add_point(0.0, 1.0);
        sk_square.add_line(a, b).unwrap();
        sk_square.add_line(b, c).unwrap();
        sk_square.add_line(c, d).unwrap();
        sk_square.add_line(d, a).unwrap();
        let s_sq = tree.add_sketch(sk_square);

        // Pad to anchor later cross-feature refs.
        let pad_id = tree.add_feature(
            Feature::Pad(crate::feature::PadParams {
                sketch: s_sq,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Anchor",
        );

        let m6 = iso_metric_table()
            .into_iter()
            .find(|s| s.designation == "M6")
            .unwrap();
        tree.add_feature(
            Feature::Hole(HoleParams {
                sketch: s_pts,
                depth_mode: HoleDepthMode::Blind { depth: 10.0 },
                drill_diameter: m6.tap_drill_diameter(),
                direction_negative: true,
                counterbore: Some(CounterboreParams {
                    diameter: 8.0,
                    depth: 2.0,
                }),
                countersink: Some(CountersinkParams {
                    diameter: 8.0,
                    angle_deg: 82.0,
                }),
                thread: Some(m6.clone()),
            }),
            "Hole",
        );
        tree.add_feature(
            Feature::Loft(LoftParams {
                profile_sketches: vec![s_sq, s_pts],
                guide_curves: vec![],
                closed: true,
                ruled: false,
            }),
            "Loft",
        );
        tree.add_feature(
            Feature::Sweep(SweepParams {
                profile_sketch: s_sq,
                path_sketch: s_pts,
                twist_angle: std::f64::consts::FRAC_PI_2,
                keep_profile_orientation: false,
            }),
            "Sweep",
        );
        tree.add_feature(
            Feature::Pipe(PipeParams {
                cross_section_sketch: s_sq,
                centerline_sketch: s_pts,
                bend_radius: 0.5,
            }),
            "Pipe",
        );
        tree.add_feature(
            Feature::Helix(HelixParams {
                profile_sketch: s_sq,
                pitch: 1.5,
                turns: 4.0,
                axis_origin: Vector3::new(0.0, 0.0, 0.0),
                axis_direction: Vector3::new(0.0, 0.0, 1.0),
                taper_angle: 0.0,
                left_handed: true,
            }),
            "Helix",
        );
        tree.add_feature(
            Feature::MultiTransform(MultiTransformParams {
                target: pad_id,
                transforms: vec![
                    TransformOp::Translate {
                        delta: Vector3::new(5.0, 0.0, 0.0),
                    },
                    TransformOp::Rotate {
                        axis: Vector3::new(0.0, 0.0, 1.0),
                        angle_rad: 1.0,
                    },
                    TransformOp::Scale { factor: 2.0 },
                    TransformOp::Mirror {
                        plane_normal: Vector3::new(1.0, 0.0, 0.0),
                    },
                ],
            }),
            "Multi",
        );
        tree.add_feature(
            Feature::DraftAngle(DraftAngleParams {
                target: pad_id,
                face_indices: vec![0, 1, 2],
                neutral_plane_normal: Vector3::new(0.0, 0.0, 1.0),
                draft_angle_deg: 5.0,
            }),
            "Draft",
        );
        tree.add_feature(
            Feature::Shell(ShellParams {
                target: pad_id,
                face_indices_to_remove: vec![],
                thickness: 0.1,
                inward_or_outward: ShellSide::Inward,
            }),
            "Shell",
        );
        tree.add_feature(
            Feature::Thickness(ThicknessParams {
                target: pad_id,
                face_index: 0,
                thickness: 0.5,
            }),
            "Thicken",
        );
        tree.add_feature(
            Feature::BooleanHistory(BooleanHistoryParams {
                operation: BoolKind::Union,
                targets: vec![pad_id, FeatureId(1)],
            }),
            "BoolHist",
        );

        let ron = ValenxProject::from_tree(&tree).to_ron().expect("serialize");
        let parsed = ValenxProject::from_ron(&ron).expect("parse");

        // 1 anchor + 10 new variants = 11 features.
        assert_eq!(parsed.feature_tree.features.len(), 11);

        // Spot-check distinctive params.
        match &parsed.feature_tree.features[1].feature {
            Feature::Hole(h) => {
                assert!(h.thread.is_some());
                assert_eq!(
                    h.thread.as_ref().unwrap().standard,
                    ThreadStandard::IsoMetric
                );
                assert!(h.counterbore.is_some());
            }
            other => panic!("expected Hole, got {other:?}"),
        }
        match &parsed.feature_tree.features[5].feature {
            Feature::Helix(h) => {
                assert!(h.left_handed);
                assert!((h.pitch - 1.5).abs() < 1e-12);
            }
            other => panic!("expected Helix, got {other:?}"),
        }
        match &parsed.feature_tree.features[6].feature {
            Feature::MultiTransform(mt) => {
                assert_eq!(mt.transforms.len(), 4);
            }
            other => panic!("expected Multi, got {other:?}"),
        }
        match &parsed.feature_tree.features[10].feature {
            Feature::BooleanHistory(bh) => {
                assert_eq!(bh.operation, BoolKind::Union);
                assert_eq!(bh.targets.len(), 2);
            }
            other => panic!("expected BooleanHistory, got {other:?}"),
        }
    }

    #[test]
    fn round_trips_phase_14_fillet_with_edge_indices() {
        // Phase 14 Task 25 — FilletParams.edge_indices round-trips
        // through RON in both Some and None forms.
        use crate::feature::{ChamferParams, FilletParams};

        let mut tree = FeatureTree::new();
        let mut sk = valenx_sketch::Sketch::new();
        sk.add_point(0.0, 0.0);
        let s = tree.add_sketch(sk);
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Base",
        );
        // Some(indices): explicit BRep edge selection.
        tree.add_feature(
            Feature::Fillet(FilletParams {
                target: pad_id,
                radius: 0.25,
                threshold_deg: 30.0,
                edge_indices: Some(vec![0, 1, 4, 7]),
            }),
            "Fillet selected",
        );
        // None: auto-by-angle, the Phase 3 default.
        tree.add_feature(
            Feature::Chamfer(ChamferParams {
                target: pad_id,
                distance: 0.1,
                threshold_deg: 45.0,
                edge_indices: None,
            }),
            "Chamfer auto",
        );

        let ron = ValenxProject::from_tree(&tree).to_ron().expect("serialize");
        let parsed = ValenxProject::from_ron(&ron).expect("parse");
        assert_eq!(parsed.feature_tree.features.len(), 3);
        match &parsed.feature_tree.features[1].feature {
            Feature::Fillet(p) => {
                assert!((p.radius - 0.25).abs() < 1e-12);
                assert_eq!(p.edge_indices.as_deref(), Some(&[0, 1, 4, 7][..]));
            }
            other => panic!("expected Fillet, got {other:?}"),
        }
        match &parsed.feature_tree.features[2].feature {
            Feature::Chamfer(p) => {
                assert!((p.distance - 0.1).abs() < 1e-12);
                assert!(
                    p.edge_indices.is_none(),
                    "edge_indices = None should round-trip as None"
                );
            }
            other => panic!("expected Chamfer, got {other:?}"),
        }
    }

    #[test]
    fn from_ron_rejects_sketch_with_out_of_range_var_handle() {
        // R33 H1: a hand-edited project whose sketch carries a Point
        // referencing variable index 999 over a single-var `vars` must
        // be rejected at load (FeatureError::SketchError /
        // sketch.corrupt_handle), NOT deserialize successfully and then
        // panic "index out of bounds" during replay on the UI thread.
        let ron = r#"(
    version: 1,
    feature_tree: (
        sketches: [
            (
                vars: [0.0],
                entities: [
                    Point((x_var: 999, y_var: 0)),
                ],
                constraints: [],
            ),
        ],
        features: [],
    ),
)"#;
        let err = ValenxProject::from_ron(ron)
            .expect_err("corrupt sketch handle must be rejected at project load");
        assert_eq!(err.code(), "feature_tree.sketch_error");
    }

    #[test]
    fn from_ron_accepts_project_with_well_formed_sketch() {
        // The same shape but with an in-range handle parses fine — the
        // validation gate does not reject legitimate documents.
        let ron = r#"(
    version: 1,
    feature_tree: (
        sketches: [
            (
                vars: [1.0, 2.0],
                entities: [
                    Point((x_var: 0, y_var: 1)),
                ],
                constraints: [],
            ),
        ],
        features: [],
    ),
)"#;
        let parsed = ValenxProject::from_ron(ron).expect("well-formed project must load");
        assert_eq!(parsed.feature_tree.sketches.len(), 1);
        assert_eq!(parsed.feature_tree.sketches[0].vars, vec![1.0, 2.0]);
    }

    #[test]
    fn loads_phase_3_ron_without_edge_indices_field() {
        // Phase 14 Task 25 backward compatibility — a RON document
        // written by Phase 3 (no edge_indices field) should parse
        // with edge_indices = None thanks to #[serde(default)].
        let ron = r#"(
    version: 1,
    feature_tree: (
        sketches: [],
        features: [
            (
                name: "F",
                suppressed: false,
                feature: Fillet((
                    target: (0),
                    radius: 0.2,
                    threshold_deg: 45.0,
                )),
            ),
        ],
    ),
)"#;
        let parsed =
            ValenxProject::from_ron(ron).expect("legacy RON without edge_indices should parse");
        assert_eq!(parsed.feature_tree.features.len(), 1);
        match &parsed.feature_tree.features[0].feature {
            Feature::Fillet(p) => {
                assert!(p.edge_indices.is_none());
                assert!((p.radius - 0.2).abs() < 1e-12);
            }
            other => panic!("expected Fillet, got {other:?}"),
        }
    }

    #[test]
    fn round_trips_full_tree_pad_pocket_mirror() {
        // Two sketches: profile for the pad + profile for the pocket.
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(valenx_sketch::Sketch::new());
        let s1 = tree.add_sketch(valenx_sketch::Sketch::new());

        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s0,
                depth: 10.0.into(),
                direction_positive: true,
            }),
            "Base Pad",
        );
        let _pocket_id = tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: s1,
                depth: 4.0.into(),
                direction_positive: false,
            }),
            "Pocket 1",
        );
        let _mirror_id = tree.add_feature(
            Feature::Mirror(MirrorParams {
                target: pad_id,
                plane_origin: Vector3::new(0.0, 0.0, 0.0),
                plane_normal: Vector3::new(1.0, 0.0, 0.0),
                keep_original: true,
            }),
            "Mirror across YZ",
        );

        // Suppress the pocket to verify the flag round-trips too.
        tree.set_suppressed(FeatureId(1), true).unwrap();

        // Persist via filesystem to cover write_to / read_from.
        let tmp = std::env::temp_dir().join("valenx_project_full_tree.valenx");
        ValenxProject::from_tree(&tree).write_to(&tmp).unwrap();
        let parsed = ValenxProject::read_from(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);

        assert_eq!(parsed.version, ValenxProject::VERSION);
        assert_eq!(parsed.feature_tree.sketches.len(), 2);
        assert_eq!(parsed.feature_tree.features.len(), 3);

        assert_eq!(parsed.feature_tree.features[0].name, "Base Pad");
        assert!(!parsed.feature_tree.features[0].suppressed);
        assert!(matches!(
            parsed.feature_tree.features[0].feature,
            Feature::Pad(_)
        ));

        assert_eq!(parsed.feature_tree.features[1].name, "Pocket 1");
        assert!(parsed.feature_tree.features[1].suppressed);
        match &parsed.feature_tree.features[1].feature {
            Feature::Pocket(p) => {
                assert_eq!(p.depth.literal(), Some(4.0));
                assert!(!p.direction_positive);
                assert_eq!(p.sketch.0, 1);
            }
            other => panic!("expected Pocket, got {other:?}"),
        }

        assert_eq!(parsed.feature_tree.features[2].name, "Mirror across YZ");
        match &parsed.feature_tree.features[2].feature {
            Feature::Mirror(m) => {
                assert_eq!(m.target, pad_id);
                assert!(m.keep_original);
                assert_eq!(m.plane_normal, Vector3::new(1.0, 0.0, 0.0));
            }
            other => panic!("expected Mirror, got {other:?}"),
        }
    }
}
