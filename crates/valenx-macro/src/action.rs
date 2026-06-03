//! Macro action enum — every user-visible Valenx UI action this v1
//! recorder knows how to capture and replay.
//!
//! Adding a new action: pick a variant name, add a struct field for
//! its parameters, then handle it in three places — `dispatcher`
//! (replay-time effect), `export::to_python` (Python emission), and
//! the UI panel that should call `record_action(...)` when the user
//! clicks the corresponding button.

use serde::{Deserialize, Serialize};
use valenx_feature_tree::feature::{Feature, FeatureId, SketchRef};

/// Logical panel identifier — matches the panel IDs the desktop shell
/// uses when calling [`crate::recorder::MacroRecorder::record`] for a
/// panel switch.
#[derive(Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum PanelId {
    /// 2D parametric sketcher.
    Sketcher,
    /// Part Design feature tree.
    PartDesign,
    /// Mesh workbench (booleans, smoothing, etc.).
    Mesh,
    /// TechDraw 2D-from-3D.
    TechDraw,
    /// Assembly + joint solver.
    Assembly,
    /// Surface (NURBS) workbench.
    Surface,
    /// Path / CAM workbench.
    Cam,
    /// Arch / BIM workbench.
    Arch,
    /// Spreadsheet workbench.
    Spreadsheet,
    /// Macro Library panel (Phase 21).
    MacroLibrary,
    /// Add-on Manager (Phase 22).
    AddonManager,
    /// FEM workbench (Phase 24).
    Fem,
    /// Any other named panel (room for community workbenches).
    Other(String),
}

impl PanelId {
    /// Stable string key — used by both the Python emitter and the RON
    /// serializer.
    pub fn as_key(&self) -> &str {
        match self {
            PanelId::Sketcher => "sketcher",
            PanelId::PartDesign => "part_design",
            PanelId::Mesh => "mesh",
            PanelId::TechDraw => "tech_draw",
            PanelId::Assembly => "assembly",
            PanelId::Surface => "surface",
            PanelId::Cam => "cam",
            PanelId::Arch => "arch",
            PanelId::Spreadsheet => "spreadsheet",
            PanelId::MacroLibrary => "macro_library",
            PanelId::AddonManager => "addon_manager",
            PanelId::Fem => "fem",
            PanelId::Other(s) => s.as_str(),
        }
    }
}

/// Kind of sketch entity created — keeps the action enum small while
/// still carrying enough info to replay deterministically.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum EntityKind {
    /// Point at `(x, y)`.
    Point {
        /// X coordinate in the sketch plane.
        x: f64,
        /// Y coordinate.
        y: f64,
    },
    /// Line from `(x1, y1)` to `(x2, y2)`.
    Line {
        /// Start X.
        x1: f64,
        /// Start Y.
        y1: f64,
        /// End X.
        x2: f64,
        /// End Y.
        y2: f64,
    },
    /// Circle at `(cx, cy)` with `radius`.
    Circle {
        /// Center X.
        cx: f64,
        /// Center Y.
        cy: f64,
        /// Radius.
        radius: f64,
    },
    /// Arc from `(x1, y1)` to `(x2, y2)` curving through `(xmid, ymid)`.
    Arc {
        /// Start X.
        x1: f64,
        /// Start Y.
        y1: f64,
        /// Mid X.
        xmid: f64,
        /// Mid Y.
        ymid: f64,
        /// End X.
        x2: f64,
        /// End Y.
        y2: f64,
    },
}

/// One user action captured by the recorder.
///
/// New variants must also be handled in `dispatcher::dispatch` and
/// `export::to_python`. The enum is `#[non_exhaustive]` so adding a
/// variant in a future phase isn't a breaking change for downstream
/// crates that match on it.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[non_exhaustive]
pub enum MacroAction {
    /// Switch the active panel.
    SwitchPanel {
        /// Target panel.
        panel_id: PanelId,
    },
    /// Add a sketch entity to the active sketcher.
    AddSketchEntity {
        /// Which entity to add.
        entity: EntityKind,
    },
    /// Solve the active sketch (run the constraint solver).
    SolveSketch {
        /// Reference to the sketch in the feature tree.
        sketch: SketchRef,
    },
    /// Add a feature to the feature tree.
    AddFeature {
        /// The feature to add (held by value).
        feature: Feature,
        /// User-visible label.
        name: String,
    },
    /// Extrude (pad) the named sketch by `depth`.
    ///
    /// Shortcut for `AddFeature(Feature::Pad(...))` for the most-common
    /// authoring step — easier to round-trip through Python.
    ExtrudeSketch {
        /// Source sketch.
        sketch: SketchRef,
        /// Extrusion depth.
        depth: f64,
    },
    /// Save the current project to `path`.
    SaveProject {
        /// Filesystem path (`.valenx`).
        path: String,
    },
    /// Load a project from `path`, replacing the current state.
    LoadProject {
        /// Filesystem path.
        path: String,
    },
    /// Run a Python script in the embedded interpreter (valenx-py).
    ///
    /// Captured for replay so a macro can chain its own Python step
    /// (useful for parametric loops the macro UI doesn't model
    /// natively).
    RunPython {
        /// Verbatim Python source.
        script: String,
    },
    /// Suppress / unsuppress a feature in the tree.
    SetFeatureSuppressed {
        /// Target feature.
        id: FeatureId,
        /// Suppressed state.
        suppressed: bool,
    },
    /// Tag for grouping a sequence of low-level events into a single
    /// "user-perceived undo step". The recorder emits this when the
    /// user clicks a button that fans out into several smaller actions.
    Step {
        /// Free-form label.
        label: String,
    },
}

impl MacroAction {
    /// Short user-facing label for the Macro Library tree-view kind
    /// column.
    pub fn kind_label(&self) -> &'static str {
        match self {
            MacroAction::SwitchPanel { .. } => "switch_panel",
            MacroAction::AddSketchEntity { .. } => "add_sketch_entity",
            MacroAction::SolveSketch { .. } => "solve_sketch",
            MacroAction::AddFeature { .. } => "add_feature",
            MacroAction::ExtrudeSketch { .. } => "extrude_sketch",
            MacroAction::SaveProject { .. } => "save_project",
            MacroAction::LoadProject { .. } => "load_project",
            MacroAction::RunPython { .. } => "run_python",
            MacroAction::SetFeatureSuppressed { .. } => "set_feature_suppressed",
            MacroAction::Step { .. } => "step",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_id_keys_are_stable() {
        assert_eq!(PanelId::Sketcher.as_key(), "sketcher");
        assert_eq!(PanelId::PartDesign.as_key(), "part_design");
        assert_eq!(PanelId::Other("foo".into()).as_key(), "foo");
    }

    #[test]
    fn action_kind_label_is_distinct_per_variant() {
        let actions: Vec<MacroAction> = vec![
            MacroAction::SwitchPanel {
                panel_id: PanelId::Mesh,
            },
            MacroAction::AddSketchEntity {
                entity: EntityKind::Point { x: 0.0, y: 0.0 },
            },
            MacroAction::SolveSketch {
                sketch: SketchRef(0),
            },
            MacroAction::SaveProject {
                path: "/tmp/p".into(),
            },
            MacroAction::RunPython {
                script: "print('hi')".into(),
            },
        ];
        let labels: Vec<&'static str> = actions.iter().map(|a| a.kind_label()).collect();
        // Distinct labels — no duplicates.
        let mut set = std::collections::HashSet::new();
        for l in &labels {
            assert!(set.insert(*l), "duplicate label {l}");
        }
    }

    #[test]
    fn ron_round_trip_action() {
        let a = MacroAction::AddSketchEntity {
            entity: EntityKind::Circle {
                cx: 1.0,
                cy: 2.0,
                radius: 3.0,
            },
        };
        let txt = ron::ser::to_string(&a).unwrap();
        let b: MacroAction = ron::from_str(&txt).unwrap();
        match b {
            MacroAction::AddSketchEntity { entity } => match entity {
                EntityKind::Circle { cx, cy, radius } => {
                    assert!((cx - 1.0).abs() < 1e-9);
                    assert!((cy - 2.0).abs() < 1e-9);
                    assert!((radius - 3.0).abs() < 1e-9);
                }
                _ => panic!("entity kind wrong"),
            },
            _ => panic!("action kind wrong"),
        }
    }
}
