//! Versioned RON persistence for the solver panel state.

use serde::{Deserialize, Serialize};

use crate::error::Solve3DError;
use crate::sketch::Sketch3D;
use crate::timeline::FeatureTimeline;

/// Schema version.
pub const VERSION: u32 = 1;

/// RON envelope storing the sketch (entities + constraints + vars).
/// The solver report is intentionally NOT persisted — re-solving on
/// load gives a deterministic, up-to-date snapshot.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PanelFile {
    /// Schema version.
    pub version: u32,
    /// The persisted sketch.
    pub sketch: Sketch3D,
}

/// Serialise.
pub fn to_ron_string(s: &Sketch3D) -> Result<String, Solve3DError> {
    let f = PanelFile {
        version: VERSION,
        sketch: s.clone(),
    };
    ron::ser::to_string(&f).map_err(|e| Solve3DError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })
}

/// Deserialise.
pub fn from_ron_str(s: &str) -> Result<PanelFile, Solve3DError> {
    let f: PanelFile = ron::from_str(s).map_err(|e| Solve3DError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })?;
    if f.version != VERSION {
        return Err(Solve3DError::BadParameter {
            name: "version",
            reason: format!("unsupported {} expected {VERSION}", f.version),
        });
    }
    Ok(f)
}

/// Schema version for a saved feature-tree document.
pub const TIMELINE_VERSION: u32 = 1;

/// RON envelope storing a parametric feature tree: the named parameters
/// (as `(name, expression)` rows) plus the [`FeatureTimeline`]. Rebuilding
/// the timeline against the parameters on load reproduces the model.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TimelineFile {
    /// Schema version.
    pub version: u32,
    /// Named parameters as `(name, expression)` rows.
    pub parameters: Vec<(String, String)>,
    /// The persisted feature timeline.
    pub timeline: FeatureTimeline,
}

/// Serialise a feature tree (parameters + timeline) to RON.
pub fn timeline_to_ron(
    parameters: &[(String, String)],
    timeline: &FeatureTimeline,
) -> Result<String, Solve3DError> {
    let f = TimelineFile {
        version: TIMELINE_VERSION,
        parameters: parameters.to_vec(),
        timeline: timeline.clone(),
    };
    ron::ser::to_string(&f).map_err(|e| Solve3DError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })
}

/// Deserialise a feature tree from RON. Rejects an unknown schema version.
pub fn timeline_from_ron(s: &str) -> Result<TimelineFile, Solve3DError> {
    let f: TimelineFile = ron::from_str(s).map_err(|e| Solve3DError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })?;
    if f.version != TIMELINE_VERSION {
        return Err(Solve3DError::BadParameter {
            name: "version",
            reason: format!("unsupported {} expected {TIMELINE_VERSION}", f.version),
        });
    }
    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::Constraint3D;
    use crate::timeline::{Feature, Op, Step};
    use crate::ParameterTable;

    #[test]
    fn round_trip_empty_sketch() {
        let s = Sketch3D::new();
        let txt = to_ron_string(&s).unwrap();
        let f = from_ron_str(&txt).unwrap();
        assert_eq!(f.version, VERSION);
        assert_eq!(f.sketch.entities.len(), 0);
    }

    #[test]
    fn round_trip_with_constraint() {
        let mut s = Sketch3D::new();
        let a = s.add_point(0.0, 0.0, 0.0);
        let b = s.add_point(1.0, 2.0, 3.0);
        s.add_constraint(Constraint3D::Coincident3 { a, b });
        let txt = to_ron_string(&s).unwrap();
        let f = from_ron_str(&txt).unwrap();
        assert_eq!(f.sketch.constraints.len(), 1);
        assert_eq!(f.sketch.entities.len(), 2);
        assert_eq!(f.sketch.vars.len(), 6);
    }

    #[test]
    fn round_trip_feature_timeline() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(
            Op::New,
            Feature::Box { dx: "size".into(), dy: "size".into(), dz: "size".into() },
        ));
        tl.push(Step::placed(
            Op::Cut,
            Feature::Cylinder { radius: "hole_r".into(), height: "hole_h".into() },
            "size / 2",
            "size / 2",
            "-0.5",
        ));
        let params = vec![
            ("size".to_string(), "1".to_string()),
            ("hole_r".to_string(), "0.25".to_string()),
            ("hole_h".to_string(), "2".to_string()),
        ];

        let txt = timeline_to_ron(&params, &tl).unwrap();
        let f = timeline_from_ron(&txt).unwrap();
        assert_eq!(f.version, TIMELINE_VERSION);
        assert_eq!(f.parameters.len(), 3);
        assert_eq!(f.timeline.len(), 2);

        // The deserialised tree rebuilds to the same punched cube.
        let mut table = ParameterTable::new();
        for (n, e) in &f.parameters {
            table.set(n, e);
        }
        let model = f.timeline.rebuild(&table).unwrap();
        assert!(model.bodies[0].faces() > 6, "round-tripped tree still punches a hole");
    }

    #[test]
    fn timeline_rejects_unknown_version() {
        let txt = timeline_to_ron(&[], &FeatureTimeline::new()).unwrap();
        assert!(timeline_from_ron(&txt).is_ok());
        let bumped = txt.replacen("version:1", "version:2", 1);
        assert_ne!(txt, bumped, "version token should be present to corrupt");
        assert!(timeline_from_ron(&bumped).is_err());
    }
}
