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
    // serde builds the sketch field-by-field, bypassing the `add_*`
    // builders' reference checks. Validate every entity/constraint
    // reference now so a crafted/corrupt file can't panic the solver's
    // raw-indexing accessors on the re-solve that follows a load.
    f.sketch.validate()?;
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
    fn from_ron_rejects_constraint_with_out_of_range_entity() {
        // serde builds the sketch field-by-field, bypassing the `add_*`
        // reference checks. A constraint referencing entities that don't
        // exist loaded fine pre-fix, then panicked the solver's
        // `entities[id.0]` accessor on re-solve. It must be rejected now.
        let mut s = Sketch3D::new();
        s.constraints.push(Constraint3D::Coincident3 {
            a: crate::entity::EntityId(0),
            b: crate::entity::EntityId(1),
        });
        let txt = to_ron_string(&s).unwrap();
        assert!(
            from_ron_str(&txt).is_err(),
            "out-of-range constraint ref must be rejected on load"
        );
    }

    #[test]
    fn from_ron_rejects_entity_with_out_of_range_var() {
        // A `Point3` whose variable indices exceed the (empty) `vars`
        // vector is structurally corrupt; `validate` rejects it on load.
        // (Were a constraint to reference this point, `Point3::read`'s
        // `vars[x_var]` would then panic on solve.)
        let mut s = Sketch3D::new();
        s.entities
            .push(crate::entity::Entity3D::Point3(crate::entity::Point3 {
                x_var: 0,
                y_var: 1,
                z_var: 2,
            }));
        let txt = to_ron_string(&s).unwrap();
        assert!(
            from_ron_str(&txt).is_err(),
            "out-of-range variable index must be rejected on load"
        );
    }

    #[test]
    fn from_ron_round_trips_valid_mixed_constraints() {
        // The too-strict direction: a VALID sketch exercising the
        // discrimination-sensitive constraints (PointInPlane -> Plane3,
        // OnPlane -> Workplane, PlaneFixed -> Plane3-or-Workplane,
        // CircleRadius -> Circle3, ArcRadius -> Arc3) must still load --
        // guards against `validate` demanding the wrong kind and rejecting
        // a legitimate file.
        let mut s = Sketch3D::new();
        let o = s.add_point(0.0, 0.0, 0.0);
        let p = s.add_point(0.0, 0.0, 1.0);
        let plane = s.add_plane(o, 0.0, 0.0, 1.0).unwrap();
        let wp = s.add_workplane(o, 0.0, 0.0, 1.0).unwrap();
        let center = s.add_point(1.0, 0.0, 0.0);
        let circle = s.add_circle(center, 0.5, 0.0, 0.0, 1.0).unwrap();
        let astart = s.add_point(1.5, 0.0, 0.0);
        let aend = s.add_point(1.0, 0.5, 0.0);
        let arc = s.add_arc(center, 0.5, 0.0, 0.0, 1.0, astart, aend).unwrap();
        s.add_constraint(Constraint3D::PointInPlane { point: p, plane });
        s.add_constraint(Constraint3D::OnPlane {
            point: p,
            workplane: wp,
        });
        s.lock_plane(plane).unwrap(); // adds a PlaneFixed constraint
        s.add_constraint(Constraint3D::CircleRadius {
            circle,
            target: 0.5,
        });
        s.add_constraint(Constraint3D::ArcRadius { arc, target: 0.5 });
        let txt = to_ron_string(&s).unwrap();
        let f = from_ron_str(&txt).expect("a valid mixed-constraint sketch must load");
        assert!(f.sketch.validate().is_ok());
    }

    #[test]
    fn from_ron_rejects_wrong_kind_constraint_ref() {
        // A constraint that wants a Line but is handed a Point: the
        // accessor's `other => panic!()` arm fires pre-fix. `validate`
        // catches the kind mismatch on load.
        let mut s = Sketch3D::new();
        let _p = s.add_point(0.0, 0.0, 0.0);
        // LineLength3 references entity 0, which is a Point3, not a Line3.
        s.constraints.push(Constraint3D::LineLength3 {
            line: crate::entity::EntityId(0),
            target: 1.0,
        });
        let txt = to_ron_string(&s).unwrap();
        assert!(
            from_ron_str(&txt).is_err(),
            "wrong-kind constraint ref must be rejected on load"
        );
    }

    #[test]
    fn round_trip_feature_timeline() {
        let mut tl = FeatureTimeline::new();
        tl.push(Step::at_origin(
            Op::New,
            Feature::Box {
                dx: "size".into(),
                dy: "size".into(),
                dz: "size".into(),
            },
        ));
        tl.push(Step::placed(
            Op::Cut,
            Feature::Cylinder {
                radius: "hole_r".into(),
                height: "hole_h".into(),
            },
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
        assert!(
            model.bodies[0].faces() > 6,
            "round-tripped tree still punches a hole"
        );
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
