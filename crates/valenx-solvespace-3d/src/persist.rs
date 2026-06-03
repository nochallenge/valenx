//! Versioned RON persistence for the solver panel state.

use serde::{Deserialize, Serialize};

use crate::error::Solve3DError;
use crate::sketch::Sketch3D;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::Constraint3D;

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
}
