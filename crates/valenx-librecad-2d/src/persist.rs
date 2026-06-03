//! Versioned RON persistence (parallel to DXF — RON is the
//! workbench-internal format, DXF is the interchange format).

use serde::{Deserialize, Serialize};

use crate::drawing::Drawing2D;
use crate::error::LibreCadError;

/// Schema version.
pub const VERSION: u32 = 1;

/// RON envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PanelFile {
    /// Schema version.
    pub version: u32,
    /// Active drawing.
    pub drawing: Drawing2D,
}

/// Serialise.
pub fn to_ron_string(d: &Drawing2D) -> Result<String, LibreCadError> {
    let f = PanelFile {
        version: VERSION,
        drawing: d.clone(),
    };
    ron::ser::to_string(&f).map_err(|e| LibreCadError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })
}

/// Deserialise.
pub fn from_ron_str(s: &str) -> Result<PanelFile, LibreCadError> {
    let f: PanelFile = ron::from_str(s).map_err(|e| LibreCadError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })?;
    if f.version != VERSION {
        return Err(LibreCadError::BadParameter {
            name: "version",
            reason: format!("unsupported {} expected {VERSION}", f.version),
        });
    }
    Ok(f)
}
