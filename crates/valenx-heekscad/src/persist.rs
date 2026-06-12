//! RON-based persistence — workbench-internal format.

use serde::{Deserialize, Serialize};

use crate::drawing::Drawing;
use crate::error::HeeksCadError;

/// Schema version.
pub const VERSION: u32 = 1;

/// RON envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PanelFile {
    /// Schema version.
    pub version: u32,
    /// Active drawing.
    pub drawing: Drawing,
}

/// Serialise.
pub fn to_ron_string(d: &Drawing) -> Result<String, HeeksCadError> {
    let f = PanelFile {
        version: VERSION,
        drawing: d.clone(),
    };
    ron::ser::to_string(&f).map_err(|e| HeeksCadError::Persist(e.to_string()))
}

/// Deserialise.
pub fn from_ron_str(s: &str) -> Result<PanelFile, HeeksCadError> {
    let f: PanelFile = ron::from_str(s).map_err(|e| HeeksCadError::Persist(e.to_string()))?;
    if f.version != VERSION {
        return Err(HeeksCadError::Persist(format!(
            "unsupported version {} expected {VERSION}",
            f.version
        )));
    }
    Ok(f)
}
