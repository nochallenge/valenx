//! Versioned RON persistence for the panel state.

use serde::{Deserialize, Serialize};

use crate::animation::Animation;
use crate::error::CamoticsError;

/// Schema version.
pub const VERSION: u32 = 1;

/// RON envelope. Only the animation spec is persisted; the report is
/// re-computed on load.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PanelFile {
    /// Schema version.
    pub version: u32,
    /// Optional animation.
    pub animation: Option<Animation>,
    /// Playback position.
    pub current_t: f64,
}

/// Serialise.
pub fn to_ron_string(
    animation: &Option<Animation>,
    current_t: f64,
) -> Result<String, CamoticsError> {
    let f = PanelFile {
        version: VERSION,
        animation: animation.clone(),
        current_t,
    };
    ron::ser::to_string(&f).map_err(|e| CamoticsError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })
}

/// Deserialise.
pub fn from_ron_str(s: &str) -> Result<PanelFile, CamoticsError> {
    let f: PanelFile = ron::from_str(s).map_err(|e| CamoticsError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })?;
    if f.version != VERSION {
        return Err(CamoticsError::BadParameter {
            name: "version",
            reason: format!("unsupported {} expected {VERSION}", f.version),
        });
    }
    Ok(f)
}
