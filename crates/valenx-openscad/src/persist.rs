//! Versioned RON persistence for the OpenSCAD panel state.

use serde::{Deserialize, Serialize};

use crate::error::OpenScadCsgError;

/// Persistence schema version.
pub const VERSION: u32 = 1;

/// RON envelope for the OpenSCAD panel — only the source text is
/// persisted; the cached AST + Solid are regenerated on load.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PanelFile {
    /// Schema version.
    pub version: u32,
    /// Current editor source.
    pub editor_text: String,
}

/// Serialise to a RON string.
pub fn to_ron_string(text: &str) -> Result<String, OpenScadCsgError> {
    let f = PanelFile {
        version: VERSION,
        editor_text: text.to_string(),
    };
    ron::ser::to_string(&f).map_err(|e| OpenScadCsgError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })
}

/// Deserialise from a RON string.
pub fn from_ron_str(s: &str) -> Result<PanelFile, OpenScadCsgError> {
    let f: PanelFile = ron::from_str(s).map_err(|e| OpenScadCsgError::BadParameter {
        name: "ron",
        reason: e.to_string(),
    })?;
    if f.version != VERSION {
        return Err(OpenScadCsgError::BadParameter {
            name: "version",
            reason: format!("unsupported version {}, expected {VERSION}", f.version),
        });
    }
    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        let src = "cube([1, 1, 1]);";
        let s = to_ron_string(src).expect("ok");
        let f = from_ron_str(&s).expect("ok");
        assert_eq!(f.editor_text, src);
        assert_eq!(f.version, VERSION);
    }
}
