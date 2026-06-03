//! Versioned RON envelope for the OpenSCAD-import panel state.

use serde::{Deserialize, Serialize};

use crate::error::OpenScadError;
use crate::panel::OpenScadImportPanelState;

/// File format version. Bumped on schema changes.
pub const VERSION: u32 = 1;

/// On-wire envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PanelFile {
    /// Format version.
    pub version: u32,
    /// Wrapped panel state.
    pub panel: OpenScadImportPanelState,
}

impl PanelFile {
    /// Wrap panel state in the current envelope.
    pub fn new(panel: OpenScadImportPanelState) -> Self {
        Self {
            version: VERSION,
            panel,
        }
    }
}

/// Serialise to a pretty RON string.
pub fn to_ron_string(p: &OpenScadImportPanelState) -> Result<String, OpenScadError> {
    let file = PanelFile::new(p.clone());
    ron::ser::to_string_pretty(&file, ron::ser::PrettyConfig::default())
        .map_err(|e| OpenScadError::Cad(format!("ron: {e}")))
}

/// Parse a panel state from a RON string. Fails on version mismatch.
pub fn from_ron_str(s: &str) -> Result<OpenScadImportPanelState, OpenScadError> {
    let file: PanelFile = ron::de::from_str(s).map_err(|e| OpenScadError::Cad(format!("ron: {e}")))?;
    if file.version != VERSION {
        return Err(OpenScadError::Cad(format!(
            "version mismatch: file = {}, expected = {}",
            file.version, VERSION
        )));
    }
    Ok(file.panel)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips() {
        let mut p = OpenScadImportPanelState::new();
        p.set_status("ok", 99);
        p.file_path_input = "cube.scad".into();
        let text = to_ron_string(&p).expect("ser");
        let back = from_ron_str(&text).expect("de");
        assert_eq!(back.file_path_input, "cube.scad");
        assert_eq!(back.last_face_count, 99);
    }

    #[test]
    fn rejects_version_mismatch() {
        let bad = "(version: 99, panel: (file_path_input: \"\", last_status: None, last_error: None, last_face_count: 0))";
        assert!(matches!(from_ron_str(bad), Err(OpenScadError::Cad(_))));
    }
}
