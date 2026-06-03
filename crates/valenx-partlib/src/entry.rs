//! [`PartEntry`] — a single installed library item.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

/// The format of the underlying part file.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum PartKind {
    /// STEP CAD exchange file (`.step`, `.stp`).
    StepFile,
    /// IGES CAD exchange file (`.iges`, `.igs`).
    IgesFile,
    /// Stereolithography mesh (`.stl`).
    StlMesh,
    /// Native NURBS surface dump (Valenx-specific RON, used by
    /// valenx-surface).
    NurbsSurface,
}

impl PartKind {
    /// Best-effort guess from a file extension. Returns `None` when
    /// the extension doesn't match a known kind.
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_ascii_lowercase().as_str() {
            "step" | "stp" => Some(Self::StepFile),
            "iges" | "igs" => Some(Self::IgesFile),
            "stl" => Some(Self::StlMesh),
            "ron" => Some(Self::NurbsSurface),
            _ => None,
        }
    }
}

/// One installed part in the library.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PartEntry {
    /// Library-unique name (typically the source filename without
    /// extension).
    pub name: String,
    /// File format.
    pub kind: PartKind,
    /// Original source URL, if any (filled in by network-fetched
    /// installs once Phase 46 v2 lands).
    pub source_url: Option<String>,
    /// Local path on disk where the part lives.
    pub local_path: PathBuf,
    /// SHA-256 hex digest of the on-disk file at install time.
    pub checksum: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kind_from_extension_matches_common_cases() {
        assert_eq!(PartKind::from_extension("step"), Some(PartKind::StepFile));
        assert_eq!(PartKind::from_extension("STL"), Some(PartKind::StlMesh));
        assert_eq!(PartKind::from_extension("igs"), Some(PartKind::IgesFile));
        assert_eq!(PartKind::from_extension("ron"), Some(PartKind::NurbsSurface));
        assert_eq!(PartKind::from_extension("docx"), None);
    }
}
