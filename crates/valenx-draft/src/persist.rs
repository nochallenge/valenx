//! RON-based persistence for draft documents.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::document::DraftDocument;
use crate::error::DraftError;

/// On-disk envelope wrapping a document with format-version metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DraftFile {
    /// Format version — bumped when on-disk schema changes.
    pub version: u32,
    /// The document payload.
    pub document: DraftDocument,
}

impl DraftFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a document with the current version tag.
    pub fn from_document(document: &DraftDocument) -> Self {
        Self {
            version: Self::VERSION,
            document: document.clone(),
        }
    }

    /// Serialize to a pretty RON string.
    pub fn to_ron(&self) -> Result<String, DraftError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| DraftError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), DraftError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, DraftError> {
        ron::from_str(s).map_err(|e| DraftError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// R29 D: bounded at [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`]
    /// (16 MiB) via the canonical helper, replacing the per-crate
    /// private `read_capped_to_string` duplicate.
    pub fn read_from(path: &Path) -> Result<Self, DraftError> {
        let s = valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?;
        Self::from_ron(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::DraftEntity;
    use crate::plane::WorkingPlane;

    #[test]
    fn round_trips_empty_document() {
        let d = DraftDocument::new(WorkingPlane::from_xy());
        let f = DraftFile::from_document(&d);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = DraftFile::from_ron(&ron).unwrap();
        assert_eq!(back.version, 1);
        assert_eq!(back.document.entity_count(), 0);
    }

    #[test]
    fn writes_to_file() {
        let d = DraftDocument::new(WorkingPlane::from_xy());
        let f = DraftFile::from_document(&d);
        let tmp = std::env::temp_dir().join("valenx_draft_test_write.ron");
        f.write_to(&tmp).unwrap();
        assert!(tmp.exists());
        let _ = std::fs::remove_file(&tmp);
    }

    /// Task 13 — exercise every entity variant through the full
    /// save→reload→compare cycle.
    #[test]
    fn round_trip_all_entity_variants() {
        let mut d = DraftDocument::new(WorkingPlane::from_xz());
        d.add_entity(DraftEntity::Line {
            start: [0.0, 0.0],
            end: [1.0, 1.0],
        });
        d.add_entity(DraftEntity::Polyline {
            points: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            closed: true,
        });
        d.add_entity(DraftEntity::Arc {
            center: [0.0, 0.0],
            radius: 2.0,
            start_angle: 0.0,
            end_angle: std::f64::consts::PI,
        });
        d.add_entity(DraftEntity::Circle {
            center: [3.0, 4.0],
            radius: 1.5,
        });
        d.add_entity(DraftEntity::Rectangle {
            min: [-1.0, -2.0],
            max: [3.0, 4.0],
        });
        d.add_entity(DraftEntity::Polygon {
            center: [0.0, 0.0],
            radius: 1.0,
            sides: 6,
        });
        d.add_entity(DraftEntity::LinearDimension {
            from: [0.0, 0.0],
            to: [10.0, 0.0],
            offset: 1.5,
        });
        d.add_entity(DraftEntity::Text {
            position: [1.0, 2.0],
            content: "label".to_string(),
            size: 0.5,
        });

        let ron = DraftFile::from_document(&d).to_ron().unwrap();
        let parsed = DraftFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.document, d);
    }

    #[test]
    fn round_trip_via_disk() {
        let mut d = DraftDocument::new(WorkingPlane::from_yz());
        d.add_entity(DraftEntity::Circle {
            center: [0.0, 0.0],
            radius: 1.0,
        });
        let tmp = std::env::temp_dir().join("valenx_draft_test_disk.ron");
        DraftFile::from_document(&d).write_to(&tmp).unwrap();
        let back = DraftFile::read_from(&tmp).unwrap();
        assert_eq!(back.document, d);
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn bad_ron_returns_ron_error() {
        let bad = "not a ron document";
        match DraftFile::from_ron(bad).unwrap_err() {
            DraftError::Ron(_) => {}
            other => panic!("wrong variant: {other:?}"),
        }
    }
}
