//! RON-based persistence for [`crate::ArchDocument`].

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::document::ArchDocument;
use crate::error::ArchError;

/// On-disk envelope wrapping an [`ArchDocument`] with format-version
/// metadata. Same shape as the persistence files in `valenx-draft`,
/// `valenx-techdraw`, etc.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArchFile {
    /// Format version — bumped when the on-disk schema changes.
    pub version: u32,
    /// The document payload.
    pub document: ArchDocument,
}

impl ArchFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a document with the current version tag.
    pub fn from_document(document: &ArchDocument) -> Self {
        Self {
            version: Self::VERSION,
            document: document.clone(),
        }
    }

    /// Serialize to a pretty RON string.
    pub fn to_ron(&self) -> Result<String, ArchError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| ArchError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic. The R27 commit claimed this was
    /// already a wrapper but the grep showed otherwise.
    pub fn write_to(&self, path: &Path) -> Result<(), ArchError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, ArchError> {
        ron::from_str(s).map_err(|e| ArchError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// R29 D: bounded at [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`]
    /// (16 MiB) via the canonical helper, replacing the per-crate
    /// private `read_capped_to_string` duplicate.
    pub fn read_from(path: &Path) -> Result<Self, ArchError> {
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
    use crate::beam::{BeamParams, BeamSection};
    use crate::column::{ColumnParams, ColumnSection};
    use crate::door::{DoorParams, DoorStyle, Side};
    use crate::entity::ArchEntity;
    use crate::roof::{RoofParams, RoofType};
    use crate::slab::SlabParams;
    use crate::space::SpaceParams;
    use crate::stair::StairParams;
    use crate::wall::WallParams;
    use crate::window::{WindowParams, WindowStyle};
    use nalgebra::Vector3;

    fn full_doc() -> ArchDocument {
        let mut d = ArchDocument::new("test-project");
        d.add_entity(ArchEntity::Wall(WallParams {
            start: Vector3::zeros(),
            end: Vector3::new(3.0, 0.0, 0.0),
            height: 2.5,
            thickness: 0.2,
            material: "Brick".into(),
        }));
        d.add_entity(ArchEntity::Slab(SlabParams {
            boundary: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(3.0, 3.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
            ],
            thickness: 0.2,
            material: "Concrete".into(),
            structural: None,
        }));
        d.add_entity(ArchEntity::Column(ColumnParams {
            base: Vector3::new(1.0, 1.0, 0.0),
            height: 2.5,
            cross_section: ColumnSection::Circular {
                radius: 0.15,
                segments: 12,
            },
            material: "Steel".into(),
            structural: None,
        }));
        d.add_entity(ArchEntity::Beam(BeamParams {
            start: Vector3::new(0.0, 0.0, 2.5),
            end: Vector3::new(3.0, 0.0, 2.5),
            cross_section: BeamSection::IBeam {
                width: 0.2,
                depth: 0.3,
                flange_thickness: 0.02,
                web_thickness: 0.01,
            },
            orientation_angle: 0.0,
            material: "Steel".into(),
            structural: None,
        }));
        d.add_entity(ArchEntity::Window(WindowParams {
            host: 1,
            position_along_wall: 1.5,
            position_height: 1.0,
            width: 0.8,
            height: 1.0,
            frame_thickness: 0.05,
            style: WindowStyle::Casement,
        }));
        d.add_entity(ArchEntity::Door(DoorParams {
            host: 1,
            position_along_wall: 2.5,
            width: 0.9,
            height: 2.1,
            style: DoorStyle::Single,
            hinge_side: Side::Left,
        }));
        d.add_entity(ArchEntity::Stair(StairParams {
            base: Vector3::zeros(),
            direction: Vector3::new(1.0, 0.0, 0.0),
            total_rise: 3.0,
            total_run: 4.0,
            num_steps: 12,
            width: 1.0,
        }));
        d.add_entity(ArchEntity::Roof(RoofParams {
            boundary: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(3.0, 3.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
            ],
            peak_height: 1.5,
            roof_type: RoofType::Gable,
        }));
        d.add_entity(ArchEntity::Space(SpaceParams {
            boundary: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(3.0, 3.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
            ],
            ceiling_height: 2.5,
            space_name: "Living".into(),
        }));
        d
    }

    #[test]
    fn round_trip_empty_document() {
        let d = ArchDocument::new("p");
        let f = ArchFile::from_document(&d);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = ArchFile::from_ron(&ron).unwrap();
        assert_eq!(back.document.count(), 0);
    }

    #[test]
    fn round_trip_full_document_via_string() {
        let d = full_doc();
        let ron = ArchFile::from_document(&d).to_ron().unwrap();
        let back = ArchFile::from_ron(&ron).unwrap();
        assert_eq!(back.document, d);
    }

    #[test]
    fn round_trip_full_document_via_disk() {
        let d = full_doc();
        let tmp = std::env::temp_dir().join("valenx_arch_persist_test.ron");
        ArchFile::from_document(&d).write_to(&tmp).unwrap();
        let back = ArchFile::read_from(&tmp).unwrap();
        assert_eq!(back.document, d);
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn bad_ron_returns_ron_error() {
        let bad = "not a ron document";
        match ArchFile::from_ron(bad).unwrap_err() {
            ArchError::Ron(_) => {}
            other => panic!("wrong variant: {other:?}"),
        }
    }
}
