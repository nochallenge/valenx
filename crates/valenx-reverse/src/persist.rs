//! RON-based persistence for point clouds.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::ReverseError;
use crate::pointcloud::PointCloud;

/// On-disk envelope wrapping a [`PointCloud`] with format-version
/// metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReverseFile {
    /// Format version.
    pub version: u32,
    /// The cloud payload.
    pub cloud: PointCloud,
}

impl ReverseFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a cloud with the current version tag.
    pub fn from_cloud(cloud: &PointCloud) -> Self {
        Self {
            version: Self::VERSION,
            cloud: cloud.clone(),
        }
    }

    /// Serialize to a pretty RON string.
    pub fn to_ron(&self) -> Result<String, ReverseError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| ReverseError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), ReverseError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, ReverseError> {
        ron::from_str(s).map_err(|e| ReverseError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// Round-23 workspace sweep: bounded at
    /// [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`] (16 MiB). Sister
    /// to the round-12 persist.rs family.
    pub fn read_from(path: &Path) -> Result<Self, ReverseError> {
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
    use nalgebra::Vector3;

    #[test]
    fn round_trip_empty_cloud() {
        let c = PointCloud::new();
        let ron = ReverseFile::from_cloud(&c).to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = ReverseFile::from_ron(&ron).unwrap();
        assert_eq!(back.version, 1);
        assert!(back.cloud.is_empty());
    }

    #[test]
    fn round_trip_with_normals() {
        let mut c = PointCloud::from_points(vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)]);
        c.normals = Some(vec![Vector3::z(), Vector3::z()]);
        let ron = ReverseFile::from_cloud(&c).to_ron().unwrap();
        let back = ReverseFile::from_ron(&ron).unwrap();
        assert_eq!(back.cloud, c);
    }

    /// Round-23 RED→GREEN: `read_from` rejects oversize `.ron`.
    #[test]
    fn read_from_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_reverse_oversize_test.ron");
        let oversize = (valenx_core::io_caps::MAX_DOC_FILE_BYTES) + 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = ReverseFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            ReverseError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
