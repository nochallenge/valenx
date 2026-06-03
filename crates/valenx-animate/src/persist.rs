//! RON-based persistence for an animation.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::animation::Animation;
use crate::error::AnimateError;

/// On-disk envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AnimateFile {
    /// Format version.
    pub version: u32,
    /// The animation payload.
    pub animation: Animation,
}

impl AnimateFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap an animation with the current version tag.
    pub fn from_animation(animation: &Animation) -> Self {
        Self {
            version: Self::VERSION,
            animation: animation.clone(),
        }
    }

    /// Pretty RON.
    pub fn to_ron(&self) -> Result<String, AnimateError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| AnimateError::Ron(e.to_string()))
    }

    /// Write to disk. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), AnimateError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse RON.
    pub fn from_ron(s: &str) -> Result<Self, AnimateError> {
        ron::from_str(s).map_err(|e| AnimateError::Ron(e.to_string()))
    }

    /// Read from disk.
    ///
    /// Round-23 workspace sweep: the read is bounded at
    /// [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`] (16 MiB) — sister
    /// to the round-12 persist.rs sweep that capped sketch / cam /
    /// arch / techdraw etc. Animations are a few dozen keyframes per
    /// channel; 16 MiB is well past anything realistic while refusing
    /// hostile multi-GB files that would OOM `String::from_utf8`.
    pub fn read_from(path: &Path) -> Result<Self, AnimateError> {
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
    use crate::animation::Keyframe;

    #[test]
    fn round_trip_empty() {
        let a = Animation::new();
        let f = AnimateFile::from_animation(&a);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = AnimateFile::from_ron(&ron).unwrap();
        assert_eq!(back.animation, a);
    }

    #[test]
    fn round_trip_with_keyframes() {
        let mut a = Animation::new();
        a.name = "test".into();
        a.push(Keyframe::at(0.0).with_joint(0, 0.0)).unwrap();
        a.push(Keyframe::at(2.0).with_joint(0, 90.0)).unwrap();
        let ron = AnimateFile::from_animation(&a).to_ron().unwrap();
        let back = AnimateFile::from_ron(&ron).unwrap();
        assert_eq!(back.animation, a);
    }

    /// Round-23 RED→GREEN: `read_from` rejects oversize `.ron`.
    /// Sister to the round-12 persist.rs family.
    #[test]
    fn read_from_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_animate_oversize_test.ron");
        let oversize = (valenx_core::io_caps::MAX_DOC_FILE_BYTES) + 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = AnimateFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            AnimateError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
