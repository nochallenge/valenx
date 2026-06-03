//! RON-based persistence for a recorded list of curve operations.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::CurvesError;
use crate::ops::{DiscretizeMode, ExtendEnd};

/// One recorded curve operation — the data the UI panel will dispatch
/// against when the user re-runs the recipe.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum CurveOp {
    /// Offset a planar curve by `d`.
    Offset {
        /// Source curve index in the file's curve list.
        curve_idx: usize,
        /// Offset distance in mm.
        d: f64,
    },
    /// Discretize a curve.
    Discretize {
        /// Source curve index.
        curve_idx: usize,
        /// Sample count.
        n: usize,
        /// Sampling mode.
        mode: DiscretizeMode,
    },
    /// Reverse a curve's parameter direction.
    Reverse {
        /// Source curve index.
        curve_idx: usize,
    },
    /// Trim a curve to a sub-range.
    Trim {
        /// Source curve index.
        curve_idx: usize,
        /// Parameter start in [0, 1].
        t_start: f64,
        /// Parameter end in [0, 1].
        t_end: f64,
    },
    /// Extend a curve past one end.
    Extend {
        /// Source curve index.
        curve_idx: usize,
        /// Extension distance.
        length: f64,
        /// Which end to extend.
        end: ExtendEnd,
    },
}

/// On-disk envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CurvesFile {
    /// Format version.
    pub version: u32,
    /// Recorded operations.
    pub ops: Vec<CurveOp>,
}

impl CurvesFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap an op list with the current version tag.
    pub fn from_ops(ops: Vec<CurveOp>) -> Self {
        Self {
            version: Self::VERSION,
            ops,
        }
    }

    /// Serialize to a pretty RON string.
    pub fn to_ron(&self) -> Result<String, CurvesError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| CurvesError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), CurvesError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, CurvesError> {
        ron::from_str(s).map_err(|e| CurvesError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// Round-23 workspace sweep: bounded at
    /// [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`] (16 MiB). Sister
    /// to the round-12 persist.rs family.
    pub fn read_from(path: &Path) -> Result<Self, CurvesError> {
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

    #[test]
    fn round_trip_empty() {
        let f = CurvesFile::from_ops(vec![]);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = CurvesFile::from_ron(&ron).unwrap();
        assert_eq!(back.ops.len(), 0);
    }

    #[test]
    fn round_trip_with_ops() {
        let ops = vec![
            CurveOp::Reverse { curve_idx: 0 },
            CurveOp::Trim {
                curve_idx: 1,
                t_start: 0.1,
                t_end: 0.9,
            },
            CurveOp::Extend {
                curve_idx: 0,
                length: 1.0,
                end: ExtendEnd::End,
            },
        ];
        let f = CurvesFile::from_ops(ops);
        let ron = f.to_ron().unwrap();
        let back = CurvesFile::from_ron(&ron).unwrap();
        assert_eq!(back.ops.len(), 3);
    }

    /// Round-23 RED→GREEN: `read_from` rejects a `.ron` file larger
    /// than `MAX_DOC_FILE_BYTES` (16 MiB) at the read-cap layer
    /// rather than slurping it into memory. Sister to the round-12
    /// sketch persist test.
    #[test]
    fn read_from_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_curves_oversize_test.ron");
        let oversize = (valenx_core::io_caps::MAX_DOC_FILE_BYTES) + 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = CurvesFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            CurvesError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
