//! RON-based persistence for lattice recipes.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::LatticeError;
use crate::lattice::Lattice;

/// On-disk envelope wrapping a [`Lattice`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LatticeFile {
    /// Format version.
    pub version: u32,
    /// The lattice recipe.
    pub lattice: Lattice,
}

impl LatticeFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a recipe.
    pub fn from_lattice(lattice: Lattice) -> Self {
        Self {
            version: Self::VERSION,
            lattice,
        }
    }

    /// Pretty RON.
    pub fn to_ron(&self) -> Result<String, LatticeError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| LatticeError::Ron(e.to_string()))
    }

    /// Write to disk. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), LatticeError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse RON.
    pub fn from_ron(s: &str) -> Result<Self, LatticeError> {
        ron::from_str(s).map_err(|e| LatticeError::Ron(e.to_string()))
    }

    /// Read from disk.
    ///
    /// R29 D: bounded at [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`]
    /// (16 MiB) via the canonical helper, replacing the per-crate
    /// private `read_capped_to_string` duplicate.
    pub fn read_from(path: &Path) -> Result<Self, LatticeError> {
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
    fn round_trip_grid() {
        let l = Lattice::Grid {
            rows: 3,
            cols: 4,
            levels: 1,
            spacing: Vector3::new(1.0, 1.0, 1.0),
        };
        let f = LatticeFile::from_lattice(l);
        let ron = f.to_ron().unwrap();
        let back = LatticeFile::from_ron(&ron).unwrap();
        assert_eq!(back.version, 1);
        match back.lattice {
            Lattice::Grid {
                rows, cols, levels, ..
            } => {
                assert_eq!((rows, cols, levels), (3, 4, 1));
            }
            _ => panic!("wrong variant"),
        }
    }

    /// R29 E: exercise the actual disk `write_to` → `read_from` path
    /// (the R28 atomic_write migration was previously covered only by
    /// the in-memory RON string round-trip). Confirms the sidecar +
    /// rename write publishes bytes that round-trip back to an equal
    /// recipe.
    #[test]
    fn write_to_read_from_round_trips_through_disk() {
        let l = Lattice::Grid {
            rows: 2,
            cols: 5,
            levels: 3,
            spacing: Vector3::new(1.5, 2.5, 0.5),
        };
        let f = LatticeFile::from_lattice(l);
        let tmp = std::env::temp_dir().join(format!(
            "valenx_lattice_disk_rt_{}.ron",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        f.write_to(&tmp).expect("write_to should publish the file");
        let back = LatticeFile::read_from(&tmp).expect("read_from should parse the written file");
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(back.version, LatticeFile::VERSION);
        match back.lattice {
            Lattice::Grid {
                rows,
                cols,
                levels,
                spacing,
            } => {
                assert_eq!((rows, cols, levels), (2, 5, 3));
                assert!((spacing - Vector3::new(1.5, 2.5, 0.5)).norm() < 1e-9);
            }
            other => panic!("wrong variant after disk round-trip: {other:?}"),
        }
    }
}
