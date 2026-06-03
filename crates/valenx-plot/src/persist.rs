//! RON-based persistence for plots.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::PlotError;
use crate::plot::Plot;

/// On-disk envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlotFile {
    /// Format version.
    pub version: u32,
    /// The plot payload.
    pub plot: Plot,
}

impl PlotFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a plot.
    pub fn from_plot(plot: &Plot) -> Self {
        Self {
            version: Self::VERSION,
            plot: plot.clone(),
        }
    }

    /// Pretty RON.
    pub fn to_ron(&self) -> Result<String, PlotError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| PlotError::Ron(e.to_string()))
    }

    /// Write to disk. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), PlotError> {
        let s = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &s)?;
        Ok(())
    }

    /// Parse RON.
    pub fn from_ron(s: &str) -> Result<Self, PlotError> {
        ron::from_str(s).map_err(|e| PlotError::Ron(e.to_string()))
    }

    /// Read from disk.
    ///
    /// Round-23 workspace sweep: bounded at
    /// [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`] (16 MiB). Sister
    /// to the round-12 persist.rs family.
    pub fn read_from(path: &Path) -> Result<Self, PlotError> {
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
    use crate::series::{Series, SeriesStyle};

    #[test]
    fn round_trip() {
        let mut p = Plot::new("Demo");
        let mut s = Series::new("a", SeriesStyle::Line);
        s.push(0.0, 0.0);
        s.push(1.0, 1.0);
        p.add_series(s);
        let f = PlotFile::from_plot(&p);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = PlotFile::from_ron(&ron).unwrap();
        assert_eq!(back.version, 1);
        assert_eq!(back.plot.title, "Demo");
    }

    /// Round-23 RED→GREEN: `read_from` rejects oversize `.ron`.
    #[test]
    fn read_from_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_plot_oversize_test.ron");
        let oversize = (valenx_core::io_caps::MAX_DOC_FILE_BYTES) + 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = PlotFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            PlotError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
