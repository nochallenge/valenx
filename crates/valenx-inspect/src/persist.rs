//! RON-based persistence for inspection reports.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::InspectError;
use crate::report::InspectReport;

/// On-disk envelope wrapping a report with format-version metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InspectFile {
    /// Format version — bumped when on-disk schema changes.
    pub version: u32,
    /// The report payload.
    pub report: InspectReport,
}

impl InspectFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a report with the current version tag.
    pub fn from_report(report: &InspectReport) -> Self {
        Self {
            version: Self::VERSION,
            report: report.clone(),
        }
    }

    /// Serialize to a pretty RON string.
    pub fn to_ron(&self) -> Result<String, InspectError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| InspectError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), InspectError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, InspectError> {
        ron::from_str(s).map_err(|e| InspectError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// Round-23 workspace sweep: bounded at
    /// [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`] (16 MiB). Sister
    /// to the round-12 persist.rs family.
    pub fn read_from(path: &Path) -> Result<Self, InspectError> {
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
    use crate::measurement::Measurement;
    use crate::report::{CheckResult, InspectReport};
    use crate::tolerance::Tolerance;
    use nalgebra::Vector3;

    #[test]
    fn round_trip_empty() {
        let r = InspectReport::with_title("Empty");
        let f = InspectFile::from_report(&r);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = InspectFile::from_ron(&ron).unwrap();
        assert_eq!(back.report.title, "Empty");
        assert!(back.report.rows.is_empty());
    }

    #[test]
    fn round_trip_with_rows() {
        let mut r = InspectReport::with_title("Run 1");
        r.add_labelled(
            "AB",
            Measurement::Distance {
                from: Vector3::zeros(),
                to: Vector3::new(5.0, 0.0, 0.0),
            },
            Tolerance::symmetric(5.0, 0.1),
            CheckResult::Pass,
        );
        let ron = InspectFile::from_report(&r).to_ron().unwrap();
        let back = InspectFile::from_ron(&ron).unwrap();
        assert_eq!(back.report, r);
    }

    #[test]
    fn bad_ron_errors() {
        match InspectFile::from_ron("nonsense").unwrap_err() {
            InspectError::Ron(_) => {}
            other => panic!("wrong variant: {other:?}"),
        }
    }

    /// Round-23 RED→GREEN: `read_from` rejects oversize `.ron`.
    #[test]
    fn read_from_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_inspect_oversize_test.ron");
        let oversize = (valenx_core::io_caps::MAX_DOC_FILE_BYTES) + 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = InspectFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            InspectError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
