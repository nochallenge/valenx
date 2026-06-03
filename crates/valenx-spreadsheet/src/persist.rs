//! RON envelope for round-tripping spreadsheet workbooks.
//!
//! Same pattern as `valenx-sketch::persist::SketchFile`,
//! `valenx-surface::persist::SurfaceFile`, etc.: a thin envelope
//! wrapping the live in-memory [`Spreadsheet`] with a `version` field
//! so we can evolve the schema without breaking older files.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::SpreadsheetError;
use crate::sheet::Spreadsheet;

/// On-disk envelope wrapping a spreadsheet workbook.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SpreadsheetFile {
    /// Format version - bumped when the on-disk schema changes.
    pub version: u32,
    /// Live workbook payload.
    pub spreadsheet: Spreadsheet,
}

impl SpreadsheetFile {
    /// Current on-disk format version. Bump and add a migration step
    /// in `from_ron` when the schema changes.
    pub const VERSION: u32 = 1;

    /// Construct an empty spreadsheet file at the current `VERSION`.
    pub fn new() -> Self {
        Self {
            version: Self::VERSION,
            spreadsheet: Spreadsheet::new(),
        }
    }

    /// Wrap an existing [`Spreadsheet`] at the current `VERSION`.
    pub fn from_spreadsheet(spreadsheet: Spreadsheet) -> Self {
        Self {
            version: Self::VERSION,
            spreadsheet,
        }
    }

    /// Serialize to a pretty-printed RON string.
    ///
    /// # Errors
    ///
    /// Returns [`SpreadsheetError::Ron`] if the underlying RON
    /// serializer fails (typically only when the model contains
    /// non-finite floats; we don't try to filter those here).
    pub fn to_ron(&self) -> Result<String, SpreadsheetError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| SpreadsheetError::Ron(e.to_string()))
    }

    /// Write to a file. Overwrites if the file exists. Round-28 H2:
    /// routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    ///
    /// # Errors
    ///
    /// IO error wrapping `std::io::Error`, or [`SpreadsheetError::Ron`]
    /// from serialization.
    pub fn write_to(&self, path: &Path) -> Result<(), SpreadsheetError> {
        valenx_core::io_caps::atomic_write_str(path, &self.to_ron()?)?;
        Ok(())
    }

    /// Parse from a RON string.
    ///
    /// # Errors
    ///
    /// Returns [`SpreadsheetError::Ron`] when the input is not valid
    /// RON or doesn't deserialize into the current schema.
    pub fn from_ron(s: &str) -> Result<Self, SpreadsheetError> {
        ron::from_str(s).map_err(|e| SpreadsheetError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// # Errors
    ///
    /// IO error wrapping `std::io::Error`, or [`SpreadsheetError::Ron`]
    /// from deserialization.
    pub fn read_from(path: &Path) -> Result<Self, SpreadsheetError> {
        // R29 D: canonical valenx_core::io_caps::read_capped_to_string at
        // MAX_DOC_FILE_BYTES (16 MiB), replacing the private dupe.
        Self::from_ron(&valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{Cell, CellRef};

    #[test]
    fn round_trips_empty() {
        let f = SpreadsheetFile::new();
        let ron = f.to_ron().unwrap();
        let parsed = SpreadsheetFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.version, 1);
        assert!(parsed.spreadsheet.sheets.is_empty());
    }

    #[test]
    fn round_trips_workbook_with_cells() {
        let mut ss = Spreadsheet::new();
        ss.add_sheet("Default");
        ss.set_cell(&CellRef::parse("Default.A1").unwrap(), Cell::Number(50.0))
            .unwrap();
        ss.set_cell(
            &CellRef::parse("Default.A2").unwrap(),
            Cell::Formula("Default.A1 * 2".into()),
        )
        .unwrap();
        ss.set_cell(
            &CellRef::parse("Default.B1").unwrap(),
            Cell::Text("wallHeight".into()),
        )
        .unwrap();

        let f = SpreadsheetFile::from_spreadsheet(ss);
        let ron = f.to_ron().unwrap();
        let parsed = SpreadsheetFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.spreadsheet.sheets.len(), 1);
        let a1 = parsed
            .spreadsheet
            .cell(&CellRef::parse("Default.A1").unwrap());
        assert!(matches!(a1, Cell::Number(n) if (*n - 50.0).abs() < 1e-12));
        let a2 = parsed
            .spreadsheet
            .cell(&CellRef::parse("Default.A2").unwrap());
        assert!(matches!(a2, Cell::Formula(s) if s == "Default.A1 * 2"));
        let b1 = parsed
            .spreadsheet
            .cell(&CellRef::parse("Default.B1").unwrap());
        assert!(matches!(b1, Cell::Text(s) if s == "wallHeight"));
    }

    #[test]
    fn writes_to_file_round_trips() {
        let mut ss = Spreadsheet::new();
        ss.add_sheet("S");
        ss.set_cell(&CellRef::parse("S.A1").unwrap(), Cell::Number(7.0))
            .unwrap();
        let f = SpreadsheetFile::from_spreadsheet(ss);

        let tmp = std::env::temp_dir().join("valenx_spreadsheet_persist_test.ron");
        f.write_to(&tmp).unwrap();
        let parsed = SpreadsheetFile::read_from(&tmp).unwrap();
        assert_eq!(parsed.spreadsheet.sheets.len(), 1);
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn malformed_ron_errors() {
        let e = SpreadsheetFile::from_ron("not_valid_ron").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.ron");
    }
}
