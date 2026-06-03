//! [`Sheet`] and [`Spreadsheet`] — the workbook data model.
//!
//! - A [`Sheet`] is a single named tab: a sparse map from
//!   `(row, col)` to [`Cell`] payload.
//! - A [`Spreadsheet`] is the whole workbook: a name-keyed map of
//!   sheets.
//!
//! Cell lookups return [`Cell::Empty`] for missing coordinates so
//! callers can pattern-match uniformly without special-casing
//! "doesn't exist yet".

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::cell::{Cell, CellRef};
use crate::error::SpreadsheetError;

/// One named sheet — `(row, col)` → cell payload.
///
/// Missing entries are [`Cell::Empty`] by convention; [`Sheet::set`]
/// erases the cell when given `Cell::Empty` instead of materialising
/// it, to keep the HashMap sparse.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Sheet {
    /// User-facing name (matches the key in [`Spreadsheet::sheets`]
    /// and the `Sheet.A1` reference prefix).
    pub name: String,
    /// Sparse `(row, col)` → cell payload map.
    pub cells: HashMap<(u32, u32), Cell>,
}

impl Sheet {
    /// Construct an empty sheet with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            cells: HashMap::new(),
        }
    }

    /// Set the cell at `(row, col)`. Inserting [`Cell::Empty`] erases
    /// the entry from the underlying map to keep it sparse.
    pub fn set(&mut self, row: u32, col: u32, cell: Cell) {
        match cell {
            Cell::Empty => {
                self.cells.remove(&(row, col));
            }
            other => {
                self.cells.insert((row, col), other);
            }
        }
    }

    /// Borrow the cell at `(row, col)`. Returns a reference to a
    /// shared [`Cell::Empty`] sentinel when the coordinate is empty.
    pub fn get(&self, row: u32, col: u32) -> &Cell {
        static EMPTY: Cell = Cell::Empty;
        self.cells.get(&(row, col)).unwrap_or(&EMPTY)
    }

    /// Iterate over all non-empty cells in arbitrary order.
    pub fn iter(&self) -> impl Iterator<Item = ((u32, u32), &Cell)> {
        self.cells.iter().map(|(&k, v)| (k, v))
    }

    /// Cell count (excludes empties because empties are never
    /// materialised in the HashMap).
    pub fn len(&self) -> usize {
        self.cells.len()
    }

    /// True when no cells have been set on this sheet.
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }
}

/// Workbook — name-keyed collection of sheets.
///
/// All cell mutations go through a [`CellRef`] so the sheet name and
/// `(row, col)` stay in sync.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Spreadsheet {
    /// All sheets, keyed by name. Mutations should go through
    /// [`Spreadsheet::add_sheet`] / [`Spreadsheet::set_cell`] so the
    /// inner [`Sheet::name`] stays in sync with the key.
    pub sheets: HashMap<String, Sheet>,
}

impl Spreadsheet {
    /// Empty workbook with no sheets.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a fresh sheet with the given name. Returns `true` when the
    /// sheet was created, `false` when a sheet with that name already
    /// existed (no-op in the latter case).
    pub fn add_sheet(&mut self, name: impl Into<String>) -> bool {
        let name = name.into();
        if self.sheets.contains_key(&name) {
            return false;
        }
        self.sheets.insert(name.clone(), Sheet::new(name));
        true
    }

    /// Remove a sheet by name. Returns `true` when something was
    /// removed.
    pub fn remove_sheet(&mut self, name: &str) -> bool {
        self.sheets.remove(name).is_some()
    }

    /// Borrow the cell at the given reference. Returns a reference to
    /// a shared [`Cell::Empty`] sentinel when the sheet or
    /// `(row, col)` is missing — callers can pattern-match uniformly.
    pub fn cell(&self, r: &CellRef) -> &Cell {
        static EMPTY: Cell = Cell::Empty;
        match self.sheets.get(&r.sheet_name) {
            Some(s) => s.get(r.row, r.col),
            None => &EMPTY,
        }
    }

    /// Set the cell at the given reference.
    ///
    /// # Errors
    ///
    /// Returns [`SpreadsheetError::BadCellRef`] when the referenced
    /// sheet doesn't exist. Callers should call [`Self::add_sheet`]
    /// first.
    pub fn set_cell(&mut self, r: &CellRef, cell: Cell) -> Result<(), SpreadsheetError> {
        let sheet = self.sheets.get_mut(&r.sheet_name).ok_or_else(|| {
            SpreadsheetError::BadCellRef(format!("no sheet named `{}`", r.sheet_name))
        })?;
        sheet.set(r.row, r.col, cell);
        Ok(())
    }

    /// Iterate over every (sheet name, sheet) pair.
    pub fn iter_sheets(&self) -> impl Iterator<Item = (&str, &Sheet)> {
        self.sheets.iter().map(|(k, v)| (k.as_str(), v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sheet_new_starts_empty() {
        let s = Sheet::new("A");
        assert_eq!(s.name, "A");
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());
    }

    #[test]
    fn sheet_set_and_get_round_trip() {
        let mut s = Sheet::new("A");
        s.set(0, 0, Cell::Number(1.5));
        assert!(matches!(s.get(0, 0), Cell::Number(n) if (*n - 1.5).abs() < 1e-12));
        assert!(matches!(s.get(1, 1), Cell::Empty));
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn sheet_setting_empty_removes_entry() {
        let mut s = Sheet::new("A");
        s.set(0, 0, Cell::Number(1.0));
        assert_eq!(s.len(), 1);
        s.set(0, 0, Cell::Empty);
        assert_eq!(s.len(), 0);
        assert!(matches!(s.get(0, 0), Cell::Empty));
    }

    #[test]
    fn sheet_iter_yields_only_non_empty() {
        let mut s = Sheet::new("A");
        s.set(0, 0, Cell::Number(1.0));
        s.set(1, 0, Cell::Text("hi".into()));
        let pairs: Vec<_> = s.iter().collect();
        assert_eq!(pairs.len(), 2);
    }

    #[test]
    fn spreadsheet_add_sheet_idempotent() {
        let mut ss = Spreadsheet::new();
        assert!(ss.add_sheet("A"));
        assert!(!ss.add_sheet("A"));
        assert_eq!(ss.sheets.len(), 1);
    }

    #[test]
    fn spreadsheet_set_cell_creates_value() {
        let mut ss = Spreadsheet::new();
        ss.add_sheet("Default");
        let r = CellRef::parse("Default.A1").unwrap();
        ss.set_cell(&r, Cell::Number(42.0)).unwrap();
        assert!(matches!(ss.cell(&r), Cell::Number(n) if (*n - 42.0).abs() < 1e-12));
    }

    #[test]
    fn spreadsheet_cell_missing_sheet_returns_empty() {
        let ss = Spreadsheet::new();
        let r = CellRef::parse("Missing.A1").unwrap();
        assert!(matches!(ss.cell(&r), Cell::Empty));
    }

    #[test]
    fn spreadsheet_set_cell_unknown_sheet_errors() {
        let mut ss = Spreadsheet::new();
        let r = CellRef::parse("Missing.A1").unwrap();
        let err = ss.set_cell(&r, Cell::Number(1.0)).unwrap_err();
        assert_eq!(err.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn spreadsheet_remove_sheet_returns_true_on_hit() {
        let mut ss = Spreadsheet::new();
        ss.add_sheet("A");
        assert!(ss.remove_sheet("A"));
        assert!(!ss.remove_sheet("A"));
    }
}
