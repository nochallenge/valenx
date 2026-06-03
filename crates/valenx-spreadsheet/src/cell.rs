//! Cell references + cell payloads.
//!
//! A [`CellRef`] is the canonical address of a cell in a workbook:
//! `Sheet.A1` where `Sheet` is the owning sheet name and `A1` is the
//! A1-style column+row coordinate.
//!
//! A [`Cell`] is the payload stored at a [`CellRef`]: an empty marker,
//! a literal number, free-form text, or a formula source string that
//! the [`crate::evaluator`] resolves to a number on demand.
//!
//! Both types round-trip cleanly through serde / RON so the
//! [`crate::persist::SpreadsheetFile`] envelope can save and load
//! workbooks without losing precision.

use serde::{Deserialize, Serialize};

use crate::error::SpreadsheetError;

/// A cell address: sheet name + zero-based `(row, col)` coordinate.
///
/// Display form is `"Sheet.A1"` — `"Sheet"` is the owning sheet's
/// name, `"A"` is column 0 (A1-style), `"1"` is row 0 (1-indexed in
/// display, 0-indexed internally to keep the maths boring).
#[derive(Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct CellRef {
    /// Name of the sheet this cell lives on (case-sensitive).
    pub sheet_name: String,
    /// Zero-based row index. Display form adds 1 (row 0 prints as
    /// `"1"`).
    pub row: u32,
    /// Zero-based column index. Display form is alphabetic A=0, B=1,
    /// ..., Z=25, AA=26, AB=27, ....
    pub col: u32,
}

impl CellRef {
    /// Construct a [`CellRef`] from a sheet name and an A1-style
    /// `"B7"` coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`SpreadsheetError::BadCellRef`] when the A1 string is
    /// empty, has no row digits, has a non-alphabetic column prefix,
    /// or overflows.
    pub fn from_a1(sheet: impl Into<String>, a1: &str) -> Result<Self, SpreadsheetError> {
        let mut chars = a1.chars();
        let mut col_str = String::new();
        let mut row_str = String::new();
        let mut saw_digit = false;
        for c in chars.by_ref() {
            if c.is_ascii_alphabetic() {
                if saw_digit {
                    return Err(SpreadsheetError::BadCellRef(format!(
                        "letters cannot follow digits in `{a1}`"
                    )));
                }
                col_str.push(c.to_ascii_uppercase());
            } else if c.is_ascii_digit() {
                saw_digit = true;
                row_str.push(c);
            } else {
                return Err(SpreadsheetError::BadCellRef(format!(
                    "unexpected character `{c}` in A1 ref `{a1}`"
                )));
            }
        }
        if col_str.is_empty() {
            return Err(SpreadsheetError::BadCellRef(format!(
                "missing column letters in `{a1}`"
            )));
        }
        if row_str.is_empty() {
            return Err(SpreadsheetError::BadCellRef(format!(
                "missing row number in `{a1}`"
            )));
        }
        let col = a1_letters_to_col(&col_str)?;
        let row: u32 = row_str.parse().map_err(|_| {
            SpreadsheetError::BadCellRef(format!("row number `{row_str}` is out of range"))
        })?;
        if row == 0 {
            return Err(SpreadsheetError::BadCellRef(format!(
                "row numbers are 1-indexed in `{a1}`"
            )));
        }
        Ok(Self {
            sheet_name: sheet.into(),
            row: row - 1,
            col,
        })
    }

    /// Parse the full `"Sheet.A1"` form.
    ///
    /// The sheet name is everything before the *last* `.` in the
    /// string; the A1 coordinate is everything after it. This lets
    /// sheet names contain `.` (e.g. `"v1.2.A1"` → sheet `"v1.2"`,
    /// cell `"A1"`) but does NOT support `.` inside A1 coords.
    ///
    /// # Errors
    ///
    /// Returns [`SpreadsheetError::BadCellRef`] when the string has no
    /// `.` separator, or when either side fails its own validation.
    pub fn parse(s: &str) -> Result<Self, SpreadsheetError> {
        let dot = s.rfind('.').ok_or_else(|| {
            SpreadsheetError::BadCellRef(format!(
                "missing `.` separator between sheet name and A1 coord in `{s}`"
            ))
        })?;
        let sheet = &s[..dot];
        let a1 = &s[dot + 1..];
        if sheet.is_empty() {
            return Err(SpreadsheetError::BadCellRef(format!(
                "empty sheet name in `{s}`"
            )));
        }
        Self::from_a1(sheet.to_string(), a1)
    }

    /// Render this ref in A1 form (no sheet prefix).
    pub fn to_a1(&self) -> String {
        format!("{}{}", col_to_a1_letters(self.col), self.row + 1)
    }
}

impl std::fmt::Display for CellRef {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}", self.sheet_name, self.to_a1())
    }
}

/// Convert a string of column letters (`"A"`, `"AA"`, ...) to the
/// 0-based column index.
///
/// Treats letters as base-26 with `A`=1 in the digit alphabet, so `"A"`
/// → 0, `"Z"` → 25, `"AA"` → 26, `"AB"` → 27, ...
fn a1_letters_to_col(letters: &str) -> Result<u32, SpreadsheetError> {
    let mut col: u32 = 0;
    for c in letters.chars() {
        if !c.is_ascii_alphabetic() {
            return Err(SpreadsheetError::BadCellRef(format!(
                "non-alphabetic character `{c}` in column letters `{letters}`"
            )));
        }
        let digit = (c.to_ascii_uppercase() as u32) - ('A' as u32) + 1;
        col = col
            .checked_mul(26)
            .and_then(|c| c.checked_add(digit))
            .ok_or_else(|| {
                SpreadsheetError::BadCellRef(format!("column index for `{letters}` overflowed"))
            })?;
    }
    Ok(col - 1)
}

/// Convert a 0-based column index to A1 letters (`0` → `"A"`, ...).
fn col_to_a1_letters(mut col: u32) -> String {
    let mut letters = Vec::new();
    loop {
        let rem = (col % 26) as u8;
        letters.push((b'A' + rem) as char);
        if col < 26 {
            break;
        }
        col = col / 26 - 1;
    }
    letters.reverse();
    letters.into_iter().collect()
}

/// Payload stored in one cell of a sheet.
///
/// Empty cells aren't materialised in the sheet's HashMap; the
/// `Empty` variant is what `Sheet::get` returns for missing
/// coordinates so callers can pattern-match uniformly.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Cell {
    /// No value entered yet.
    #[default]
    Empty,
    /// Literal numeric value — evaluator returns it as-is.
    Number(f64),
    /// Free-form text — evaluator errors if asked to use it as a
    /// number.
    Text(String),
    /// Formula source string (e.g. `"=A1+B2"` or `"A1+B2"` — the
    /// leading `=` is stripped by the parser if present). Stored as
    /// text so the user can edit it; resolved to a number by the
    /// evaluator on demand.
    Formula(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn column_letters_round_trip() {
        for col in [0u32, 1, 25, 26, 27, 51, 52, 701, 702, 1000] {
            let letters = col_to_a1_letters(col);
            let back = a1_letters_to_col(&letters).unwrap();
            assert_eq!(back, col, "round-trip failed for col {col} ({letters})");
        }
    }

    #[test]
    fn parse_default_a1_basic() {
        let r = CellRef::parse("Default.A1").unwrap();
        assert_eq!(r.sheet_name, "Default");
        assert_eq!(r.row, 0);
        assert_eq!(r.col, 0);
        assert_eq!(r.to_a1(), "A1");
        assert_eq!(r.to_string(), "Default.A1");
    }

    #[test]
    fn parse_my_sheet_aa15() {
        let r = CellRef::parse("MySheet.AA15").unwrap();
        assert_eq!(r.sheet_name, "MySheet");
        assert_eq!(r.col, 26);
        assert_eq!(r.row, 14);
        assert_eq!(r.to_a1(), "AA15");
        assert_eq!(r.to_string(), "MySheet.AA15");
    }

    #[test]
    fn parse_z9_boundary() {
        let r = CellRef::parse("Sheet.Z9").unwrap();
        assert_eq!(r.col, 25);
        assert_eq!(r.row, 8);
    }

    #[test]
    fn from_a1_lowercase_normalises() {
        let r = CellRef::from_a1("S", "ab12").unwrap();
        assert_eq!(r.col, 27);
        assert_eq!(r.row, 11);
    }

    #[test]
    fn parse_missing_dot_errors() {
        let e = CellRef::parse("BadRef").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn parse_empty_sheet_errors() {
        let e = CellRef::parse(".A1").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn parse_missing_row_errors() {
        let e = CellRef::parse("Sheet.AA").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn parse_missing_col_errors() {
        let e = CellRef::parse("Sheet.7").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn parse_letters_after_digits_errors() {
        let e = CellRef::parse("Sheet.A1B").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn parse_zero_row_errors() {
        let e = CellRef::parse("Sheet.A0").unwrap_err();
        assert_eq!(e.code(), "spreadsheet.bad_cell_ref");
    }

    #[test]
    fn parse_dotted_sheet_name() {
        // Sheet name `"v1.2"` with A1 coord `"B3"`.
        let r = CellRef::parse("v1.2.B3").unwrap();
        assert_eq!(r.sheet_name, "v1.2");
        assert_eq!(r.col, 1);
        assert_eq!(r.row, 2);
    }

    #[test]
    fn cell_default_is_empty() {
        let c: Cell = Cell::default();
        assert!(matches!(c, Cell::Empty));
    }
}
