//! Spreadsheet workbench error taxonomy.
//!
//! Same stable `code()` / `category()` pattern as `valenx-sketch`,
//! `valenx-cad`, `valenx-surface`, etc. — every variant maps to a
//! kebab-cased identifier that LLM / scripting layers can branch on
//! without parsing the human-readable message.

use thiserror::Error;

/// Errors raised by spreadsheet cell parsing, formula parsing, or
/// evaluation.
#[derive(Debug, Error)]
pub enum SpreadsheetError {
    /// A `"Sheet.A1"` style cell reference could not be parsed — bad
    /// alphabetic column, missing row number, malformed sheet prefix,
    /// etc.
    #[error("bad cell reference `{0}`")]
    BadCellRef(String),

    /// Formula parser bailed out at `position` characters into `input`
    /// with `reason` explaining what was expected.
    #[error("parse error at position {position} in `{input}`: {reason}")]
    ParseError {
        /// Original formula source.
        input: String,
        /// Zero-based character offset where parsing failed.
        position: usize,
        /// Human-readable explanation of what the parser expected.
        reason: String,
    },

    /// Evaluator hit something it could not compute (text cell where a
    /// number was expected, unknown function, type mismatch, etc.).
    #[error("evaluation error: {0}")]
    EvaluationError(String),

    /// Cell A1 references B1 which references A1 — recursion detected
    /// during evaluation.
    #[error("circular reference involving `{0}`")]
    CircularReference(String),

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON serialise / parse error.
    #[error("ron: {0}")]
    Ron(String),
}

impl SpreadsheetError {
    /// Stable kebab-cased identifier; never changes across versions.
    ///
    /// Mirrors `SketchError::code()` / `SurfaceError::code()` so the
    /// same scripting layer can branch on errors from any workbench.
    pub fn code(&self) -> &'static str {
        match self {
            SpreadsheetError::BadCellRef(_) => "spreadsheet.bad_cell_ref",
            SpreadsheetError::ParseError { .. } => "spreadsheet.parse_error",
            SpreadsheetError::EvaluationError(_) => "spreadsheet.evaluation_error",
            SpreadsheetError::CircularReference(_) => "spreadsheet.circular_reference",
            SpreadsheetError::Io(_) => "spreadsheet.io",
            SpreadsheetError::Ron(_) => "spreadsheet.ron",
        }
    }

    /// Coarse classification — `"input"`, `"evaluation"`, or `"io"`.
    /// Useful when the caller wants to retry on IO failures but route
    /// input/evaluation failures directly to the user.
    pub fn category(&self) -> &'static str {
        match self {
            SpreadsheetError::BadCellRef(_) | SpreadsheetError::ParseError { .. } => "input",
            SpreadsheetError::EvaluationError(_) | SpreadsheetError::CircularReference(_) => {
                "evaluation"
            }
            SpreadsheetError::Io(_) | SpreadsheetError::Ron(_) => "io",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code() {
        let cases: Vec<(SpreadsheetError, &str, &str)> = vec![
            (
                SpreadsheetError::BadCellRef("BadRef".into()),
                "spreadsheet.bad_cell_ref",
                "input",
            ),
            (
                SpreadsheetError::ParseError {
                    input: "1+".into(),
                    position: 2,
                    reason: "expected expression".into(),
                },
                "spreadsheet.parse_error",
                "input",
            ),
            (
                SpreadsheetError::EvaluationError("text in numeric context".into()),
                "spreadsheet.evaluation_error",
                "evaluation",
            ),
            (
                SpreadsheetError::CircularReference("Sheet.A1".into()),
                "spreadsheet.circular_reference",
                "evaluation",
            ),
            (
                SpreadsheetError::Ron("malformed".into()),
                "spreadsheet.ron",
                "io",
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn io_error_wraps_std() {
        let e: SpreadsheetError = std::io::Error::other("disk gone").into();
        assert_eq!(e.code(), "spreadsheet.io");
        assert_eq!(e.category(), "io");
    }
}
