//! 4-line FASTQ reader + writer.
//!
//! The format: header `@<name> [description]`, sequence line,
//! `+` separator (with optional repeated header), quality line.
//! This implementation supports the strict 4-line variant only —
//! line-wrapped sequences are rare in modern FASTQ and not worth
//! the parser complexity.

use std::io;

use thiserror::Error;

use crate::fastq::{FastqError, FastqRecord};

/// Errors raised by [`read_str`].
#[derive(Debug, Error)]
pub enum ReadError {
    /// An underlying IO error (only surfaced by callers that wrap
    /// [`read_str`] with a streaming source).
    #[error("io: {0}")]
    Io(#[from] io::Error),
    /// FASTQ record shape was wrong (missing `@` header, missing `+`
    /// separator, fewer than 4 lines, etc.).
    #[error("record at line {line}: {msg}")]
    BadRecord {
        /// 1-based line number of the offending header.
        line: usize,
        /// Human-readable explanation.
        msg: String,
    },
    /// Record fields were well-shaped but failed [`FastqRecord::new`]
    /// validation (e.g. quality length mismatch).
    #[error("record at line {line}: {err}")]
    Record {
        /// 1-based line number of the offending header.
        line: usize,
        /// The underlying [`FastqError`].
        err: FastqError,
    },
}

/// Parse a 4-line FASTQ-shaped string into `FastqRecord`s. Empty
/// lines between records are skipped; quality / sequence length
/// mismatches and missing `+` separators are reported as
/// `BadRecord` / `Record` errors with the offending line number.
pub fn read_str(s: &str) -> Result<Vec<FastqRecord>, ReadError> {
    let mut out = Vec::new();
    let mut lines = s.lines().enumerate();
    while let Some((i, header)) = lines.next() {
        if header.is_empty() {
            continue;
        }
        let header = header
            .strip_prefix('@')
            .ok_or_else(|| ReadError::BadRecord {
                line: i + 1,
                msg: format!("expected '@' at start, got `{header}`"),
            })?;
        let (name, description) = match header.split_once(char::is_whitespace) {
            Some((n, d)) => (n.to_string(), Some(d.to_string())),
            None => (header.to_string(), None),
        };
        let (_, seq) = lines.next().ok_or_else(|| ReadError::BadRecord {
            line: i + 1,
            msg: "missing sequence line".into(),
        })?;
        let (plus_idx, plus) = lines.next().ok_or_else(|| ReadError::BadRecord {
            line: i + 1,
            msg: "missing '+' separator".into(),
        })?;
        if !plus.starts_with('+') {
            return Err(ReadError::BadRecord {
                line: plus_idx + 1,
                msg: format!("expected '+' separator, got `{plus}`"),
            });
        }
        let (_, qual) = lines.next().ok_or_else(|| ReadError::BadRecord {
            line: i + 1,
            msg: "missing quality line".into(),
        })?;
        let rec = FastqRecord::new(
            name,
            description,
            seq.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .map_err(|err| ReadError::Record { line: i + 1, err })?;
        out.push(rec);
    }
    Ok(out)
}

/// Render a slice of `FastqRecord`s back to FASTQ text in the
/// strict 4-line shape (header / sequence / `+` / quality).
pub fn write_string(records: &[FastqRecord]) -> Result<String, std::fmt::Error> {
    use std::fmt::Write;
    let mut s = String::new();
    for rec in records {
        match &rec.description {
            Some(d) => writeln!(s, "@{} {}", rec.name, d)?,
            None => writeln!(s, "@{}", rec.name)?,
        }
        s.push_str(std::str::from_utf8(&rec.sequence).unwrap_or(""));
        s.push('\n');
        s.push_str("+\n");
        s.push_str(std::str::from_utf8(&rec.quality).unwrap_or(""));
        s.push('\n');
    }
    Ok(s)
}
