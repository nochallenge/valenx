//! FASTA format reader + writer.
//!
//! The format spec is informal but well-converged:
//! - Records start with `>` followed by an arbitrary identifier line
//!   (everything after `>` up to newline becomes the sequence name).
//! - Lines until the next `>` (or EOF) are the sequence body.
//! - Comments start with `;` (rare; we tolerate but don't preserve).
//! - Whitespace within a body line is stripped (some sources mix in
//!   numbers / spaces for readability).
//!
//! We don't preserve original line widths on write — sequences emit
//! with `LINE_WIDTH = 60` per the historical convention.

use thiserror::Error;

use crate::sequence::{Sequence, SequenceError};
use crate::Alphabet;

const LINE_WIDTH: usize = 60;

/// Errors raised by [`read`].
#[derive(Debug, Error)]
pub enum FastaError {
    /// A record body failed [`Sequence`] alphabet validation.
    #[error(transparent)]
    Sequence(#[from] SequenceError),
}

/// Parse a FASTA-shaped string into `Sequence`s. Every record's body
/// is validated against `alphabet`.
pub fn read(text: &str, alphabet: Alphabet) -> Result<Vec<Sequence>, FastaError> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut body = String::new();

    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        if let Some(stripped) = line.strip_prefix('>') {
            if let Some(prev_name) = name.take() {
                out.push(Sequence::new(prev_name, alphabet, &body)?);
                body.clear();
            }
            name = Some(stripped.trim().to_string());
        } else {
            for c in line.chars() {
                if !c.is_whitespace() && !c.is_ascii_digit() {
                    body.push(c);
                }
            }
        }
    }
    if let Some(final_name) = name {
        out.push(Sequence::new(final_name, alphabet, &body)?);
    }
    Ok(out)
}

/// Render a slice of `Sequence`s back to FASTA text. Body lines wrap
/// at 60 characters per the historical FASTA convention.
pub fn write(seqs: &[Sequence]) -> String {
    let mut out = String::new();
    for s in seqs {
        out.push('>');
        out.push_str(&s.name);
        out.push('\n');
        let body = s.as_str();
        let mut i = 0;
        while i < body.len() {
            let end = (i + LINE_WIDTH).min(body.len());
            out.push_str(&body[i..end]);
            out.push('\n');
            i = end;
        }
    }
    out
}
