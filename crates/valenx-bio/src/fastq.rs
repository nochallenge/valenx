//! FASTQ records: name + (optional) description + sequence bytes
//! + per-base Phred+33 quality bytes.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A single FASTQ read.
///
/// Quality bytes are stored verbatim in their Phred+33 encoding —
/// callers that want raw integer scores must subtract 33 themselves.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FastqRecord {
    /// Read name (the text after `@` on the title line).
    pub name: String,
    /// Optional free-text description that followed the name.
    pub description: Option<String>,
    /// Read residues, verbatim from the FASTQ file.
    pub sequence: Vec<u8>,
    /// Per-base Phred+33 quality bytes — same length as `sequence`.
    pub quality: Vec<u8>,
}

/// Construction errors from [`FastqRecord::new`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum FastqError {
    /// `sequence` and `quality` have different lengths — every base
    /// must have exactly one quality byte.
    #[error("sequence length {seq} != quality length {qual}")]
    LengthMismatch {
        /// Number of sequence bytes.
        seq: usize,
        /// Number of quality bytes.
        qual: usize,
    },
    /// The supplied read name is the empty string.
    #[error("name must not be empty")]
    EmptyName,
}

impl FastqRecord {
    /// Build a `FastqRecord`, validating that `sequence.len() ==
    /// quality.len()` and the name is non-empty.
    pub fn new(
        name: String,
        description: Option<String>,
        sequence: Vec<u8>,
        quality: Vec<u8>,
    ) -> Result<Self, FastqError> {
        if name.is_empty() {
            return Err(FastqError::EmptyName);
        }
        if sequence.len() != quality.len() {
            return Err(FastqError::LengthMismatch {
                seq: sequence.len(),
                qual: quality.len(),
            });
        }
        Ok(Self {
            name,
            description,
            sequence,
            quality,
        })
    }

    /// Read length in bases.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// `true` if the read has zero bases.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Minimum Phred score across the read, with the +33 offset
    /// already removed. `None` when the read is empty.
    pub fn min_quality(&self) -> Option<u8> {
        self.quality.iter().min().map(|q| q.saturating_sub(33))
    }

    /// Maximum Phred score across the read, with the +33 offset
    /// already removed. `None` when the read is empty.
    pub fn max_quality(&self) -> Option<u8> {
        self.quality.iter().max().map(|q| q.saturating_sub(33))
    }
}
