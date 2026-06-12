//! FASTA reader and writer.
//!
//! A FASTA record is a `>` header line followed by one or more
//! sequence lines (which the reader concatenates, ignoring
//! whitespace). The reader is alphabet-agnostic — the caller supplies
//! the [`SeqKind`] since FASTA carries no type marker. A streaming
//! iterator ([`FastaReader`]) is provided for large multi-record
//! files.

use crate::alphabet::SeqKind;
use crate::error::{BioseqError, Result};
use crate::record::SeqRecord;
use crate::seq::Seq;

/// Splits a FASTA header (`>id description...`) into `(id, description)`.
/// The `>` must already be stripped.
fn split_header(header: &str) -> (String, String) {
    let h = header.trim();
    match h.split_once(char::is_whitespace) {
        Some((id, rest)) => (id.to_string(), rest.trim().to_string()),
        None => (h.to_string(), String::new()),
    }
}

/// Parses an entire FASTA string into a list of [`SeqRecord`]s.
///
/// Every record's sequence is validated against `kind`. Returns
/// [`BioseqError::Parse`] if the text contains sequence data before
/// the first `>` header, or [`BioseqError::Alphabet`] on an illegal
/// residue.
pub fn parse(text: &str, kind: SeqKind) -> Result<Vec<SeqRecord>> {
    let mut records = Vec::new();
    let mut cur_header: Option<String> = None;
    let mut cur_seq = String::new();

    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim_end();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            // Flush the previous record.
            if let Some(h) = cur_header.take() {
                records.push(build_record(&h, &cur_seq, kind)?);
                cur_seq.clear();
            }
            cur_header = Some(rest.to_string());
        } else if line.starts_with(';') {
            // Old-style FASTA comment line — ignore.
            continue;
        } else {
            if cur_header.is_none() {
                return Err(BioseqError::parse(
                    "fasta",
                    format!(
                        "sequence data before first '>' header (line {})",
                        lineno + 1
                    ),
                ));
            }
            cur_seq.push_str(line.trim());
        }
    }
    if let Some(h) = cur_header.take() {
        records.push(build_record(&h, &cur_seq, kind)?);
    }
    Ok(records)
}

/// Parses a FASTA string expected to hold exactly one record.
///
/// Returns [`BioseqError::Parse`] if the file holds zero or more than
/// one record.
pub fn parse_single(text: &str, kind: SeqKind) -> Result<SeqRecord> {
    let mut recs = parse(text, kind)?;
    match recs.len() {
        1 => Ok(recs.pop().expect("len checked")),
        n => Err(BioseqError::parse(
            "fasta",
            format!("expected exactly one record, found {n}"),
        )),
    }
}

fn build_record(header: &str, seq_str: &str, kind: SeqKind) -> Result<SeqRecord> {
    let (id, description) = split_header(header);
    let seq = Seq::new(kind, seq_str)?;
    Ok(SeqRecord {
        name: id.clone(),
        id,
        description,
        seq,
        features: Vec::new(),
        annotations: Default::default(),
        references: Vec::new(),
    })
}

/// A streaming, allocation-light FASTA record iterator.
///
/// Splits the input into per-record chunks lazily so a multi-gigabyte
/// FASTA does not have to be fully materialized as a `Vec`. Each
/// [`Iterator::next`] yields one `Result<SeqRecord>`.
pub struct FastaReader<'a> {
    lines: std::str::Lines<'a>,
    kind: SeqKind,
    pending_header: Option<String>,
    done: bool,
}

impl<'a> FastaReader<'a> {
    /// Builds a streaming reader over `text`.
    pub fn new(text: &'a str, kind: SeqKind) -> Self {
        FastaReader {
            lines: text.lines(),
            kind,
            pending_header: None,
            done: false,
        }
    }
}

impl Iterator for FastaReader<'_> {
    type Item = Result<SeqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        // Find a header if we don't already hold one.
        let mut header = self.pending_header.take();
        let mut seq = String::new();
        loop {
            let line = match self.lines.next() {
                Some(l) => l.trim_end(),
                None => {
                    self.done = true;
                    break;
                }
            };
            if line.is_empty() || line.starts_with(';') {
                continue;
            }
            if let Some(rest) = line.strip_prefix('>') {
                if header.is_none() {
                    header = Some(rest.to_string());
                } else {
                    // Start of the next record — stash and stop.
                    self.pending_header = Some(rest.to_string());
                    break;
                }
            } else {
                if header.is_none() {
                    self.done = true;
                    return Some(Err(BioseqError::parse(
                        "fasta",
                        "sequence data before first '>' header",
                    )));
                }
                seq.push_str(line.trim());
            }
        }
        let header = header?;
        Some(build_record(&header, &seq, self.kind))
    }
}

/// Options for [`write()`] / [`write_record`].
#[derive(Copy, Clone, Debug)]
pub struct FastaWriteOptions {
    /// Residues per sequence line. `0` means no wrapping (one long
    /// line). Default `70`.
    pub line_width: usize,
}

impl Default for FastaWriteOptions {
    fn default() -> Self {
        FastaWriteOptions { line_width: 70 }
    }
}

/// Serializes one record to a FASTA string.
pub fn write_record(rec: &SeqRecord, opts: FastaWriteOptions) -> String {
    let mut out = String::new();
    out.push('>');
    out.push_str(&rec.id);
    if !rec.description.is_empty() {
        out.push(' ');
        out.push_str(&rec.description);
    }
    out.push('\n');
    let bytes = rec.seq.as_bytes();
    if opts.line_width == 0 {
        out.push_str(rec.seq.as_str());
        out.push('\n');
    } else {
        for chunk in bytes.chunks(opts.line_width) {
            out.push_str(std::str::from_utf8(chunk).expect("residues are ASCII"));
            out.push('\n');
        }
    }
    out
}

/// Serializes a list of records to a multi-record FASTA string.
pub fn write(records: &[SeqRecord], opts: FastaWriteOptions) -> String {
    let mut out = String::new();
    for rec in records {
        out.push_str(&write_record(rec, opts));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const MULTI: &str = ">seq1 first sequence\nACGTACGT\nACGT\n>seq2 second\nTTTT\n";

    #[test]
    fn parse_multi_record() {
        let recs = parse(MULTI, SeqKind::Dna).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "seq1");
        assert_eq!(recs[0].description, "first sequence");
        assert_eq!(recs[0].seq.as_str(), "ACGTACGTACGT");
        assert_eq!(recs[1].id, "seq2");
        assert_eq!(recs[1].seq.as_str(), "TTTT");
    }

    #[test]
    fn parse_single_record() {
        let r = parse_single(">only\nACGT\n", SeqKind::Dna).unwrap();
        assert_eq!(r.id, "only");
        assert!(parse_single(MULTI, SeqKind::Dna).is_err());
        assert!(parse_single("", SeqKind::Dna).is_err());
    }

    #[test]
    fn header_with_no_description() {
        let r = parse_single(">bare\nACGT\n", SeqKind::Dna).unwrap();
        assert_eq!(r.id, "bare");
        assert_eq!(r.description, "");
    }

    #[test]
    fn data_before_header_is_error() {
        assert!(parse("ACGT\n>seq\nACGT\n", SeqKind::Dna).is_err());
    }

    #[test]
    fn invalid_residue_is_error() {
        assert!(parse(">s\nACGZ\n", SeqKind::Dna).is_err());
    }

    #[test]
    fn streaming_reader_matches_parse() {
        let streamed: Vec<SeqRecord> = FastaReader::new(MULTI, SeqKind::Dna)
            .map(|r| r.unwrap())
            .collect();
        let batched = parse(MULTI, SeqKind::Dna).unwrap();
        assert_eq!(streamed, batched);
    }

    #[test]
    fn write_wraps_lines() {
        let seq = Seq::new(SeqKind::Dna, "ACGTACGTACGT").unwrap();
        let rec = SeqRecord::new("s", seq).with_description("d");
        let out = write_record(&rec, FastaWriteOptions { line_width: 4 });
        assert_eq!(out, ">s d\nACGT\nACGT\nACGT\n");
    }

    #[test]
    fn write_no_wrap() {
        let seq = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let rec = SeqRecord::new("s", seq);
        let out = write_record(&rec, FastaWriteOptions { line_width: 0 });
        assert_eq!(out, ">s\nACGTACGT\n");
    }

    #[test]
    fn write_then_parse_roundtrip() {
        let recs = parse(MULTI, SeqKind::Dna).unwrap();
        let text = write(&recs, FastaWriteOptions::default());
        let reparsed = parse(&text, SeqKind::Dna).unwrap();
        assert_eq!(recs, reparsed);
    }
}
