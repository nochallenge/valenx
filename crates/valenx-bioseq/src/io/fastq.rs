//! FASTQ reader and writer.
//!
//! A FASTQ record is four lines: `@id description`, the sequence, a
//! `+` separator line, and the quality string. The reader parses the
//! quality string into per-base Phred scores via
//! [`crate::io::quality`]; a [`FastqRecord`] carries the sequence as a
//! [`SeqRecord`] plus the decoded `Vec<u8>` of Phred scores.

use crate::alphabet::SeqKind;
use crate::error::{BioseqError, Result};
use crate::io::quality::{self, QualityEncoding};
use crate::record::SeqRecord;
use crate::seq::Seq;

/// A FASTQ record — a sequence plus its per-base Phred quality scores.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqRecord {
    /// The sequence with id / description.
    pub record: SeqRecord,
    /// Per-base Phred quality scores; same length as the sequence.
    pub quality: Vec<u8>,
}

impl FastqRecord {
    /// Mean Phred quality across all bases.
    pub fn mean_quality(&self) -> f64 {
        quality::mean_quality(&self.quality)
    }

    /// Number of bases (and quality values).
    pub fn len(&self) -> usize {
        self.record.seq.len()
    }

    /// `true` if the record holds no bases.
    pub fn is_empty(&self) -> bool {
        self.record.seq.is_empty()
    }
}

fn split_header(header: &str) -> (String, String) {
    match header.split_once(char::is_whitespace) {
        Some((id, rest)) => (id.to_string(), rest.trim().to_string()),
        None => (header.to_string(), String::new()),
    }
}

/// Parses an entire FASTQ string into a list of [`FastqRecord`]s.
///
/// `kind` declares the sequence alphabet (FASTQ carries no type
/// marker); `encoding` selects the Phred / Solexa offset. Returns
/// [`BioseqError::Parse`] on a malformed record (wrong line count,
/// missing `@` / `+`, or a sequence/quality length mismatch).
pub fn parse(
    text: &str,
    kind: SeqKind,
    encoding: QualityEncoding,
) -> Result<Vec<FastqRecord>> {
    let lines: Vec<&str> = text.lines().map(|l| l.trim_end()).collect();
    // Drop a trailing empty line if present.
    let lines: Vec<&str> = lines
        .into_iter()
        .filter(|l| !l.is_empty())
        .collect();
    if lines.len() % 4 != 0 {
        return Err(BioseqError::parse(
            "fastq",
            format!(
                "line count {} is not a multiple of 4 (records are 4-line)",
                lines.len()
            ),
        ));
    }
    let mut records = Vec::with_capacity(lines.len() / 4);
    for (rec_idx, chunk) in lines.chunks(4).enumerate() {
        records.push(parse_record(chunk, rec_idx, kind, encoding)?);
    }
    Ok(records)
}

fn parse_record(
    chunk: &[&str],
    rec_idx: usize,
    kind: SeqKind,
    encoding: QualityEncoding,
) -> Result<FastqRecord> {
    let header = chunk[0].strip_prefix('@').ok_or_else(|| {
        BioseqError::parse(
            "fastq",
            format!("record {rec_idx}: header line must start with '@'"),
        )
    })?;
    let seq_line = chunk[1];
    let sep = chunk[2];
    if !sep.starts_with('+') {
        return Err(BioseqError::parse(
            "fastq",
            format!("record {rec_idx}: separator line must start with '+'"),
        ));
    }
    let qual_line = chunk[3];
    if seq_line.len() != qual_line.len() {
        return Err(BioseqError::parse(
            "fastq",
            format!(
                "record {rec_idx}: sequence length {} != quality length {}",
                seq_line.len(),
                qual_line.len()
            ),
        ));
    }
    let (id, description) = split_header(header);
    let seq = Seq::new(kind, seq_line)?;
    let quality = quality::decode(qual_line.as_bytes(), encoding)?;
    let record = SeqRecord {
        name: id.clone(),
        id,
        description,
        seq,
        features: Vec::new(),
        annotations: Default::default(),
        references: Vec::new(),
    };
    Ok(FastqRecord { record, quality })
}

/// A streaming FASTQ record iterator for large files.
pub struct FastqReader<'a> {
    lines: std::str::Lines<'a>,
    kind: SeqKind,
    encoding: QualityEncoding,
    index: usize,
}

impl<'a> FastqReader<'a> {
    /// Builds a streaming reader.
    pub fn new(text: &'a str, kind: SeqKind, encoding: QualityEncoding) -> Self {
        FastqReader {
            lines: text.lines(),
            kind,
            encoding,
            index: 0,
        }
    }

    /// Pulls the next non-empty line.
    fn next_nonempty(&mut self) -> Option<&'a str> {
        for l in self.lines.by_ref() {
            let t = l.trim_end();
            if !t.is_empty() {
                return Some(t);
            }
        }
        None
    }
}

impl Iterator for FastqReader<'_> {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let l0 = self.next_nonempty()?;
        // Once a record starts, all four lines must be present.
        let l1 = self.next_nonempty();
        let l2 = self.next_nonempty();
        let l3 = self.next_nonempty();
        let (l1, l2, l3) = match (l1, l2, l3) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => {
                return Some(Err(BioseqError::parse(
                    "fastq",
                    "truncated record (fewer than 4 lines)",
                )))
            }
        };
        let idx = self.index;
        self.index += 1;
        Some(parse_record(&[l0, l1, l2, l3], idx, self.kind, self.encoding))
    }
}

/// Serializes one FASTQ record to a string.
pub fn write_record(rec: &FastqRecord, encoding: QualityEncoding) -> String {
    let mut out = String::new();
    out.push('@');
    out.push_str(&rec.record.id);
    if !rec.record.description.is_empty() {
        out.push(' ');
        out.push_str(&rec.record.description);
    }
    out.push('\n');
    out.push_str(rec.record.seq.as_str());
    out.push('\n');
    out.push_str("+\n");
    let qual = quality::encode(&rec.quality, encoding);
    out.push_str(std::str::from_utf8(&qual).expect("quality chars are ASCII"));
    out.push('\n');
    out
}

/// Serializes a list of FASTQ records to a string.
pub fn write(records: &[FastqRecord], encoding: QualityEncoding) -> String {
    let mut out = String::new();
    for rec in records {
        out.push_str(&write_record(rec, encoding));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str =
        "@read1 desc\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTT\n+\n!!!!\n";

    #[test]
    fn parse_two_records() {
        let recs = parse(SAMPLE, SeqKind::Dna, QualityEncoding::Phred33).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].record.id, "read1");
        assert_eq!(recs[0].record.description, "desc");
        assert_eq!(recs[0].record.seq.as_str(), "ACGTACGT");
        // 'I' = 73 -> Q40.
        assert_eq!(recs[0].quality, vec![40; 8]);
        // '!' = 33 -> Q0.
        assert_eq!(recs[1].quality, vec![0; 4]);
    }

    #[test]
    fn mean_quality_works() {
        let recs = parse(SAMPLE, SeqKind::Dna, QualityEncoding::Phred33).unwrap();
        assert_eq!(recs[0].mean_quality(), 40.0);
        assert_eq!(recs[1].mean_quality(), 0.0);
    }

    #[test]
    fn length_mismatch_is_error() {
        let bad = "@r\nACGT\n+\nIII\n"; // 4 bases, 3 quals
        assert!(parse(bad, SeqKind::Dna, QualityEncoding::Phred33).is_err());
    }

    #[test]
    fn wrong_line_count_is_error() {
        let bad = "@r\nACGT\n+\n"; // 3 lines
        assert!(parse(bad, SeqKind::Dna, QualityEncoding::Phred33).is_err());
    }

    #[test]
    fn missing_at_marker_is_error() {
        let bad = "r\nACGT\n+\nIIII\n";
        assert!(parse(bad, SeqKind::Dna, QualityEncoding::Phred33).is_err());
    }

    #[test]
    fn missing_plus_marker_is_error() {
        let bad = "@r\nACGT\nx\nIIII\n";
        assert!(parse(bad, SeqKind::Dna, QualityEncoding::Phred33).is_err());
    }

    #[test]
    fn write_then_parse_roundtrip() {
        let recs = parse(SAMPLE, SeqKind::Dna, QualityEncoding::Phred33).unwrap();
        let text = write(&recs, QualityEncoding::Phred33);
        let reparsed = parse(&text, SeqKind::Dna, QualityEncoding::Phred33).unwrap();
        assert_eq!(recs, reparsed);
    }

    #[test]
    fn streaming_reader_matches_batch() {
        let streamed: Vec<FastqRecord> =
            FastqReader::new(SAMPLE, SeqKind::Dna, QualityEncoding::Phred33)
                .map(|r| r.unwrap())
                .collect();
        let batched = parse(SAMPLE, SeqKind::Dna, QualityEncoding::Phred33).unwrap();
        assert_eq!(streamed, batched);
    }
}
