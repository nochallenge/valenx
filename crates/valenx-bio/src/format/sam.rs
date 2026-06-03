//! Minimal SAM-text reader.
//!
//! Reads the 11 mandatory fields per record (qname / flag / rname /
//! pos / mapq / cigar / rnext / pnext / tlen / seq / qual) and
//! collects header lines verbatim. Optional tags (TAG:TYPE:VALUE)
//! are stored as raw strings — the canonical type does not parse
//! them. BAM (binary BGZF-encoded SAM) is intentionally out of scope;
//! convert with `samtools view -h` first.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A parsed SAM file — header lines (kept verbatim) plus
/// [`SamRecord`]s.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sam {
    /// SAM header lines (those beginning with `@`), kept verbatim in
    /// source order.
    pub header: Vec<String>,
    /// Alignment records in file order.
    pub records: Vec<SamRecord>,
}

/// One SAM alignment record. Field names mirror the SAM v1 spec.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SamRecord {
    /// Read name (`QNAME`).
    pub qname: String,
    /// Bitwise flag (`FLAG`) — see [`SamRecord::is_mapped`].
    pub flag: u32,
    /// Reference sequence name (`RNAME`).
    pub rname: String,
    /// 1-based leftmost mapping position (`POS`).
    pub pos: u64,
    /// Mapping quality (`MAPQ`).
    pub mapq: u8,
    /// CIGAR string describing the alignment.
    pub cigar: String,
    /// Reference name of mate / next read (`RNEXT`).
    pub rnext: String,
    /// 1-based position of mate / next read (`PNEXT`).
    pub pnext: u64,
    /// Observed template length (`TLEN`).
    pub tlen: i64,
    /// Read bases (`SEQ`).
    pub seq: String,
    /// Per-base ASCII quality (`QUAL`).
    pub qual: String,
    /// Optional tag fields (`TAG:TYPE:VALUE`), one entry per tag.
    pub tags: Vec<String>,
}

impl SamRecord {
    /// SAM flag bit 4 = "segment unmapped".
    pub fn is_mapped(&self) -> bool {
        (self.flag & 0x4) == 0
    }
}

/// Errors raised by [`read_str`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum SamError {
    /// A non-header line was malformed (fewer than 11 tab fields or
    /// a numeric field that did not parse).
    #[error("line {line}: {msg}")]
    Bad {
        /// 1-based line number of the offending record.
        line: usize,
        /// Short human-readable explanation.
        msg: String,
    },
    /// A record violated the SAM spec semantically (e.g. SEQ and QUAL
    /// of unequal length). Distinct from `Bad` so downstream tools can
    /// telemetry "spec violation" separately from "couldn't parse".
    #[error("malformed SAM: {0}")]
    Malformed(String),
}

/// Parse a SAM-format string into a [`Sam`] container.
///
/// Header lines (those beginning with `@`) are captured verbatim;
/// alignment records are split on tab and the first 11 fields are
/// type-checked. Trailing fields become [`SamRecord::tags`].
///
/// # Errors
///
/// Returns [`SamError::Bad`] when a record line has fewer than 11
/// fields or when a numeric field fails to parse.
pub fn read_str(s: &str) -> Result<Sam, SamError> {
    let mut sam = Sam::default();
    for (i, line) in s.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        if line.starts_with('@') {
            sam.header.push(line.to_string());
            continue;
        }
        let line_no = i + 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return Err(SamError::Bad {
                line: line_no,
                msg: format!("expected at least 11 fields, got {}", fields.len()),
            });
        }
        let parse_u32 = |x: &str| -> Result<u32, SamError> {
            x.parse().map_err(|_| SamError::Bad {
                line: line_no,
                msg: format!("not an integer: `{x}`"),
            })
        };
        let parse_u64 = |x: &str| -> Result<u64, SamError> {
            x.parse().map_err(|_| SamError::Bad {
                line: line_no,
                msg: format!("not an integer: `{x}`"),
            })
        };
        let parse_u8 = |x: &str| -> Result<u8, SamError> {
            x.parse().map_err(|_| SamError::Bad {
                line: line_no,
                msg: format!("not a small integer: `{x}`"),
            })
        };
        let parse_i64 = |x: &str| -> Result<i64, SamError> {
            x.parse().map_err(|_| SamError::Bad {
                line: line_no,
                msg: format!("not an integer: `{x}`"),
            })
        };
        let seq = fields[9];
        let qual = fields[10];
        // SAM v1 §1.4: when SEQ ≠ "*" and QUAL ≠ "*", QUAL must have
        // the same length as SEQ. A silent length drift would let
        // downstream `seq.bytes().zip(qual.bytes())` truncate to the
        // shorter of the two and quietly drop bases — refuse here so
        // the malformed file is caught at parse time.
        if seq != "*" && qual != "*" && seq.len() != qual.len() {
            return Err(SamError::Malformed(format!(
                "line {line_no}: SEQ length {} != QUAL length {} (spec violation)",
                seq.len(),
                qual.len()
            )));
        }
        sam.records.push(SamRecord {
            qname: fields[0].to_string(),
            flag: parse_u32(fields[1])?,
            rname: fields[2].to_string(),
            pos: parse_u64(fields[3])?,
            mapq: parse_u8(fields[4])?,
            cigar: fields[5].to_string(),
            rnext: fields[6].to_string(),
            pnext: parse_u64(fields[7])?,
            tlen: parse_i64(fields[8])?,
            seq: seq.to_string(),
            qual: qual.to_string(),
            tags: fields[11..].iter().map(|s| s.to_string()).collect(),
        });
    }
    Ok(sam)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn happy_path_parses_header_and_one_record() {
        let sam_text = "@HD\tVN:1.6\tSO:coordinate\n\
                        @SQ\tSN:chr1\tLN:1000\n\
                        rec1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n";
        let sam = read_str(sam_text).unwrap();
        assert_eq!(sam.header.len(), 2);
        assert_eq!(sam.records.len(), 1);
        let r = &sam.records[0];
        assert_eq!(r.qname, "rec1");
        assert_eq!(r.flag, 0);
        assert_eq!(r.rname, "chr1");
        assert_eq!(r.pos, 100);
        assert_eq!(r.seq, "ACGTACGT");
        assert_eq!(r.qual, "IIIIIIII");
    }

    #[test]
    fn missing_qual_is_allowed_when_star() {
        let sam_text = "@HD\tVN:1.6\n\
                        rec1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\t*\n";
        let sam = read_str(sam_text).unwrap();
        assert_eq!(sam.records.len(), 1);
        assert_eq!(sam.records[0].qual, "*");
    }

    #[test]
    fn seq_qual_length_mismatch_rejected() {
        // SAM v1 §1.4: when both SEQ and QUAL are non-"*", they must have
        // the same length. Round-3 added the validation; this test pins it.
        let header = "@HD\tVN:1.6\n";
        let record = "rec1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGT\tIIII\n"; // seq=8, qual=4
        let result = read_str(&format!("{header}{record}"));
        assert!(matches!(result, Err(SamError::Malformed(_))));
        let err = result.unwrap_err();
        if let SamError::Malformed(s) = err {
            assert!(s.contains("SEQ"), "msg: {s}");
            assert!(s.contains("QUAL"), "msg: {s}");
        } else {
            panic!("expected Malformed, got {err:?}");
        }
    }
}
