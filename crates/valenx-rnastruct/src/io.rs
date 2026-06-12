//! Connectivity-table (ct) and bpseq structure file I/O.
//!
//! Two long-standing plain-text formats pair an RNA sequence with a
//! secondary structure:
//!
//! - **ct** (connectivity table) — the mfold / RNAstructure native
//!   format. A header line `<n>  <title>` followed by one line per
//!   base: `index  base  prev  next  pair  index2`. `pair` is the
//!   1-based partner index or `0` for unpaired.
//! - **bpseq** — the Gutell-lab / Comparative RNA Web format. Optional
//!   `#`-comment lines, then one line per base: `index  base  pair`.
//!
//! Both are 1-based on disk; this module converts to/from the 0-based
//! [`Structure`] partner array. Dot-bracket I/O lives on
//! [`Structure`] itself ([`Structure::from_dot_bracket`] /
//! [`Structure::to_dot_bracket`]).

use crate::error::{Result, RnaStructError};
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// A parsed structure file: a sequence plus its secondary structure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StructureRecord {
    /// An optional title / id (the ct header text, or empty).
    pub title: String,
    /// The RNA sequence.
    pub seq: RnaSeq,
    /// The secondary structure (same length as `seq`).
    pub structure: Structure,
}

impl StructureRecord {
    /// Builds a record, checking the sequence and structure lengths
    /// agree.
    ///
    /// # Errors
    /// [`RnaStructError::Structure`] on a length mismatch.
    pub fn new(title: impl Into<String>, seq: RnaSeq, structure: Structure) -> Result<Self> {
        if seq.len() != structure.len() {
            return Err(RnaStructError::structure(format!(
                "sequence length {} != structure length {}",
                seq.len(),
                structure.len()
            )));
        }
        Ok(StructureRecord {
            title: title.into(),
            seq,
            structure,
        })
    }
}

// ---------------------------------------------------------------------
// ct format
// ---------------------------------------------------------------------

/// Parses a ct-format string into a [`StructureRecord`].
///
/// # Errors
/// [`RnaStructError::Parse`] on a malformed header / data line, a
/// non-RNA base, or an out-of-range / asymmetric partner index.
pub fn read_ct(text: &str) -> Result<StructureRecord> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines
        .next()
        .ok_or_else(|| RnaStructError::parse("ct", "empty file"))?;
    let mut hfields = header.split_whitespace();
    let n: usize = hfields
        .next()
        .ok_or_else(|| RnaStructError::parse("ct", "missing length in header"))?
        .parse()
        .map_err(|_| RnaStructError::parse("ct", "header length is not a number"))?;
    let title: String = hfields.collect::<Vec<_>>().join(" ");

    // The header length `n` is attacker-controlled. Both allocations
    // below (`Vec::with_capacity(n)` and `vec![None; n]`) happen BEFORE
    // any data line is read, so an unbounded `n` (e.g. `999999999999`
    // in a 12-byte file) attempts a multi-GB allocation → OOM. A ct
    // record stores one base per data line and each data line is ≥1
    // byte, so a structure can never have more bases than there are
    // bytes in the file — `n > text.len()` is a sound, tight upper
    // bound that rejects the bomb before allocating.
    if n > text.len() {
        return Err(RnaStructError::parse(
            "ct",
            format!("declared length {n} exceeds file size {} bytes", text.len()),
        ));
    }

    let mut bases: Vec<u8> = Vec::with_capacity(n);
    let mut partner: Vec<Option<usize>> = vec![None; n];

    for (row, line) in lines.enumerate() {
        if row >= n {
            break; // ignore trailing junk
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 5 {
            return Err(RnaStructError::parse(
                "ct",
                format!("data line {} has fewer than 5 fields", row + 1),
            ));
        }
        let idx: usize = f[0]
            .parse()
            .map_err(|_| RnaStructError::parse("ct", "bad base index"))?;
        if idx != row + 1 {
            return Err(RnaStructError::parse(
                "ct",
                format!("base index {idx} out of order (expected {})", row + 1),
            ));
        }
        let base = f[1].as_bytes()[0];
        bases.push(base);
        let pair: usize = f[4]
            .parse()
            .map_err(|_| RnaStructError::parse("ct", "bad pair index"))?;
        if pair != 0 {
            if pair > n {
                return Err(RnaStructError::parse(
                    "ct",
                    format!("pair index {pair} exceeds length {n}"),
                ));
            }
            partner[row] = Some(pair - 1);
        }
    }
    if bases.len() != n {
        return Err(RnaStructError::parse(
            "ct",
            format!("header says {n} bases but found {}", bases.len()),
        ));
    }
    finish_record(title, &bases, partner, "ct")
}

/// Serialises a [`StructureRecord`] to ct format.
pub fn write_ct(record: &StructureRecord) -> String {
    let n = record.seq.len();
    let mut out = format!("{n}  {}\n", record.title);
    let bytes = record.seq.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        let prev = i; // 1-based (i-1)+1 == i; 0 for the first
        let next = if i + 1 < n { i + 2 } else { 0 };
        let pair = record.structure.partner(i).map(|p| p + 1).unwrap_or(0);
        // index  base  prev  next  pair  index2
        out.push_str(&format!(
            "{}  {}  {}  {}  {}  {}\n",
            i + 1,
            b as char,
            prev,
            next,
            pair,
            i + 1
        ));
    }
    out
}

// ---------------------------------------------------------------------
// bpseq format
// ---------------------------------------------------------------------

/// Parses a bpseq-format string into a [`StructureRecord`].
///
/// # Errors
/// [`RnaStructError::Parse`] on a malformed line, a non-RNA base, or
/// an out-of-range / asymmetric partner index.
pub fn read_bpseq(text: &str) -> Result<StructureRecord> {
    let mut bases: Vec<u8> = Vec::new();
    let mut pairs: Vec<usize> = Vec::new(); // 1-based, 0 = unpaired
    let mut title = String::new();

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('#') {
            // first comment line doubles as a title
            if title.is_empty() {
                title = rest.trim().to_string();
            }
            continue;
        }
        let f: Vec<&str> = trimmed.split_whitespace().collect();
        if f.len() < 3 {
            return Err(RnaStructError::parse(
                "bpseq",
                "data line has fewer than 3 fields",
            ));
        }
        let idx: usize = f[0]
            .parse()
            .map_err(|_| RnaStructError::parse("bpseq", "bad base index"))?;
        if idx != bases.len() + 1 {
            return Err(RnaStructError::parse(
                "bpseq",
                format!("base index {idx} out of order"),
            ));
        }
        bases.push(f[1].as_bytes()[0]);
        let pair: usize = f[2]
            .parse()
            .map_err(|_| RnaStructError::parse("bpseq", "bad pair index"))?;
        pairs.push(pair);
    }

    let n = bases.len();
    if n == 0 {
        return Err(RnaStructError::parse("bpseq", "no data lines"));
    }
    let mut partner: Vec<Option<usize>> = vec![None; n];
    for (i, &p) in pairs.iter().enumerate() {
        if p != 0 {
            if p > n {
                return Err(RnaStructError::parse(
                    "bpseq",
                    format!("pair index {p} exceeds length {n}"),
                ));
            }
            partner[i] = Some(p - 1);
        }
    }
    finish_record(title, &bases, partner, "bpseq")
}

/// Serialises a [`StructureRecord`] to bpseq format.
pub fn write_bpseq(record: &StructureRecord) -> String {
    let mut out = String::new();
    if !record.title.is_empty() {
        out.push_str(&format!("# {}\n", record.title));
    }
    let bytes = record.seq.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        let pair = record.structure.partner(i).map(|p| p + 1).unwrap_or(0);
        out.push_str(&format!("{} {} {}\n", i + 1, b as char, pair));
    }
    out
}

/// Shared tail of both readers: validates the partner array's
/// symmetry, builds the [`RnaSeq`] / [`Structure`], and assembles the
/// record.
fn finish_record(
    title: String,
    bases: &[u8],
    partner: Vec<Option<usize>>,
    format: &'static str,
) -> Result<StructureRecord> {
    // Symmetry: if i points to j, j must point back to i.
    for (i, &p) in partner.iter().enumerate() {
        if let Some(j) = p {
            if partner.get(j).copied().flatten() != Some(i) {
                return Err(RnaStructError::parse(
                    format,
                    format!("partner table is not symmetric at base {}", i + 1),
                ));
            }
        }
    }
    let seq = RnaSeq::parse(bases).map_err(|e| RnaStructError::parse(format, e.to_string()))?;
    let structure = Structure::from_partner(partner)
        .map_err(|e| RnaStructError::parse(format, e.to_string()))?;
    StructureRecord::new(title, seq, structure)
        .map_err(|e| RnaStructError::parse(format, e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_record() -> StructureRecord {
        let seq = RnaSeq::parse("GGGAAACCC").unwrap();
        let structure = Structure::from_dot_bracket("(((...)))").unwrap();
        StructureRecord::new("demo", seq, structure).unwrap()
    }

    #[test]
    fn ct_roundtrip() {
        let rec = sample_record();
        let text = write_ct(&rec);
        let back = read_ct(&text).unwrap();
        assert_eq!(back.seq, rec.seq);
        assert_eq!(back.structure, rec.structure);
        assert_eq!(back.title, "demo");
    }

    #[test]
    fn bpseq_roundtrip() {
        let rec = sample_record();
        let text = write_bpseq(&rec);
        let back = read_bpseq(&text).unwrap();
        assert_eq!(back.seq, rec.seq);
        assert_eq!(back.structure, rec.structure);
        assert_eq!(back.title, "demo");
    }

    #[test]
    fn ct_header_count_must_match() {
        let bad = "5  bad\n1 G 0 2 0 1\n2 C 1 0 0 2\n";
        assert!(read_ct(bad).is_err());
    }

    #[test]
    fn ct_rejects_length_exceeding_file_size() {
        // R32 L1: read_ct parsed `n` from the header with no upper
        // bound and eagerly did `vec![None; n]` (and
        // `Vec::with_capacity(n)`) BEFORE reading any data line. A tiny
        // file declaring a huge `n` (e.g. `999999999999  x`) attempts a
        // multi-GB allocation → OOM. A structure can't have more bases
        // than there are bytes on disk, so `n > text.len()` is rejected
        // up front. The error must name the file-size bound (proving the
        // cap fired before the later "found 0 bases" mismatch path).
        let text = "50  x\n"; // n=50 declared, file is 6 bytes
        let err = read_ct(text).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("file size") || msg.contains("exceeds"),
            "expected a file-size-bound rejection, got: {msg}"
        );
    }

    #[test]
    fn ct_rejects_asymmetric() {
        // base 1 pairs base 3 but base 3 pairs nobody
        let bad = "3  x\n1 G 0 2 3 1\n2 A 1 3 0 2\n3 C 2 0 0 3\n";
        assert!(read_ct(bad).is_err());
    }

    #[test]
    fn bpseq_comments_and_title() {
        let text = "# my title\n1 G 9\n2 G 8\n3 G 7\n4 A 0\n5 A 0\n6 A 0\n7 C 3\n8 C 2\n9 C 1\n";
        let rec = read_bpseq(text).unwrap();
        assert_eq!(rec.title, "my title");
        assert_eq!(rec.structure.n_pairs(), 3);
    }

    #[test]
    fn bpseq_rejects_garbage() {
        assert!(read_bpseq("not a bpseq file at all").is_err());
        assert!(read_bpseq("").is_err());
    }

    #[test]
    fn record_length_mismatch_rejected() {
        let seq = RnaSeq::parse("GGG").unwrap();
        let structure = Structure::empty(5);
        assert!(StructureRecord::new("x", seq, structure).is_err());
    }

    #[test]
    fn ct_rejects_empty_file() {
        assert!(read_ct("").is_err());
        assert!(read_ct("   \n  \n").is_err());
    }

    #[test]
    fn ct_rejects_non_numeric_header_length() {
        // The header's first field must parse as a count.
        assert!(read_ct("notanumber  title\n1 G 0 2 0 1\n").is_err());
    }

    #[test]
    fn ct_rejects_short_data_line() {
        // A data line needs at least 5 whitespace-separated fields.
        let bad = "2  x\n1 G 0 2\n2 C 1 0 0 2\n";
        assert!(read_ct(bad).is_err());
    }

    #[test]
    fn ct_rejects_out_of_order_index() {
        // Base indices must run 1, 2, 3, … in order.
        let bad = "2  x\n2 G 0 2 0 2\n1 C 1 0 0 1\n";
        assert!(read_ct(bad).is_err());
    }

    #[test]
    fn ct_rejects_pair_index_beyond_length() {
        // A partner index larger than the sequence length is invalid.
        let bad = "2  x\n1 G 0 2 9 1\n2 C 1 0 0 2\n";
        assert!(read_ct(bad).is_err());
    }

    #[test]
    fn ct_rejects_non_numeric_pair_field() {
        let bad = "2  x\n1 G 0 2 zz 1\n2 C 1 0 0 2\n";
        assert!(read_ct(bad).is_err());
    }

    #[test]
    fn ct_ignores_trailing_junk_lines() {
        // Lines beyond the header count are silently dropped.
        let text = "3  ok\n1 G 0 2 0 1\n2 A 1 3 0 2\n3 C 2 0 0 3\nGARBAGE TRAILER\n";
        let rec = read_ct(text).unwrap();
        assert_eq!(rec.seq.len(), 3);
    }

    #[test]
    fn bpseq_rejects_short_data_line() {
        // A bpseq data line needs at least 3 fields.
        assert!(read_bpseq("1 G\n").is_err());
    }

    #[test]
    fn bpseq_rejects_out_of_order_index() {
        assert!(read_bpseq("2 G 0\n1 C 0\n").is_err());
    }

    #[test]
    fn bpseq_rejects_pair_index_beyond_length() {
        let bad = "1 G 9\n2 C 0\n";
        assert!(read_bpseq(bad).is_err());
    }

    #[test]
    fn bpseq_write_omits_empty_title() {
        // The `write_bpseq` no-title branch: an empty title produces no
        // leading `#` comment line.
        let seq = RnaSeq::parse("GGGAAACCC").unwrap();
        let structure = Structure::from_dot_bracket("(((...)))").unwrap();
        let rec = StructureRecord::new("", seq, structure).unwrap();
        let text = write_bpseq(&rec);
        assert!(!text.starts_with('#'), "no title => no comment line");
        // It still round-trips.
        let back = read_bpseq(&text).unwrap();
        assert_eq!(back.structure.n_pairs(), 3);
        assert_eq!(back.title, "");
    }
}
