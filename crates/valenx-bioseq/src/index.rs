//! FASTA `.fai`-style indexing and sequence hashing / deduplication.
//!
//! [`FastaIndex`] builds the same record→offset table that `samtools
//! faidx` writes, enabling random access into a large FASTA without
//! reparsing it. The hashing helpers ([`sequence_hash`],
//! [`deduplicate`]) support collapsing identical sequences — common
//! when merging FASTA files or removing PCR duplicates by sequence.

use crate::error::{BioseqError, Result};
use crate::record::SeqRecord;
use std::collections::HashMap;

/// One `.fai` index entry — the five columns `samtools faidx` writes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FaiEntry {
    /// Record name (the FASTA id).
    pub name: String,
    /// Sequence length in residues.
    pub length: usize,
    /// Byte offset of the first residue of the sequence in the file.
    pub offset: usize,
    /// Number of residues per sequence line.
    pub line_bases: usize,
    /// Number of bytes per sequence line (residues + the newline).
    pub line_width: usize,
}

/// A `.fai`-style FASTA index — an ordered list of [`FaiEntry`] plus a
/// name→index lookup.
#[derive(Clone, Debug, Default)]
pub struct FastaIndex {
    entries: Vec<FaiEntry>,
    by_name: HashMap<String, usize>,
}

impl FastaIndex {
    /// Builds an index by scanning a FASTA string.
    ///
    /// Tracks byte offsets so the resulting index matches the on-disk
    /// `.fai` produced by `samtools faidx` for the same file (assuming
    /// a uniform line width per record — the standard FASTA layout).
    /// Returns [`BioseqError::Parse`] if the text starts with sequence
    /// data before any `>` header.
    pub fn build(text: &str) -> Result<FastaIndex> {
        let mut entries: Vec<FaiEntry> = Vec::new();
        let mut by_name: HashMap<String, usize> = HashMap::new();

        let bytes = text.as_bytes();
        let mut byte_pos = 0usize;
        let mut cur: Option<FaiEntry> = None;
        let mut first_line_seen = false;

        for line in text.split_inclusive('\n') {
            let trimmed = line.trim_end_matches(['\n', '\r']);
            let line_len = line.len();
            if let Some(rest) = trimmed.strip_prefix('>') {
                // Flush the previous record.
                if let Some(e) = cur.take() {
                    by_name.insert(e.name.clone(), entries.len());
                    entries.push(e);
                }
                let name = rest
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                cur = Some(FaiEntry {
                    name,
                    length: 0,
                    offset: byte_pos + line_len, // first seq byte
                    line_bases: 0,
                    line_width: 0,
                });
                first_line_seen = true;
            } else if !trimmed.is_empty() {
                if cur.is_none() && !first_line_seen {
                    return Err(BioseqError::parse(
                        "fasta",
                        "sequence data before first '>' header",
                    ));
                }
                if let Some(e) = cur.as_mut() {
                    if e.line_bases == 0 {
                        e.line_bases = trimmed.len();
                        e.line_width = line_len;
                    }
                    e.length += trimmed.len();
                }
            }
            byte_pos += line_len;
        }
        if let Some(e) = cur.take() {
            by_name.insert(e.name.clone(), entries.len());
            entries.push(e);
        }
        let _ = bytes;
        Ok(FastaIndex { entries, by_name })
    }

    /// Number of indexed records.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// `true` if the index has no records.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// All entries, in file order.
    pub fn entries(&self) -> &[FaiEntry] {
        &self.entries
    }

    /// Looks up an entry by record name.
    pub fn get(&self, name: &str) -> Option<&FaiEntry> {
        self.by_name.get(name).map(|&i| &self.entries[i])
    }

    /// Serializes the index to the tab-separated `.fai` text format
    /// (`name\tlength\toffset\tlinebases\tlinewidth`).
    pub fn to_fai_string(&self) -> String {
        let mut out = String::new();
        for e in &self.entries {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                e.name, e.length, e.offset, e.line_bases, e.line_width
            ));
        }
        out
    }

    /// Parses a `.fai` text file back into a [`FastaIndex`].
    ///
    /// Returns [`BioseqError::Parse`] on a malformed line.
    pub fn from_fai_string(text: &str) -> Result<FastaIndex> {
        let mut entries = Vec::new();
        let mut by_name = HashMap::new();
        for (lineno, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() != 5 {
                return Err(BioseqError::parse(
                    "fai",
                    format!("line {} has {} columns, expected 5", lineno + 1, cols.len()),
                ));
            }
            let parse_num = |s: &str, what: &'static str| -> Result<usize> {
                s.parse::<usize>()
                    .map_err(|_| BioseqError::parse("fai", format!("bad {what}: `{s}`")))
            };
            let entry = FaiEntry {
                name: cols[0].to_string(),
                length: parse_num(cols[1], "length")?,
                offset: parse_num(cols[2], "offset")?,
                line_bases: parse_num(cols[3], "line_bases")?,
                line_width: parse_num(cols[4], "line_width")?,
            };
            by_name.insert(entry.name.clone(), entries.len());
            entries.push(entry);
        }
        Ok(FastaIndex { entries, by_name })
    }
}

/// A 64-bit FNV-1a hash of a sequence's residues (case-insensitive).
///
/// FNV-1a is a fast, well-distributed non-cryptographic hash — exactly
/// what is needed to bucket identical sequences for deduplication.
pub fn sequence_hash(seq: &crate::seq::Seq) -> u64 {
    fnv1a(seq.as_bytes())
}

/// FNV-1a over a byte slice (residues uppercased so case does not
/// affect the hash).
fn fnv1a(bytes: &[u8]) -> u64 {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut h = OFFSET;
    for &b in bytes {
        h ^= b.to_ascii_uppercase() as u64;
        h = h.wrapping_mul(PRIME);
    }
    h
}

/// Deduplicates a slice of records by sequence content.
///
/// Returns the records whose sequence is seen for the first time, in
/// input order. Two records with identical residues (ignoring id /
/// description / case) collapse to the first occurrence.
pub fn deduplicate(records: &[SeqRecord]) -> Vec<SeqRecord> {
    let mut seen: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut out: Vec<SeqRecord> = Vec::new();
    for rec in records {
        let h = sequence_hash(&rec.seq);
        let bucket = seen.entry(h).or_default();
        // Guard against hash collisions: confirm byte equality.
        let is_dup = bucket
            .iter()
            .any(|&i| out[i].seq.as_bytes() == rec.seq.as_bytes());
        if !is_dup {
            bucket.push(out.len());
            out.push(rec.clone());
        }
    }
    out
}

/// Groups records by identical sequence content. Returns one `Vec` of
/// records per distinct sequence; each inner `Vec` preserves input
/// order. Useful for reporting which ids share a sequence.
pub fn group_by_sequence(records: &[SeqRecord]) -> Vec<Vec<SeqRecord>> {
    let mut groups: Vec<Vec<SeqRecord>> = Vec::new();
    let mut index: HashMap<u64, Vec<usize>> = HashMap::new();
    for rec in records {
        let h = sequence_hash(&rec.seq);
        let candidates = index.entry(h).or_default();
        let mut placed = false;
        for &gi in candidates.iter() {
            if groups[gi][0].seq.as_bytes() == rec.seq.as_bytes() {
                groups[gi].push(rec.clone());
                placed = true;
                break;
            }
        }
        if !placed {
            candidates.push(groups.len());
            groups.push(vec![rec.clone()]);
        }
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{Seq, SeqKind};

    const FASTA: &str = ">seq1 first\nACGTACGT\nACGT\n>seq2 second\nTTTTTTTT\n";

    #[test]
    fn build_index_counts_records() {
        let idx = FastaIndex::build(FASTA).unwrap();
        assert_eq!(idx.len(), 2);
        assert!(!idx.is_empty());
    }

    #[test]
    fn index_entry_fields() {
        let idx = FastaIndex::build(FASTA).unwrap();
        let e1 = idx.get("seq1").unwrap();
        assert_eq!(e1.name, "seq1");
        assert_eq!(e1.length, 12); // ACGTACGT + ACGT
        assert_eq!(e1.line_bases, 8); // first line had 8 residues
        // Offset: ">seq1 first\n" is 12 bytes -> first residue at 12.
        assert_eq!(e1.offset, 12);
        let e2 = idx.get("seq2").unwrap();
        assert_eq!(e2.length, 8);
    }

    #[test]
    fn data_before_header_errors() {
        assert!(FastaIndex::build("ACGT\n>x\nACGT\n").is_err());
    }

    #[test]
    fn fai_string_roundtrip() {
        let idx = FastaIndex::build(FASTA).unwrap();
        let fai = idx.to_fai_string();
        let parsed = FastaIndex::from_fai_string(&fai).unwrap();
        assert_eq!(parsed.len(), idx.len());
        assert_eq!(parsed.get("seq1"), idx.get("seq1"));
        assert_eq!(parsed.get("seq2"), idx.get("seq2"));
    }

    #[test]
    fn malformed_fai_errors() {
        assert!(FastaIndex::from_fai_string("only\ttwo\n").is_err());
        assert!(FastaIndex::from_fai_string("n\tnotanumber\t0\t8\t9\n").is_err());
    }

    #[test]
    fn hash_is_case_insensitive_and_stable() {
        let a = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let b = Seq::new(SeqKind::Dna, "acgt").unwrap();
        assert_eq!(sequence_hash(&a), sequence_hash(&b));
        // Different content -> (almost surely) different hash.
        let c = Seq::new(SeqKind::Dna, "TTTT").unwrap();
        assert_ne!(sequence_hash(&a), sequence_hash(&c));
    }

    #[test]
    fn deduplicate_collapses_identical_sequences() {
        let recs = vec![
            SeqRecord::new("a", Seq::new(SeqKind::Dna, "ACGT").unwrap()),
            SeqRecord::new("b", Seq::new(SeqKind::Dna, "ACGT").unwrap()),
            SeqRecord::new("c", Seq::new(SeqKind::Dna, "TTTT").unwrap()),
        ];
        let dedup = deduplicate(&recs);
        // a and b share a sequence -> 2 unique records.
        assert_eq!(dedup.len(), 2);
        assert_eq!(dedup[0].id, "a");
        assert_eq!(dedup[1].id, "c");
    }

    #[test]
    fn group_by_sequence_buckets_records() {
        let recs = vec![
            SeqRecord::new("a", Seq::new(SeqKind::Dna, "ACGT").unwrap()),
            SeqRecord::new("b", Seq::new(SeqKind::Dna, "TTTT").unwrap()),
            SeqRecord::new("c", Seq::new(SeqKind::Dna, "ACGT").unwrap()),
        ];
        let groups = group_by_sequence(&recs);
        assert_eq!(groups.len(), 2);
        // The ACGT group has a and c.
        let acgt = groups.iter().find(|g| g.len() == 2).unwrap();
        assert_eq!(acgt[0].id, "a");
        assert_eq!(acgt[1].id, "c");
    }
}
