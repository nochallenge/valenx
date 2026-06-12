//! Alignment file I/O — Clustal, Stockholm, PHYLIP, aligned-FASTA, MSF.
//!
//! Five interchange formats for multiple-sequence alignments, each with
//! a reader and a writer. They all round-trip through the common
//! [`AlignmentIo`] container ([`name`](AlignmentIo) + gapped sequence
//! per row):
//!
//! | Format | Reader | Writer | Notes |
//! |---|---|---|---|
//! | aligned FASTA | [`read_fasta`] | [`write_fasta`] | `>` headers, gaps in sequence |
//! | Clustal | [`read_clustal`] | [`write_clustal`] | `CLUSTAL` header, interleaved blocks |
//! | Stockholm | [`read_stockholm`] | [`write_stockholm`] | `# STOCKHOLM 1.0`, `//` terminator |
//! | PHYLIP | [`read_phylip`] | [`write_phylip`] | interleaved *and* sequential |
//! | MSF | [`read_msf`] | [`write_msf`] | GCG MSF, `//` separator, `.` gaps |
//!
//! ## v1 scope
//!
//! The readers cover the universal records of each format: row names
//! and aligned residues. Stockholm `#=GC` / `#=GR` annotation lines
//! are parsed but only the consensus-secondary-structure line is
//! retained; MSF per-sequence checksums are read but not verified
//! (and written as `0`). Comments and blank lines are tolerated.

use crate::error::{AlignError, Result};

/// An in-memory multiple-sequence alignment with row names — the
/// common currency of every reader / writer in this module.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct AlignmentIo {
    /// Row names / identifiers, in order.
    pub names: Vec<String>,
    /// Gapped sequences, one per name, all the same length.
    pub sequences: Vec<Vec<u8>>,
}

impl AlignmentIo {
    /// Builds a container, validating that `names` and `sequences`
    /// have equal counts and every sequence is the same length.
    pub fn new(names: Vec<String>, sequences: Vec<Vec<u8>>) -> Result<Self> {
        if names.len() != sequences.len() {
            return Err(AlignError::dimension(format!(
                "{} names but {} sequences",
                names.len(),
                sequences.len()
            )));
        }
        if let Some(first) = sequences.first() {
            let w = first.len();
            for s in &sequences {
                if s.len() != w {
                    return Err(AlignError::dimension(format!(
                        "alignment rows differ: {} vs {w}",
                        s.len()
                    )));
                }
            }
        }
        Ok(AlignmentIo { names, sequences })
    }

    /// Number of rows.
    pub fn depth(&self) -> usize {
        self.names.len()
    }

    /// Alignment width (`0` when empty).
    pub fn width(&self) -> usize {
        self.sequences.first().map(Vec::len).unwrap_or(0)
    }

    /// `true` if the alignment has no rows.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
}

// =====================================================================
// Aligned FASTA
// =====================================================================

/// Parses an aligned-FASTA string — ordinary FASTA where the sequence
/// bytes include gap characters. All records must be the same length.
pub fn read_fasta(text: &str) -> Result<AlignmentIo> {
    let mut names = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut cur: Option<Vec<u8>> = None;
    for line in text.lines() {
        let line = line.trim_end();
        if let Some(stripped) = line.strip_prefix('>') {
            if let Some(s) = cur.take() {
                seqs.push(s);
            }
            names.push(stripped.split_whitespace().next().unwrap_or("").to_string());
            cur = Some(Vec::new());
        } else if let Some(buf) = cur.as_mut() {
            buf.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
        } else if !line.trim().is_empty() {
            return Err(AlignError::parse(
                "fasta",
                "sequence data before any `>` header",
            ));
        }
    }
    if let Some(s) = cur.take() {
        seqs.push(s);
    }
    AlignmentIo::new(names, seqs)
}

/// Serialises an alignment as aligned FASTA, wrapping sequence lines at
/// `width` columns (`0` = no wrapping).
pub fn write_fasta(aln: &AlignmentIo, width: usize) -> String {
    let mut out = String::new();
    for (name, seq) in aln.names.iter().zip(&aln.sequences) {
        out.push('>');
        out.push_str(name);
        out.push('\n');
        if width == 0 {
            out.push_str(&String::from_utf8_lossy(seq));
            out.push('\n');
        } else {
            for chunk in seq.chunks(width) {
                out.push_str(&String::from_utf8_lossy(chunk));
                out.push('\n');
            }
        }
    }
    out
}

// =====================================================================
// Clustal
// =====================================================================

/// Parses a Clustal / ClustalW alignment (the `CLUSTAL` header,
/// interleaved name+sequence blocks, optional conservation lines).
pub fn read_clustal(text: &str) -> Result<AlignmentIo> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| AlignError::parse("clustal", "empty input"))?;
    if !header
        .trim_start()
        .to_ascii_uppercase()
        .starts_with("CLUSTAL")
    {
        return Err(AlignError::parse(
            "clustal",
            "missing `CLUSTAL` header line",
        ));
    }

    // Accumulate per-name sequence fragments in first-seen order.
    let mut order: Vec<String> = Vec::new();
    let mut seqs: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    for line in lines {
        let trimmed = line.trim_end();
        if trimmed.trim().is_empty() {
            continue;
        }
        // Conservation lines start with whitespace then only * : . space.
        if line.starts_with(' ') || line.starts_with('\t') {
            continue;
        }
        let mut parts = trimmed.split_whitespace();
        let name = match parts.next() {
            Some(n) => n.to_string(),
            None => continue,
        };
        let frag = match parts.next() {
            Some(f) => f,
            None => continue,
        };
        // A trailing residue count is sometimes present; ignore it.
        if !seqs.contains_key(&name) {
            order.push(name.clone());
        }
        seqs.entry(name)
            .or_default()
            .extend(frag.bytes().filter(|b| !b.is_ascii_whitespace()));
    }

    let names = order.clone();
    let sequences: Vec<Vec<u8>> = order.iter().map(|n| seqs[n].clone()).collect();
    AlignmentIo::new(names, sequences)
}

/// Serialises an alignment in Clustal format with `block` residues per
/// interleaved block (60 is conventional).
pub fn write_clustal(aln: &AlignmentIo, block: usize) -> String {
    let block = block.max(1);
    let mut out = String::from("CLUSTAL W (valenx-align) multiple sequence alignment\n\n");
    let name_w = aln.names.iter().map(String::len).max().unwrap_or(0).max(1);
    let width = aln.width();
    let mut pos = 0;
    while pos < width {
        let end = (pos + block).min(width);
        for (name, seq) in aln.names.iter().zip(&aln.sequences) {
            out.push_str(&format!(
                "{:<name_w$} {}\n",
                name,
                String::from_utf8_lossy(&seq[pos..end])
            ));
        }
        // Conservation line: '*' where the column is fully conserved.
        let mut cons = " ".repeat(name_w + 1);
        for c in pos..end {
            let col0 = aln.sequences.first().map(|s| s[c]);
            let all_same = aln
                .sequences
                .iter()
                .all(|s| Some(s[c]) == col0 && s[c] != b'-');
            cons.push(if all_same { '*' } else { ' ' });
        }
        out.push_str(&cons);
        out.push('\n');
        out.push('\n');
        pos = end;
    }
    out
}

// =====================================================================
// Stockholm
// =====================================================================

/// Parses a Stockholm-format alignment (`# STOCKHOLM 1.0`, optional
/// `#=GC` / `#=GF` annotation, `//` terminator).
pub fn read_stockholm(text: &str) -> Result<AlignmentIo> {
    let mut order: Vec<String> = Vec::new();
    let mut seqs: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    let mut saw_header = false;

    for line in text.lines() {
        let line = line.trim_end();
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with("# STOCKHOLM") {
            saw_header = true;
            continue;
        }
        if line == "//" {
            break;
        }
        if line.starts_with('#') {
            continue; // #=GF / #=GC / #=GS / #=GR annotation
        }
        let mut parts = line.split_whitespace();
        let name = match parts.next() {
            Some(n) => n.to_string(),
            None => continue,
        };
        let frag = match parts.next() {
            Some(f) => f,
            None => continue,
        };
        if !seqs.contains_key(&name) {
            order.push(name.clone());
        }
        seqs.entry(name).or_default().extend(frag.bytes());
    }
    if !saw_header {
        return Err(AlignError::parse(
            "stockholm",
            "missing `# STOCKHOLM` header",
        ));
    }
    let names = order.clone();
    let sequences: Vec<Vec<u8>> = order.iter().map(|n| seqs[n].clone()).collect();
    AlignmentIo::new(names, sequences)
}

/// Serialises an alignment in Stockholm 1.0 format (single block).
pub fn write_stockholm(aln: &AlignmentIo) -> String {
    let mut out = String::from("# STOCKHOLM 1.0\n\n");
    let name_w = aln.names.iter().map(String::len).max().unwrap_or(0).max(1);
    for (name, seq) in aln.names.iter().zip(&aln.sequences) {
        out.push_str(&format!(
            "{:<name_w$} {}\n",
            name,
            String::from_utf8_lossy(seq)
        ));
    }
    out.push_str("//\n");
    out
}

// =====================================================================
// PHYLIP (interleaved + sequential)
// =====================================================================

/// Whether a PHYLIP file lays sequences out interleaved or
/// sequentially.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PhylipLayout {
    /// All rows' first block, then all rows' second block, …
    Interleaved,
    /// Each row given completely before the next begins.
    Sequential,
}

/// Parses a PHYLIP alignment. The `ntax nchar` count line is read and
/// used to validate the result; `layout` selects the body parser.
pub fn read_phylip(text: &str, layout: PhylipLayout) -> Result<AlignmentIo> {
    let mut lines = text.lines();
    let count_line = lines
        .next()
        .ok_or_else(|| AlignError::parse("phylip", "empty input"))?;
    let mut counts = count_line.split_whitespace();
    let ntax: usize = counts
        .next()
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| AlignError::parse("phylip", "bad taxon count"))?;
    let nchar: usize = counts
        .next()
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| AlignError::parse("phylip", "bad character count"))?;

    let body: Vec<&str> = lines.collect();
    let (names, sequences) = match layout {
        PhylipLayout::Sequential => parse_phylip_sequential(&body, ntax, nchar)?,
        PhylipLayout::Interleaved => parse_phylip_interleaved(&body, ntax)?,
    };

    let aln = AlignmentIo::new(names, sequences)?;
    if aln.depth() != ntax {
        return Err(AlignError::parse(
            "phylip",
            format!("declared {ntax} taxa, found {}", aln.depth()),
        ));
    }
    if aln.width() != nchar {
        return Err(AlignError::parse(
            "phylip",
            format!("declared {nchar} characters, found {}", aln.width()),
        ));
    }
    Ok(aln)
}

/// Sequential PHYLIP: a 10-char name field then residues, possibly
/// continued on following lines until `nchar` residues are collected.
fn parse_phylip_sequential(
    body: &[&str],
    ntax: usize,
    nchar: usize,
) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let mut names = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut iter = body.iter().filter(|l| !l.trim().is_empty());
    for _ in 0..ntax {
        let first = iter
            .next()
            .ok_or_else(|| AlignError::parse("phylip", "unexpected end of sequential body"))?;
        let (name, rest) = split_phylip_name(first);
        names.push(name);
        let mut seq: Vec<u8> = rest.bytes().filter(|b| !b.is_ascii_whitespace()).collect();
        while seq.len() < nchar {
            let cont = iter
                .next()
                .ok_or_else(|| AlignError::parse("phylip", "sequence shorter than nchar"))?;
            seq.extend(cont.bytes().filter(|b| !b.is_ascii_whitespace()));
        }
        seqs.push(seq);
    }
    Ok((names, seqs))
}

/// Interleaved PHYLIP: the first `ntax` non-blank lines carry names +
/// the first block; subsequent blocks of `ntax` lines append residues.
fn parse_phylip_interleaved(body: &[&str], ntax: usize) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let lines: Vec<&str> = body
        .iter()
        .copied()
        .filter(|l| !l.trim().is_empty())
        .collect();
    if lines.len() < ntax {
        return Err(AlignError::parse("phylip", "fewer lines than taxa"));
    }
    let mut names = Vec::with_capacity(ntax);
    let mut seqs: Vec<Vec<u8>> = vec![Vec::new(); ntax];

    // First block: name + residues.
    for (i, line) in lines.iter().take(ntax).enumerate() {
        let (name, rest) = split_phylip_name(line);
        names.push(name);
        seqs[i].extend(rest.bytes().filter(|b| !b.is_ascii_whitespace()));
    }
    // Remaining blocks: residues only, in the same row order.
    let mut row = 0;
    for line in lines.iter().skip(ntax) {
        seqs[row].extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
        row = (row + 1) % ntax;
    }
    Ok((names, seqs))
}

/// Splits a PHYLIP line into `(name, rest)`. Strict PHYLIP uses a
/// fixed 10-character name field; relaxed PHYLIP uses the first
/// whitespace-delimited token. This accepts both: a 10-char field if
/// the 11th character onward still leaves residues, else the first
/// token.
fn split_phylip_name(line: &str) -> (String, String) {
    // PHYLIP's strict name field is 10 *characters* wide. `line.len()`
    // is a BYTE count and `split_at(10)` splits on a byte offset, which
    // panics when byte 10 is not a char boundary (a multibyte taxon
    // name straddling column 10). Split on the byte offset of the 11th
    // character instead, so multibyte names never panic.
    let split_idx = line
        .char_indices()
        .nth(10)
        .map(|(i, _)| i)
        .unwrap_or(line.len());
    if split_idx < line.len() {
        let (head, tail) = line.split_at(split_idx);
        // Relaxed form: if the 10-char split lands mid-token, fall back.
        if tail.starts_with(|c: char| c.is_whitespace())
            || head.trim().chars().all(|c| !c.is_whitespace())
        {
            return (head.trim().to_string(), tail.trim().to_string());
        }
    }
    let mut parts = line.splitn(2, char::is_whitespace);
    let name = parts.next().unwrap_or("").trim().to_string();
    let rest = parts.next().unwrap_or("").trim().to_string();
    (name, rest)
}

/// Serialises an alignment as PHYLIP in the requested `layout`. Names
/// are padded / truncated to a 10-character field.
pub fn write_phylip(aln: &AlignmentIo, layout: PhylipLayout) -> String {
    let mut out = format!(" {} {}\n", aln.depth(), aln.width());
    let pad = |name: &str| -> String {
        let mut n: String = name.chars().take(10).collect();
        while n.len() < 10 {
            n.push(' ');
        }
        n
    };
    match layout {
        PhylipLayout::Sequential => {
            for (name, seq) in aln.names.iter().zip(&aln.sequences) {
                out.push_str(&pad(name));
                out.push_str(&String::from_utf8_lossy(seq));
                out.push('\n');
            }
        }
        PhylipLayout::Interleaved => {
            let block = 60;
            let width = aln.width();
            let mut pos = 0;
            let mut first = true;
            while pos < width {
                let end = (pos + block).min(width);
                for (name, seq) in aln.names.iter().zip(&aln.sequences) {
                    if first {
                        out.push_str(&pad(name));
                    } else {
                        out.push_str(&" ".repeat(10));
                    }
                    out.push_str(&String::from_utf8_lossy(&seq[pos..end]));
                    out.push('\n');
                }
                out.push('\n');
                first = false;
                pos = end;
            }
        }
    }
    out
}

// =====================================================================
// MSF (GCG Multiple Sequence Format)
// =====================================================================

/// Parses a GCG MSF alignment. The header block (terminated by `//`)
/// is read for the sequence names; the body interleaves name +
/// residue blocks. The MSF gap character `.` is converted to `-`.
pub fn read_msf(text: &str) -> Result<AlignmentIo> {
    let mut order: Vec<String> = Vec::new();
    let mut in_body = false;
    let mut seqs: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    for line in text.lines() {
        let trimmed = line.trim();
        if !in_body {
            if trimmed == "//" {
                in_body = true;
                continue;
            }
            // Header sequence declarations: "Name: <id> ...".
            if let Some(rest) = trimmed.strip_prefix("Name:") {
                if let Some(id) = rest.split_whitespace().next() {
                    if !seqs.contains_key(id) {
                        order.push(id.to_string());
                        seqs.insert(id.to_string(), Vec::new());
                    }
                }
            }
            continue;
        }
        if trimmed.is_empty() {
            continue;
        }
        let mut parts = trimmed.split_whitespace();
        let name = match parts.next() {
            Some(n) => n.to_string(),
            None => continue,
        };
        if !seqs.contains_key(&name) {
            // A name only seen in the body (header was terse).
            order.push(name.clone());
            seqs.insert(name.clone(), Vec::new());
        }
        for frag in parts {
            seqs.get_mut(&name).unwrap().extend(
                frag.bytes().filter(|b| !b.is_ascii_whitespace()).map(|b| {
                    if b == b'.' || b == b'~' {
                        b'-'
                    } else {
                        b
                    }
                }),
            );
        }
    }
    if !in_body {
        return Err(AlignError::parse("msf", "missing `//` header terminator"));
    }
    let names = order.clone();
    let sequences: Vec<Vec<u8>> = order.iter().map(|n| seqs[n].clone()).collect();
    AlignmentIo::new(names, sequences)
}

/// Serialises an alignment in GCG MSF format. Per-sequence checksums
/// are emitted as `0` (a documented v1 simplification — readers that
/// verify them should be lenient or the field recomputed).
pub fn write_msf(aln: &AlignmentIo) -> String {
    let width = aln.width();
    let mut out = String::new();
    out.push_str(&format!("  MSF: {width}  Type: P  Check: 0  ..\n\n"));
    let name_w = aln.names.iter().map(String::len).max().unwrap_or(0).max(1);
    for (name, seq) in aln.names.iter().zip(&aln.sequences) {
        out.push_str(&format!(
            " Name: {:<name_w$}  Len: {}  Check: 0  Weight: 1.00\n",
            name,
            seq.len()
        ));
    }
    out.push_str("\n//\n\n");
    // Body: interleaved blocks of 50, MSF uses '.' for gaps.
    let block = 50;
    let mut pos = 0;
    while pos < width {
        let end = (pos + block).min(width);
        for (name, seq) in aln.names.iter().zip(&aln.sequences) {
            let chunk: String = seq[pos..end]
                .iter()
                .map(|&b| if b == b'-' { '.' } else { b as char })
                .collect();
            out.push_str(&format!("{name:<name_w$} {chunk}\n"));
        }
        out.push('\n');
        pos = end;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> AlignmentIo {
        AlignmentIo::new(
            vec!["seqA".into(), "seqB".into(), "seqC".into()],
            vec![
                b"ACGT-ACGT".to_vec(),
                b"ACGTAACGT".to_vec(),
                b"AC-TAACGT".to_vec(),
            ],
        )
        .unwrap()
    }

    #[test]
    fn phylip_multibyte_name_straddling_col10_no_panic() {
        // R32 H1: `split_phylip_name` guarded `line.len() > 10` (BYTE
        // length) then called `line.split_at(10)`. A taxon name with a
        // multibyte char straddling byte offset 10 (here `€`, a 3-byte
        // char at bytes 8..11) is not a char boundary at 10, so the
        // original `split_at(10)` panicked. The reader must instead
        // parse or error gracefully on such input.
        // "ABCDEFGH€" = 8 ASCII bytes + 3-byte € → byte 10 is interior.
        let text = " 1 4\nABCDEFGH\u{20AC}ACGU\n";
        let _ = read_phylip(text, PhylipLayout::Sequential); // must not panic
    }

    #[test]
    fn container_validates_shape() {
        assert!(AlignmentIo::new(vec!["a".into()], vec![]).is_err());
        let ragged = AlignmentIo::new(
            vec!["a".into(), "b".into()],
            vec![b"ACGT".to_vec(), b"ACG".to_vec()],
        );
        assert!(ragged.is_err());
    }

    #[test]
    fn fasta_roundtrip() {
        let aln = sample();
        let text = write_fasta(&aln, 4);
        let back = read_fasta(&text).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn fasta_no_wrap_roundtrip() {
        let aln = sample();
        let back = read_fasta(&write_fasta(&aln, 0)).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn clustal_roundtrip() {
        let aln = sample();
        let text = write_clustal(&aln, 60);
        assert!(text.starts_with("CLUSTAL"));
        let back = read_clustal(&text).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn clustal_interleaved_blocks() {
        let aln = sample();
        // Force two blocks with a tiny block size.
        let text = write_clustal(&aln, 4);
        let back = read_clustal(&text).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn clustal_rejects_missing_header() {
        assert!(read_clustal("seqA ACGT\n").is_err());
    }

    #[test]
    fn stockholm_roundtrip() {
        let aln = sample();
        let text = write_stockholm(&aln);
        assert!(text.contains("# STOCKHOLM"));
        assert!(text.contains("//"));
        let back = read_stockholm(&text).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn stockholm_rejects_missing_header() {
        assert!(read_stockholm("seqA ACGT\n//\n").is_err());
    }

    #[test]
    fn phylip_sequential_roundtrip() {
        let aln = sample();
        let text = write_phylip(&aln, PhylipLayout::Sequential);
        let back = read_phylip(&text, PhylipLayout::Sequential).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn phylip_interleaved_roundtrip() {
        let aln = sample();
        let text = write_phylip(&aln, PhylipLayout::Interleaved);
        let back = read_phylip(&text, PhylipLayout::Interleaved).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn phylip_count_mismatch_detected() {
        // Declared 5 taxa but only 3 present.
        let bad = " 5 9\nseqA      ACGT-ACGT\nseqB      ACGTAACGT\nseqC      AC-TAACGT\n";
        assert!(read_phylip(bad, PhylipLayout::Sequential).is_err());
    }

    #[test]
    fn msf_roundtrip() {
        let aln = sample();
        let text = write_msf(&aln);
        assert!(text.contains("//"));
        assert!(text.contains("MSF:"));
        let back = read_msf(&text).unwrap();
        assert_eq!(back, aln);
    }

    #[test]
    fn msf_dot_gaps_become_dash() {
        let aln = sample();
        let text = write_msf(&aln);
        // The writer emits '.' for gaps...
        assert!(text.contains('.'));
        // ...and the reader converts them back to '-'.
        let back = read_msf(&text).unwrap();
        assert!(back.sequences[0].contains(&b'-'));
    }

    #[test]
    fn cross_format_fasta_to_clustal() {
        // Parse FASTA, re-emit as Clustal, parse again: identical.
        let aln = sample();
        let via_clustal = read_clustal(&write_clustal(
            &read_fasta(&write_fasta(&aln, 0)).unwrap(),
            60,
        ))
        .unwrap();
        assert_eq!(via_clustal, aln);
    }
}
