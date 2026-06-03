//! GenBank flat-file reader and writer.
//!
//! Parses the INSDC GenBank format — `LOCUS` / `DEFINITION` /
//! `ACCESSION` / `VERSION` / `KEYWORDS` / `SOURCE` / `ORGANISM` /
//! `REFERENCE` / `FEATURES` / `ORIGIN` — into a [`SeqRecord`], and
//! emits a valid GenBank flat file from one.
//!
//! Feature locations (including `complement(...)`, `join(...)`,
//! `order(...)`, and `n^n+1` between-bases) are handled by
//! [`crate::io::locstr`]. Feature qualifiers may span multiple lines
//! and may contain embedded `=` and `/` characters inside quoted
//! values — the parser only treats a `/` at the qualifier column
//! (21) as the start of a new qualifier.
//!
//! `REFERENCE` blocks are parsed into [`Reference`] entries and
//! round-tripped on write. The literature fields are preserved as
//! free text (AUTHORS / CONSRTM / TITLE / JOURNAL / PUBMED / MEDLINE
//! / REMARK).

use crate::alphabet::SeqKind;
use crate::error::{BioseqError, Result};
use crate::io::locstr;
use crate::record::{Reference, SeqFeature, SeqRecord};
use crate::seq::{Seq, Topology};
use std::collections::BTreeMap;

/// Parses a GenBank flat file holding a single record.
///
/// The `LOCUS` line's molecule-type field selects the [`SeqKind`]
/// (`DNA`/`RNA`/`AA`) and the topology (`circular`/`linear`). Returns
/// [`BioseqError::Parse`] if no `ORIGIN` sequence is found.
pub fn parse(text: &str) -> Result<SeqRecord> {
    let lines: Vec<&str> = text.lines().collect();
    let mut id = String::new();
    let mut definition = String::new();
    let mut accession = String::new();
    let mut version = String::new();
    let mut keywords = String::new();
    let mut kind = SeqKind::Dna;
    let mut topology = Topology::Linear;
    let mut annotations: BTreeMap<String, String> = BTreeMap::new();
    let mut features: Vec<SeqFeature> = Vec::new();
    let mut references: Vec<Reference> = Vec::new();
    let mut origin = String::new();

    let mut i = 0;
    while i < lines.len() {
        let line = lines[i];
        if line.starts_with("LOCUS") {
            let (k, t, locus_id) = parse_locus(line);
            kind = k;
            topology = t;
            id = locus_id;
            i += 1;
        } else if line.starts_with("DEFINITION") {
            let (text, next) = collect_continuation(&lines, i, 12);
            definition = strip_field("DEFINITION", &text);
            i = next;
        } else if let Some(rest) = line.strip_prefix("ACCESSION") {
            accession = rest.trim().to_string();
            i += 1;
        } else if let Some(rest) = line.strip_prefix("VERSION") {
            version = rest.trim().to_string();
            i += 1;
        } else if line.starts_with("KEYWORDS") {
            let (text, next) = collect_continuation(&lines, i, 12);
            keywords = strip_field("KEYWORDS", &text).trim_end_matches('.').to_string();
            i = next;
        } else if let Some(rest) = line.strip_prefix("SOURCE") {
            let src = rest.trim().to_string();
            if !src.is_empty() {
                annotations.insert("source".to_string(), src);
            }
            i += 1;
        } else if let Some(rest) = line.trim_start().strip_prefix("ORGANISM") {
            let org = rest.trim().to_string();
            if !org.is_empty() {
                annotations.insert("organism".to_string(), org);
            }
            i += 1;
        } else if line.starts_with("REFERENCE") {
            let (reference, next) = parse_reference(&lines, i);
            references.push(reference);
            i = next;
        } else if line.starts_with("FEATURES") {
            let (parsed, next) = parse_features(&lines, i + 1)?;
            features = parsed;
            i = next;
        } else if line.starts_with("ORIGIN") {
            let (seq, next) = parse_origin(&lines, i + 1);
            origin = seq;
            i = next;
        } else if line.starts_with("//") {
            break;
        } else {
            i += 1;
        }
    }

    if origin.is_empty() {
        return Err(BioseqError::parse(
            "genbank",
            "no ORIGIN sequence block found",
        ));
    }
    if !accession.is_empty() {
        annotations.insert("accession".to_string(), accession.clone());
    }
    if !version.is_empty() {
        annotations.insert("version".to_string(), version);
    }
    if !keywords.is_empty() {
        annotations.insert("keywords".to_string(), keywords);
    }
    let final_id = if !accession.is_empty() {
        accession
    } else {
        id.clone()
    };
    let seq = Seq::with_topology(kind, origin, topology)?;
    Ok(SeqRecord {
        id: final_id,
        name: id,
        description: definition,
        seq,
        features,
        annotations,
        references,
    })
}

/// Parses the `LOCUS` line → `(kind, topology, locus_name)`.
fn parse_locus(line: &str) -> (SeqKind, Topology, String) {
    let rest = &line["LOCUS".len()..];
    let toks: Vec<&str> = rest.split_whitespace().collect();
    let name = toks.first().copied().unwrap_or("").to_string();
    let mut kind = SeqKind::Dna;
    let mut topology = Topology::Linear;
    for t in &toks {
        let lower = t.to_ascii_lowercase();
        if lower == "circular" {
            topology = Topology::Circular;
        } else if lower == "linear" {
            topology = Topology::Linear;
        } else if lower.contains("rna") {
            kind = SeqKind::Rna;
        } else if lower == "aa" || lower == "prt" {
            kind = SeqKind::Protein;
        } else if lower.contains("dna") {
            kind = SeqKind::Dna;
        }
    }
    (kind, topology, name)
}

/// Collects a field plus any continuation lines (indented to `indent`).
/// Returns the joined text and the index of the next non-continuation
/// line.
fn collect_continuation(lines: &[&str], start: usize, indent: usize) -> (String, usize) {
    let mut text = lines[start].to_string();
    let mut i = start + 1;
    while i < lines.len() {
        let l = lines[i];
        // A continuation line is blank-prefixed at least `indent` cols
        // and not itself a new top-level keyword.
        if l.len() > indent && l.as_bytes()[..indent].iter().all(|&b| b == b' ') {
            text.push(' ');
            text.push_str(l.trim());
            i += 1;
        } else {
            break;
        }
    }
    (text, i)
}

/// Drops a leading keyword from a collected field.
fn strip_field(keyword: &str, text: &str) -> String {
    text.trim_start()
        .strip_prefix(keyword)
        .unwrap_or(text)
        .trim()
        .to_string()
}

/// Parses one `REFERENCE` block starting at `start`. Returns the
/// reference and the line index just past it.
fn parse_reference(lines: &[&str], start: usize) -> (Reference, usize) {
    // The REFERENCE line itself: `REFERENCE   1  (bases 1 to 1859)`
    let header = lines[start];
    let body = strip_field("REFERENCE", header);
    // Pick the first integer for the reference number; everything
    // inside parens (bases ... to ...) goes to `bases`.
    let mut number: usize = 0;
    let mut bases = String::new();
    if let Some(open) = body.find('(') {
        let head = &body[..open];
        for tok in head.split_whitespace() {
            if let Ok(n) = tok.parse::<usize>() {
                number = n;
                break;
            }
        }
        if let Some(close) = body[open + 1..].find(')') {
            let inner = body[open + 1..open + 1 + close].trim();
            // `bases 1 to 1859` -> `1..1859`.
            bases = canonicalize_bases(inner);
        }
    } else {
        for tok in body.split_whitespace() {
            if let Ok(n) = tok.parse::<usize>() {
                number = n;
                break;
            }
        }
    }
    let mut reference = Reference {
        number,
        bases,
        ..Reference::default()
    };

    // Sub-fields are indented 2 columns: AUTHORS, CONSRTM, TITLE,
    // JOURNAL, PUBMED, MEDLINE, REMARK.
    let mut i = start + 1;
    while i < lines.len() {
        let line = lines[i];
        if line.starts_with("REFERENCE") || (!line.starts_with(' ') && !line.is_empty()) {
            break;
        }
        if line.is_empty() {
            i += 1;
            continue;
        }
        let trimmed = line.trim_start();
        // Recognize sub-field keywords.
        let kw_end = trimmed
            .find([' ', '\t'])
            .unwrap_or(trimmed.len());
        let kw = &trimmed[..kw_end];
        match kw {
            "AUTHORS" | "CONSRTM" | "TITLE" | "JOURNAL" | "PUBMED" | "MEDLINE" | "REMARK" => {
                let (text, next) = collect_continuation(lines, i, 12);
                let value = strip_field(kw, &text);
                match kw {
                    "AUTHORS" => reference.authors = value,
                    "CONSRTM" => reference.consortium = value,
                    "TITLE" => reference.title = value,
                    "JOURNAL" => reference.journal = value,
                    "PUBMED" => reference.pubmed = value,
                    "MEDLINE" => reference.medline = value,
                    "REMARK" => reference.remark = value,
                    _ => unreachable!(),
                }
                i = next;
            }
            _ => {
                // Unknown sub-field — skip to keep going.
                i += 1;
            }
        }
    }
    (reference, i)
}

/// Converts a GenBank `bases X to Y` phrasing into `X..Y`.
fn canonicalize_bases(s: &str) -> String {
    let trimmed = s.trim_start().trim_start_matches("bases ").trim();
    // Replace ` to ` with `..` (allowing `to`).
    if let Some((a, b)) = trimmed.split_once(" to ") {
        format!("{}..{}", a.trim(), b.trim())
    } else {
        trimmed.to_string()
    }
}

/// Parses the FEATURES table starting at line `start`. Returns the
/// features and the index just past the table.
fn parse_features(lines: &[&str], start: usize) -> Result<(Vec<SeqFeature>, usize)> {
    let mut features = Vec::new();
    let mut i = start;
    // Feature key starts at column 5; qualifiers at column 21.
    while i < lines.len() {
        let line = lines[i];
        if line.is_empty() {
            i += 1;
            continue;
        }
        // A new top-level section (ORIGIN, etc.) or `//` ends the table.
        if !line.starts_with(' ') {
            break;
        }
        let trimmed = line.trim_start();
        let indent = line.len() - trimmed.len();
        if indent < 21 && !trimmed.starts_with('/') {
            // This is a feature key + location line.
            let mut parts = trimmed.splitn(2, char::is_whitespace);
            let key = parts.next().unwrap_or("").to_string();
            let mut loc_str = parts.next().unwrap_or("").trim().to_string();
            i += 1;
            // The location may wrap onto continuation lines that are
            // indented to the qualifier column but do not start `/`.
            while i < lines.len() {
                let cont = lines[i];
                let ct = cont.trim_start();
                let cind = cont.len() - ct.len();
                if cind >= 21 && !ct.starts_with('/') && !loc_complete(&loc_str) {
                    loc_str.push_str(ct);
                    i += 1;
                } else {
                    break;
                }
            }
            let location = locstr::parse_location(&loc_str)?;
            let mut feature = SeqFeature::new(key, location);
            // Collect qualifiers.
            while i < lines.len() {
                let q = lines[i];
                if q.is_empty() {
                    i += 1;
                    continue;
                }
                let qt = q.trim_start();
                let qind = q.len() - qt.len();
                if !q.starts_with(' ') {
                    break;
                }
                if qind < 21 && !qt.starts_with('/') {
                    // Next feature.
                    break;
                }
                if qt.starts_with('/') {
                    let (key, value, next) = parse_qualifier(lines, i);
                    feature.qualifiers.insert(key, value);
                    i = next;
                } else {
                    // Stray continuation we didn't expect — skip.
                    i += 1;
                }
            }
            features.push(feature);
        } else {
            i += 1;
        }
    }
    Ok((features, i))
}

/// Heuristic: are the parentheses in a location string balanced (so it
/// is syntactically complete)?
fn loc_complete(s: &str) -> bool {
    let opens = s.matches('(').count();
    let closes = s.matches(')').count();
    opens == closes
}

/// Parses one `/key=value` qualifier, collecting continuation lines.
/// Returns `(key, value, next_index)`.
///
/// Continuation handling: a quoted value runs until a closing `"`,
/// which may be many lines later. An unquoted value (`/foo=123` or a
/// bare flag `/pseudo`) is one line. The parser only treats a `/` at
/// the qualifier column (21) as the start of a new qualifier; embedded
/// `/` and `=` inside a quoted value are preserved verbatim.
fn parse_qualifier(lines: &[&str], start: usize) -> (String, String, usize) {
    let first = lines[start].trim_start();
    let body = first.strip_prefix('/').unwrap_or(first);
    let (key, mut value, quoted) = match body.split_once('=') {
        Some((k, v)) => {
            let v = v.trim();
            if let Some(stripped) = v.strip_prefix('"') {
                (k.to_string(), stripped.to_string(), true)
            } else {
                (k.to_string(), v.to_string(), false)
            }
        }
        // A bare flag qualifier like `/pseudo`.
        None => (body.to_string(), String::new(), false),
    };
    let mut i = start + 1;
    if quoted {
        // Continuation lines until a closing quote.
        if let Some(stripped) = value.strip_suffix('"') {
            // Single-line quoted value.
            value = stripped.to_string();
        } else {
            while i < lines.len() {
                let l = lines[i];
                let lt = l.trim_start();
                let lind = l.len() - lt.len();
                if lind < 21 {
                    // Out of the qualifier column — end of the feature.
                    break;
                }
                // A `/` at the qualifier column with no preceding open
                // quote belongs to the NEXT qualifier — but only if we
                // are not in the middle of a quoted value (we are; we
                // got here precisely because the open `"` was not
                // closed on the prior line). So keep going.
                value.push(' ');
                if let Some(end) = lt.strip_suffix('"') {
                    value.push_str(end);
                    i += 1;
                    break;
                } else {
                    value.push_str(lt);
                    i += 1;
                }
            }
        }
    }
    (key, value, i)
}

/// Parses the ORIGIN sequence block. Returns the concatenated residues
/// and the index just past `//`.
fn parse_origin(lines: &[&str], start: usize) -> (String, usize) {
    let mut seq = String::new();
    let mut i = start;
    while i < lines.len() {
        let line = lines[i];
        if line.starts_with("//") {
            i += 1;
            break;
        }
        for tok in line.split_whitespace() {
            // Skip the leading base-count number on each ORIGIN line.
            if tok.chars().all(|c| c.is_ascii_digit()) {
                continue;
            }
            seq.push_str(tok);
        }
        i += 1;
    }
    (seq.to_ascii_uppercase(), i)
}

/// Serializes a [`SeqRecord`] to a GenBank flat file.
///
/// Emits `LOCUS` / `DEFINITION` / `ACCESSION` / `VERSION` / `KEYWORDS`
/// / `SOURCE` + `ORGANISM` / one `REFERENCE` per [`Reference`] /
/// `FEATURES` / `ORIGIN`. The `LOCUS` molecule type and topology are
/// taken from the record's [`Seq`].
pub fn write(rec: &SeqRecord) -> String {
    let mut out = String::new();
    let mol = match rec.seq.kind() {
        SeqKind::Dna => "DNA",
        SeqKind::Rna => "RNA",
        SeqKind::Protein => "AA",
    };
    let topo = if rec.seq.is_circular() {
        "circular"
    } else {
        "linear"
    };
    let locus_name = if rec.name.is_empty() {
        &rec.id
    } else {
        &rec.name
    };
    out.push_str(&format!(
        "LOCUS       {:<16} {} bp    {} {} UNK\n",
        locus_name,
        rec.seq.len(),
        mol,
        topo
    ));
    let definition = if rec.description.is_empty() {
        "."
    } else {
        &rec.description
    };
    out.push_str(&format!("DEFINITION  {definition}\n"));
    if let Some(acc) = rec.annotations.get("accession") {
        out.push_str(&format!("ACCESSION   {acc}\n"));
    } else {
        out.push_str(&format!("ACCESSION   {}\n", rec.id));
    }
    if let Some(ver) = rec.annotations.get("version") {
        out.push_str(&format!("VERSION     {ver}\n"));
    }
    if let Some(kw) = rec.annotations.get("keywords") {
        out.push_str(&format!("KEYWORDS    {kw}.\n"));
    }
    if let Some(org) = rec.annotations.get("organism") {
        out.push_str(&format!("SOURCE      {org}\n"));
        out.push_str(&format!("  ORGANISM  {org}\n"));
    }
    for reference in &rec.references {
        write_reference(&mut out, reference);
    }
    // FEATURES table.
    out.push_str("FEATURES             Location/Qualifiers\n");
    for f in &rec.features {
        let loc = locstr::write_location(&f.location);
        out.push_str(&format!("     {:<15} {}\n", f.feature_type, loc));
        for (k, v) in &f.qualifiers {
            if v.is_empty() {
                out.push_str(&format!("                     /{k}\n"));
            } else if is_unquoted_qualifier(k) {
                out.push_str(&format!("                     /{k}={v}\n"));
            } else {
                out.push_str(&format!("                     /{k}=\"{v}\"\n"));
            }
        }
    }
    // ORIGIN — 60 bases per line, 6 blocks of 10, with a 1-based count.
    out.push_str("ORIGIN\n");
    let lower = rec.seq.as_str().to_ascii_lowercase();
    let lower_bytes = lower.as_bytes();
    let mut pos = 0;
    while pos < lower_bytes.len() {
        let end = (pos + 60).min(lower_bytes.len());
        out.push_str(&format!("{:>9}", pos + 1));
        let mut col = pos;
        while col < end {
            let block_end = (col + 10).min(end);
            out.push(' ');
            out.push_str(std::str::from_utf8(&lower_bytes[col..block_end]).unwrap());
            col = block_end;
        }
        out.push('\n');
        pos = end;
    }
    out.push_str("//\n");
    out
}

/// Writes one [`Reference`] block in canonical GenBank form.
fn write_reference(out: &mut String, r: &Reference) {
    if r.bases.is_empty() {
        out.push_str(&format!("REFERENCE   {}\n", r.number));
    } else {
        // Emit `1  (bases 1 to 1859)` if the bases field looks like a
        // range; otherwise quote it verbatim.
        let bases_phrase = if let Some((a, b)) = r.bases.split_once("..") {
            format!("(bases {a} to {b})")
        } else {
            format!("({})", r.bases)
        };
        out.push_str(&format!(
            "REFERENCE   {}  {}\n",
            r.number, bases_phrase
        ));
    }
    if !r.authors.is_empty() {
        out.push_str(&format!("  AUTHORS   {}\n", r.authors));
    }
    if !r.consortium.is_empty() {
        out.push_str(&format!("  CONSRTM   {}\n", r.consortium));
    }
    if !r.title.is_empty() {
        out.push_str(&format!("  TITLE     {}\n", r.title));
    }
    if !r.journal.is_empty() {
        out.push_str(&format!("  JOURNAL   {}\n", r.journal));
    }
    if !r.pubmed.is_empty() {
        out.push_str(&format!("  PUBMED    {}\n", r.pubmed));
    }
    if !r.medline.is_empty() {
        out.push_str(&format!("  MEDLINE   {}\n", r.medline));
    }
    if !r.remark.is_empty() {
        out.push_str(&format!("  REMARK    {}\n", r.remark));
    }
}

/// The INSDC feature-table spec marks a small set of qualifiers as
/// taking an unquoted value (`/codon_start=1`, `/transl_table=11`,
/// etc.). Everything else is double-quoted.
fn is_unquoted_qualifier(key: &str) -> bool {
    matches!(
        key,
        "codon_start"
            | "transl_table"
            | "translation"
            | "anticodon"
            | "estimated_length"
            | "number"
            | "rpt_unit_range"
            | "tag_peptide"
            | "compare"
            | "citation"
            | "direction"
            | "rpt_type"
            | "transl_except"
            | "mod_base"
    )
    .then_some(true)
    .map(|b| {
        // The list above is conservative; many INSDC viewers accept the
        // value quoted, so we ONLY emit unquoted for the small handful
        // where the spec forbids the quotes outright.
        match key {
            "codon_start" | "transl_table" | "number" | "estimated_length" | "citation"
            | "direction" => b,
            _ => false,
        }
    })
    .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Location;

    // NOTE: a raw string literal — GenBank is column-sensitive (feature
    // keys at column 5, qualifiers at column 21), so the leading
    // whitespace of every line is load-bearing and must be preserved.
    const SAMPLE: &str = r#"LOCUS       TESTSEQ                30 bp    DNA     linear   UNK
DEFINITION  A small synthetic test sequence.
ACCESSION   TEST001
VERSION     TEST001.1
SOURCE      synthetic construct
  ORGANISM  synthetic construct
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="synthetic construct"
     CDS             1..30
                     /gene="testGene"
                     /product="test protein"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;

    #[test]
    fn parse_basic_record() {
        let rec = parse(SAMPLE).unwrap();
        assert_eq!(rec.id, "TEST001");
        assert_eq!(rec.name, "TESTSEQ");
        assert_eq!(rec.description, "A small synthetic test sequence.");
        assert_eq!(rec.seq.len(), 30);
        assert_eq!(rec.seq.as_str(), "ATGAAAGGGTTTCCCGGGAAATTTTAGCTA");
        assert_eq!(rec.seq.kind(), SeqKind::Dna);
    }

    #[test]
    fn parse_features_and_qualifiers() {
        let rec = parse(SAMPLE).unwrap();
        assert_eq!(rec.features.len(), 2);
        let cds = rec
            .features
            .iter()
            .find(|f| f.feature_type == "CDS")
            .unwrap();
        assert_eq!(cds.qualifier("gene"), Some("testGene"));
        assert_eq!(cds.qualifier("product"), Some("test protein"));
        assert_eq!(cds.location, Location::single(0, 30));
    }

    #[test]
    fn parse_annotations() {
        let rec = parse(SAMPLE).unwrap();
        assert_eq!(
            rec.annotations.get("organism").map(String::as_str),
            Some("synthetic construct")
        );
        assert_eq!(
            rec.annotations.get("version").map(String::as_str),
            Some("TEST001.1")
        );
    }

    #[test]
    fn circular_topology_from_locus() {
        let gb = "LOCUS       PLASMID    10 bp    DNA     circular UNK\n\
DEFINITION  .\n\
ORIGIN\n        1 atgcatgcat\n//\n";
        let rec = parse(gb).unwrap();
        assert!(rec.seq.is_circular());
    }

    #[test]
    fn missing_origin_is_error() {
        let gb = "LOCUS       X    5 bp    DNA     linear UNK\nDEFINITION  .\n//\n";
        assert!(parse(gb).is_err());
    }

    #[test]
    fn complement_join_location_parses() {
        // Raw string — the FEATURES block is column-sensitive.
        let gb = r#"LOCUS       X    30 bp    DNA     linear UNK
DEFINITION  .
FEATURES             Location/Qualifiers
     CDS             complement(join(1..6,20..24))
                     /gene="spliced"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;
        let rec = parse(gb).unwrap();
        let cds = &rec.features[0];
        assert!(cds.location.is_reverse());
        assert_eq!(cds.location.spans().len(), 2);
    }

    #[test]
    fn write_then_parse_roundtrip() {
        let rec = parse(SAMPLE).unwrap();
        let text = write(&rec);
        let reparsed = parse(&text).unwrap();
        assert_eq!(rec.seq.as_str(), reparsed.seq.as_str());
        assert_eq!(rec.description, reparsed.description);
        assert_eq!(rec.features.len(), reparsed.features.len());
        assert_eq!(
            rec.features[1].qualifier("gene"),
            reparsed.features[1].qualifier("gene")
        );
    }

    #[test]
    fn written_file_has_required_sections() {
        let rec = parse(SAMPLE).unwrap();
        let text = write(&rec);
        assert!(text.starts_with("LOCUS"));
        assert!(text.contains("DEFINITION"));
        assert!(text.contains("FEATURES"));
        assert!(text.contains("ORIGIN"));
        assert!(text.trim_end().ends_with("//"));
    }

    // -------------------------------------------------------------------
    // Reference-block coverage.
    // -------------------------------------------------------------------

    const SAMPLE_REF: &str = r#"LOCUS       REFSAMPLE              30 bp    DNA     linear   UNK
DEFINITION  Reference block fixture.
ACCESSION   REF001
VERSION     REF001.1
KEYWORDS    synthetic; test.
SOURCE      synthetic construct
  ORGANISM  synthetic construct
REFERENCE   1  (bases 1 to 30)
  AUTHORS   Doe,J. and Roe,R.
  TITLE     A synthetic test article
  JOURNAL   J. Test. Biol. 1, 1-10 (2024)
  PUBMED    12345678
REFERENCE   2  (bases 1 to 15)
  CONSRTM   Synthetic Sequencing Consortium
  TITLE     Second reference for the partial range
  JOURNAL   Unpublished
  REMARK    Submitted by the consortium.
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="synthetic construct"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;

    #[test]
    fn reference_blocks_parse() {
        let rec = parse(SAMPLE_REF).unwrap();
        assert_eq!(rec.references.len(), 2);
        let r0 = &rec.references[0];
        assert_eq!(r0.number, 1);
        assert_eq!(r0.bases, "1..30");
        assert_eq!(r0.authors, "Doe,J. and Roe,R.");
        assert_eq!(r0.title, "A synthetic test article");
        assert_eq!(r0.journal, "J. Test. Biol. 1, 1-10 (2024)");
        assert_eq!(r0.pubmed, "12345678");
        let r1 = &rec.references[1];
        assert_eq!(r1.number, 2);
        assert_eq!(r1.bases, "1..15");
        assert_eq!(r1.consortium, "Synthetic Sequencing Consortium");
        assert_eq!(r1.remark, "Submitted by the consortium.");
    }

    #[test]
    fn keywords_field_parses() {
        let rec = parse(SAMPLE_REF).unwrap();
        assert_eq!(
            rec.annotations.get("keywords").map(String::as_str),
            Some("synthetic; test")
        );
    }

    #[test]
    fn reference_blocks_round_trip() {
        let rec = parse(SAMPLE_REF).unwrap();
        let text = write(&rec);
        let reparsed = parse(&text).unwrap();
        assert_eq!(reparsed.references.len(), 2);
        assert_eq!(reparsed.references, rec.references);
    }

    // -------------------------------------------------------------------
    // Full qualifier handling: multi-line, embedded `=`/`/`.
    // -------------------------------------------------------------------

    #[test]
    fn qualifier_with_embedded_slash_and_equals_in_quoted_value_round_trips() {
        // Most GenBank /note= qualifiers contain free text; the parser
        // must not split on `=` or `/` characters that are inside a
        // quoted value.
        let gb = r#"LOCUS       QUALTST                30 bp    DNA     linear   UNK
DEFINITION  .
ACCESSION   QT01
FEATURES             Location/Qualifiers
     misc_feature    1..30
                     /note="a/b=c and a long second
                     line of text including = and / characters"
                     /label="x/y/z"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;
        let rec = parse(gb).unwrap();
        let f = &rec.features[0];
        let note = f.qualifier("note").unwrap();
        assert!(note.starts_with("a/b=c"), "got: {note}");
        assert!(note.contains("characters"), "got: {note}");
        // Embedded `=`/`/` survive intact.
        assert!(note.contains('='));
        assert!(note.contains('/'));
        // The second qualifier is not lost.
        assert_eq!(f.qualifier("label"), Some("x/y/z"));
    }

    #[test]
    fn order_location_parses_and_round_trips() {
        let gb = r#"LOCUS       ORDTST                 30 bp    DNA     linear   UNK
DEFINITION  .
ACCESSION   OT01
FEATURES             Location/Qualifiers
     primer_bind     order(1..6,20..24)
                     /label="paired_primer"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;
        let rec = parse(gb).unwrap();
        assert!(matches!(rec.features[0].location, Location::Order(_)));
        // Round-trip preserves the `order` operator (no silent
        // collapse to `join`).
        let text = write(&rec);
        assert!(text.contains("order(1..6,20..24)"));
    }

    #[test]
    fn between_bases_location_parses_and_round_trips() {
        let gb = r#"LOCUS       BTWTST                 30 bp    DNA     linear   UNK
DEFINITION  .
ACCESSION   BT01
FEATURES             Location/Qualifiers
     misc_recomb     10^11
                     /note="recombination point"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;
        let rec = parse(gb).unwrap();
        assert!(matches!(rec.features[0].location, Location::Between { .. }));
        let text = write(&rec);
        assert!(text.contains("10^11"));
    }

    #[test]
    fn cross_record_location_raises_typed_error() {
        let gb = r#"LOCUS       XRTST                  30 bp    DNA     linear   UNK
DEFINITION  .
ACCESSION   XR01
FEATURES             Location/Qualifiers
     misc_feature    J00194.1:1..10
                     /note="cross-record"
ORIGIN
        1 atgaaagggt ttcccgggaa attttagcta
//
"#;
        let err = parse(gb).unwrap_err();
        match err {
            BioseqError::CrossRecordLocation { accession, .. } => {
                assert_eq!(accession, "J00194.1");
            }
            other => panic!("expected CrossRecordLocation, got {other:?}"),
        }
    }
}
