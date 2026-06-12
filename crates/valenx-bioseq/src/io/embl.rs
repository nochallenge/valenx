//! EMBL flat-file reader and writer.
//!
//! The EMBL format is the European counterpart of GenBank. Each line
//! begins with a two-letter line-type code:
//!
//! - `ID` — identification (id, topology, molecule type, length)
//! - `AC` — accession number(s)
//! - `DE` — description
//! - `KW` — keywords
//! - `OS` — organism / source
//! - `OC` — organism classification (skipped)
//! - `RN`/`RC`/`RP`/`RA`/`RT`/`RL`/`RX` — reference block
//! - `FT` — feature table (key + location + qualifiers)
//! - `SQ` — sequence header, followed by sequence lines
//! - `//` — end of record
//!
//! Reference blocks are parsed into [`Reference`] entries and emitted
//! by [`write()`]. The parser preserves multi-line `RA`/`RT`/`RL` text
//! verbatim and recognizes `RX PUBMED;...` / `RX MEDLINE;...`
//! cross-references.

use crate::alphabet::SeqKind;
use crate::error::{BioseqError, Result};
use crate::io::locstr;
use crate::record::{Reference, SeqFeature, SeqRecord};
use crate::seq::{Seq, Topology};
use std::collections::BTreeMap;

/// Parses an EMBL flat file holding a single record.
///
/// Returns [`BioseqError::Parse`] if no `SQ` sequence block is found.
pub fn parse(text: &str) -> Result<SeqRecord> {
    let lines: Vec<&str> = text.lines().collect();
    let mut id = String::new();
    let mut accession = String::new();
    let mut description = String::new();
    let mut keywords = String::new();
    let mut kind = SeqKind::Dna;
    let mut topology = Topology::Linear;
    let mut annotations: BTreeMap<String, String> = BTreeMap::new();
    let mut features: Vec<SeqFeature> = Vec::new();
    let mut references: Vec<Reference> = Vec::new();
    let mut current_ref: Option<Reference> = None;
    let mut origin = String::new();
    let mut in_sequence = false;

    let mut i = 0;
    while i < lines.len() {
        let line = lines[i];
        if line.starts_with("//") {
            break;
        }
        // The two-letter line code lives in columns 1-2. `line.get(..2)`
        // is char-boundary-safe: it returns `None` when fewer than two
        // bytes are present *or* when byte index 2 falls inside a
        // multi-byte UTF-8 char (a malformed line — real EMBL codes are
        // always two ASCII letters). Either way the line is not a valid
        // record line, so skip it rather than panicking on a raw byte
        // slice.
        let code = match line.get(..2) {
            Some(c) => c,
            None => {
                i += 1;
                continue;
            }
        };
        let body = line.get(5..).unwrap_or("").trim_end();

        match code {
            "ID" => {
                let (parsed_id, k, t) = parse_id(body);
                id = parsed_id;
                kind = k;
                topology = t;
                i += 1;
            }
            "AC" => {
                // Accession line: semicolon-separated, first is primary.
                if accession.is_empty() {
                    if let Some(first) = body.split(';').next() {
                        accession = first.trim().to_string();
                    }
                }
                i += 1;
            }
            "DE" => {
                if !description.is_empty() {
                    description.push(' ');
                }
                description.push_str(body.trim());
                i += 1;
            }
            "KW" => {
                if !keywords.is_empty() {
                    keywords.push(' ');
                }
                keywords.push_str(body.trim_end_matches('.').trim());
                i += 1;
            }
            "OS" => {
                annotations
                    .entry("organism".to_string())
                    .or_insert_with(|| body.trim().to_string());
                i += 1;
            }
            "RN" => {
                // Push any prior reference and start a new one.
                if let Some(r) = current_ref.take() {
                    references.push(r);
                }
                // `RN   [1]` -> number = 1
                let mut r = Reference::default();
                let trimmed = body.trim().trim_start_matches('[').trim_end_matches(']');
                if let Ok(n) = trimmed.parse::<usize>() {
                    r.number = n;
                }
                current_ref = Some(r);
                i += 1;
            }
            "RC" => {
                if let Some(r) = current_ref.as_mut() {
                    if !r.remark.is_empty() {
                        r.remark.push(' ');
                    }
                    r.remark.push_str(body.trim());
                }
                i += 1;
            }
            "RP" => {
                if let Some(r) = current_ref.as_mut() {
                    if !r.bases.is_empty() {
                        r.bases.push(' ');
                    }
                    r.bases.push_str(body.trim().trim_end_matches('.'));
                }
                i += 1;
            }
            "RA" => {
                if let Some(r) = current_ref.as_mut() {
                    if !r.authors.is_empty() {
                        r.authors.push(' ');
                    }
                    r.authors.push_str(body.trim().trim_end_matches(';'));
                }
                i += 1;
            }
            "RG" => {
                if let Some(r) = current_ref.as_mut() {
                    if !r.consortium.is_empty() {
                        r.consortium.push(' ');
                    }
                    r.consortium.push_str(body.trim().trim_end_matches(';'));
                }
                i += 1;
            }
            "RT" => {
                if let Some(r) = current_ref.as_mut() {
                    // EMBL RT format: `RT   "title text";` — strip the
                    // trailing `;`, then strip surrounding quotes.
                    let mut text = body.trim().to_string();
                    if let Some(stripped) = text.strip_suffix(';') {
                        text = stripped.to_string();
                    }
                    text = text.trim().to_string();
                    if let Some(rest) = text.strip_prefix('"') {
                        text = rest.to_string();
                    }
                    if let Some(rest) = text.strip_suffix('"') {
                        text = rest.to_string();
                    }
                    let text = text.trim();
                    if !text.is_empty() {
                        if !r.title.is_empty() {
                            r.title.push(' ');
                        }
                        r.title.push_str(text);
                    }
                }
                i += 1;
            }
            "RL" => {
                if let Some(r) = current_ref.as_mut() {
                    if !r.journal.is_empty() {
                        r.journal.push(' ');
                    }
                    r.journal.push_str(body.trim().trim_end_matches('.'));
                }
                i += 1;
            }
            "RX" => {
                if let Some(r) = current_ref.as_mut() {
                    let trimmed = body.trim().trim_end_matches('.');
                    if let Some(rest) = trimmed.strip_prefix("PUBMED;") {
                        r.pubmed = rest.trim().to_string();
                    } else if let Some(rest) = trimmed.strip_prefix("MEDLINE;") {
                        r.medline = rest.trim().to_string();
                    } else if let Some(rest) = trimmed.strip_prefix("DOI;") {
                        // Stash a DOI in the journal field if no
                        // journal is set, otherwise as a remark; this
                        // is a pragmatic choice since `Reference`
                        // does not have a dedicated DOI slot.
                        let doi = rest.trim();
                        if r.journal.is_empty() {
                            r.journal = format!("DOI:{doi}");
                        } else if !r.remark.contains(doi) {
                            if !r.remark.is_empty() {
                                r.remark.push(' ');
                            }
                            r.remark.push_str(&format!("DOI:{doi}"));
                        }
                    }
                }
                i += 1;
            }
            "FT" => {
                let (parsed, next) = parse_feature_table(&lines, i)?;
                features = parsed;
                i = next;
            }
            "SQ" => {
                in_sequence = true;
                i += 1;
            }
            "  " if in_sequence => {
                // Sequence data line: residues with a trailing count.
                for tok in line.split_whitespace() {
                    if tok.chars().all(|c| c.is_ascii_digit()) {
                        continue;
                    }
                    origin.push_str(tok);
                }
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }
    // Flush the trailing reference, if any.
    if let Some(r) = current_ref.take() {
        references.push(r);
    }

    if origin.is_empty() {
        return Err(BioseqError::parse("embl", "no SQ sequence block found"));
    }
    let final_id = if !accession.is_empty() {
        accession.clone()
    } else {
        id.clone()
    };
    if !accession.is_empty() {
        annotations.insert("accession".to_string(), accession);
    }
    if !keywords.is_empty() {
        annotations.insert("keywords".to_string(), keywords);
    }
    let seq = Seq::with_topology(kind, origin.to_ascii_uppercase(), topology)?;
    Ok(SeqRecord {
        id: final_id,
        name: id,
        description,
        seq,
        features,
        annotations,
        references,
    })
}

/// Serializes a [`SeqRecord`] to an EMBL flat file.
///
/// Emits `ID` / `AC` / `DE` / `KW` / `OS` / `RN/RA/RT/RL/RX` per
/// reference / `FT` / `SQ`. The `ID` molecule type and topology are
/// taken from the record's [`Seq`].
pub fn write(rec: &SeqRecord) -> String {
    let mut out = String::new();
    let mol = match rec.seq.kind() {
        SeqKind::Dna => "DNA",
        SeqKind::Rna => "mRNA",
        SeqKind::Protein => "AA",
    };
    let topo = if rec.seq.is_circular() {
        "circular"
    } else {
        "linear"
    };
    let id_name = if rec.name.is_empty() {
        &rec.id
    } else {
        &rec.name
    };
    out.push_str(&format!(
        "ID   {}; SV 1; {}; {}; STD; UNC; {} BP.\nXX\n",
        id_name,
        topo,
        mol,
        rec.seq.len()
    ));
    if let Some(acc) = rec.annotations.get("accession") {
        out.push_str(&format!("AC   {acc};\nXX\n"));
    }
    if !rec.description.is_empty() {
        out.push_str(&format!("DE   {}\nXX\n", rec.description));
    }
    if let Some(kw) = rec.annotations.get("keywords") {
        out.push_str(&format!("KW   {kw}.\nXX\n"));
    }
    if let Some(org) = rec.annotations.get("organism") {
        out.push_str(&format!("OS   {org}\nXX\n"));
    }
    for r in &rec.references {
        out.push_str(&format!("RN   [{}]\n", r.number));
        if !r.bases.is_empty() {
            out.push_str(&format!("RP   {}\n", r.bases));
        }
        if !r.remark.is_empty() {
            out.push_str(&format!("RC   {}\n", r.remark));
        }
        if !r.authors.is_empty() {
            out.push_str(&format!("RA   {};\n", r.authors));
        }
        if !r.consortium.is_empty() {
            out.push_str(&format!("RG   {};\n", r.consortium));
        }
        if !r.title.is_empty() {
            out.push_str(&format!("RT   \"{}\";\n", r.title));
        }
        if !r.journal.is_empty() {
            out.push_str(&format!("RL   {}.\n", r.journal));
        }
        if !r.pubmed.is_empty() {
            out.push_str(&format!("RX   PUBMED; {}.\n", r.pubmed));
        }
        if !r.medline.is_empty() {
            out.push_str(&format!("RX   MEDLINE; {}.\n", r.medline));
        }
        out.push_str("XX\n");
    }
    // FT feature table.
    out.push_str("FH   Key             Location/Qualifiers\n");
    out.push_str("FH\n");
    for f in &rec.features {
        let loc = locstr::write_location(&f.location);
        out.push_str(&format!("FT   {:<16}{}\n", f.feature_type, loc));
        for (k, v) in &f.qualifiers {
            if v.is_empty() {
                out.push_str(&format!("FT                   /{k}\n"));
            } else {
                out.push_str(&format!("FT                   /{k}=\"{v}\"\n"));
            }
        }
    }
    out.push_str("XX\n");
    // SQ sequence block — 60 bp per line, 6 blocks of 10, with a
    // right-aligned 1-based count.
    let lower = rec.seq.as_str().to_ascii_lowercase();
    let lower_bytes = lower.as_bytes();
    let (na, nc, ng, nt, no) = nucleotide_counts(lower_bytes);
    out.push_str(&format!(
        "SQ   Sequence {} BP; {} A; {} C; {} G; {} T; {} other;\n",
        rec.seq.len(),
        na,
        nc,
        ng,
        nt,
        no
    ));
    let mut pos = 0;
    while pos < lower_bytes.len() {
        let end = (pos + 60).min(lower_bytes.len());
        out.push_str("     ");
        let mut col = pos;
        while col < end {
            let block_end = (col + 10).min(end);
            out.push_str(std::str::from_utf8(&lower_bytes[col..block_end]).unwrap());
            out.push(' ');
            col = block_end;
        }
        out.push_str(&format!("{end:>10}\n"));
        pos = end;
    }
    out.push_str("//\n");
    out
}

/// Counts `A C G T` and `other` (everything else) in a residue slice.
fn nucleotide_counts(bytes: &[u8]) -> (usize, usize, usize, usize, usize) {
    let mut a = 0;
    let mut c = 0;
    let mut g = 0;
    let mut t = 0;
    let mut o = 0;
    for &b in bytes {
        match b.to_ascii_uppercase() {
            b'A' => a += 1,
            b'C' => c += 1,
            b'G' => g += 1,
            b'T' | b'U' => t += 1,
            _ => o += 1,
        }
    }
    (a, c, g, t, o)
}

/// Parses the `ID` line body → `(id, kind, topology)`.
///
/// Modern EMBL `ID` lines are semicolon-delimited:
/// `ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.`
fn parse_id(body: &str) -> (String, SeqKind, Topology) {
    let mut kind = SeqKind::Dna;
    let mut topology = Topology::Linear;
    let parts: Vec<&str> = body.split(';').map(str::trim).collect();
    let id = parts
        .first()
        .map(|s| s.split_whitespace().next().unwrap_or("").to_string())
        .unwrap_or_default();
    for p in &parts {
        let lower = p.to_ascii_lowercase();
        if lower == "circular" {
            topology = Topology::Circular;
        } else if lower == "linear" {
            topology = Topology::Linear;
        } else if lower.contains("rna") {
            kind = SeqKind::Rna;
        } else if lower.contains("protein") || lower == "aa" {
            kind = SeqKind::Protein;
        } else if lower.contains("dna") {
            kind = SeqKind::Dna;
        }
    }
    (id, kind, topology)
}

/// Parses the `FT` feature-table lines starting at index `start`.
fn parse_feature_table(lines: &[&str], start: usize) -> Result<(Vec<SeqFeature>, usize)> {
    let mut features: Vec<SeqFeature> = Vec::new();
    let mut i = start;
    // Pending feature being assembled.
    let mut cur_key: Option<String> = None;
    let mut cur_loc = String::new();
    let mut cur_quals: BTreeMap<String, String> = BTreeMap::new();
    let mut pending_qual: Option<(String, String, bool)> = None;

    let flush = |key: &mut Option<String>,
                 loc: &mut String,
                 quals: &mut BTreeMap<String, String>,
                 out: &mut Vec<SeqFeature>|
     -> Result<()> {
        if let Some(k) = key.take() {
            let location = locstr::parse_location(loc.trim())?;
            let mut f = SeqFeature::new(k, location);
            f.qualifiers = std::mem::take(quals);
            out.push(f);
        }
        loc.clear();
        Ok(())
    };

    while i < lines.len() {
        let line = lines[i];
        if !line.starts_with("FT") {
            break;
        }
        let body = line.get(5..).unwrap_or("");
        // Feature-key lines have a non-space at column 5 (index 5);
        // continuation/qualifier lines are indented further.
        let body_trimmed = body.trim_start();
        let key_line = !body.is_empty() && !body.starts_with(' ') && !body_trimmed.starts_with('/');

        if key_line {
            // Close the qualifier and feature in progress.
            if let Some((k, v, _)) = pending_qual.take() {
                cur_quals.insert(k, v);
            }
            flush(&mut cur_key, &mut cur_loc, &mut cur_quals, &mut features)?;
            let mut sp = body_trimmed.splitn(2, char::is_whitespace);
            cur_key = Some(sp.next().unwrap_or("").to_string());
            cur_loc = sp.next().unwrap_or("").trim().to_string();
        } else if body_trimmed.starts_with('/') {
            // New qualifier; close any pending one.
            if let Some((k, v, _)) = pending_qual.take() {
                cur_quals.insert(k, v);
            }
            let q = body_trimmed.strip_prefix('/').unwrap_or(body_trimmed);
            match q.split_once('=') {
                Some((k, v)) => {
                    let v = v.trim();
                    if let Some(rest) = v.strip_prefix('"') {
                        if let Some(done) = rest.strip_suffix('"') {
                            cur_quals.insert(k.to_string(), done.to_string());
                        } else {
                            pending_qual = Some((k.to_string(), rest.to_string(), true));
                        }
                    } else {
                        cur_quals.insert(k.to_string(), v.to_string());
                    }
                }
                None => {
                    cur_quals.insert(q.to_string(), String::new());
                }
            }
        } else {
            // Continuation: of a quoted qualifier, or of the location.
            let cont = body_trimmed;
            if let Some((k, mut v, _)) = pending_qual.take() {
                v.push(' ');
                if let Some(done) = cont.strip_suffix('"') {
                    v.push_str(done);
                    cur_quals.insert(k, v);
                } else {
                    v.push_str(cont);
                    pending_qual = Some((k, v, true));
                }
            } else if cur_key.is_some() {
                // Location continuation.
                cur_loc.push_str(cont);
            }
        }
        i += 1;
    }
    if let Some((k, v, _)) = pending_qual.take() {
        cur_quals.insert(k, v);
    }
    flush(&mut cur_key, &mut cur_loc, &mut cur_quals, &mut features)?;
    Ok((features, i))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Location;

    // NOTE: a raw string literal — EMBL is column-sensitive (the two-
    // letter line code in columns 1-2, sequence lines indented with
    // spaces), so the leading whitespace of every line is load-bearing.
    const SAMPLE: &str = r#"ID   TEST01; SV 1; linear; DNA; STD; SYN; 30 BP.
AC   TEST01;
DE   A small synthetic EMBL test sequence.
OS   synthetic construct
FT   source          1..30
FT                   /organism="synthetic construct"
FT   CDS             1..30
FT                   /gene="emblGene"
FT                   /product="embl test protein"
SQ   Sequence 30 BP; 9 A; 6 C; 9 G; 6 T; 0 other;
     atgaaagggt ttcccgggaa attttagcta                                  30
//
"#;

    #[test]
    fn parse_basic_record() {
        let rec = parse(SAMPLE).unwrap();
        assert_eq!(rec.id, "TEST01");
        assert_eq!(rec.description, "A small synthetic EMBL test sequence.");
        assert_eq!(rec.seq.len(), 30);
        assert_eq!(rec.seq.as_str(), "ATGAAAGGGTTTCCCGGGAAATTTTAGCTA");
    }

    #[test]
    fn parse_features() {
        let rec = parse(SAMPLE).unwrap();
        assert_eq!(rec.features.len(), 2);
        let cds = rec
            .features
            .iter()
            .find(|f| f.feature_type == "CDS")
            .unwrap();
        assert_eq!(cds.qualifier("gene"), Some("emblGene"));
        assert_eq!(cds.qualifier("product"), Some("embl test protein"));
        assert_eq!(cds.location, Location::single(0, 30));
    }

    #[test]
    fn parse_organism_annotation() {
        let rec = parse(SAMPLE).unwrap();
        assert_eq!(
            rec.annotations.get("organism").map(String::as_str),
            Some("synthetic construct")
        );
    }

    #[test]
    fn malformed_multibyte_line_code_does_not_panic() {
        // A line whose second byte is *inside* a multi-byte UTF-8 char.
        // `a` is one byte, `€` is three, so byte index 2 lands mid-char.
        // The old `&line[..2]` sliced on a raw byte length and panicked
        // ("byte index 2 is not a char boundary"). A real EMBL line code
        // is always two ASCII chars, so such a line is malformed and must
        // be skipped gracefully — never abort the parser. The record
        // still parses on the strength of its valid lines.
        let embl = "a\u{20AC}foo bar\nID   X; SV 1; linear; DNA; STD; SYN; 4 BP.\nSQ   Sequence 4 BP;\n     acgt                                                              4\n//\n";
        let rec = parse(embl).expect("malformed leading line must be skipped, not panic");
        assert_eq!(rec.id, "X");
        assert_eq!(rec.seq.len(), 4);
    }

    #[test]
    fn circular_topology() {
        let embl = r#"ID   P1; SV 1; circular; DNA; STD; SYN; 10 BP.
SQ   Sequence 10 BP;
     atgcatgcat                                                        10
//
"#;
        let rec = parse(embl).unwrap();
        assert!(rec.seq.is_circular());
    }

    #[test]
    fn missing_sequence_is_error() {
        let embl = "ID   X; SV 1; linear; DNA; STD; SYN; 5 BP.\nDE   no seq.\n//\n";
        assert!(parse(embl).is_err());
    }

    #[test]
    fn rna_molecule_type() {
        let embl = r#"ID   R1; SV 1; linear; mRNA; STD; SYN; 6 BP.
SQ   Sequence 6 BP;
     augcau                                                             6
//
"#;
        let rec = parse(embl).unwrap();
        assert_eq!(rec.seq.kind(), SeqKind::Rna);
    }

    // -------------------------------------------------------------------
    // Reference blocks (RN/RC/RP/RA/RT/RL/RX).
    // -------------------------------------------------------------------

    const SAMPLE_REF_EMBL: &str = r#"ID   REF01; SV 1; linear; DNA; STD; SYN; 30 BP.
AC   REF01;
DE   Reference block fixture.
KW   synthetic; test.
OS   synthetic construct
RN   [1]
RP   1-30
RA   Doe,J., Roe,R.;
RT   "A synthetic test article";
RL   J. Test. Biol. 1:1-10(2024).
RX   PUBMED; 12345678.
RN   [2]
RC   second reference for the partial range
RG   Synthetic Sequencing Consortium;
RT   "Partial submission";
RL   Unpublished.
FT   source          1..30
FT                   /organism="synthetic construct"
SQ   Sequence 30 BP; 9 A; 6 C; 9 G; 6 T; 0 other;
     atgaaagggt ttcccgggaa attttagcta                                  30
//
"#;

    #[test]
    fn embl_reference_blocks_parse() {
        let rec = parse(SAMPLE_REF_EMBL).unwrap();
        assert_eq!(rec.references.len(), 2);
        let r0 = &rec.references[0];
        assert_eq!(r0.number, 1);
        assert!(r0.bases.contains("1-30"));
        assert!(r0.authors.contains("Doe"));
        assert_eq!(r0.title, "A synthetic test article");
        assert!(r0.journal.contains("J. Test. Biol."));
        assert_eq!(r0.pubmed, "12345678");

        let r1 = &rec.references[1];
        assert_eq!(r1.number, 2);
        assert!(r1.remark.contains("partial range"));
        assert_eq!(r1.consortium, "Synthetic Sequencing Consortium");
    }

    #[test]
    fn embl_keywords_parse() {
        let rec = parse(SAMPLE_REF_EMBL).unwrap();
        let kw = rec.annotations.get("keywords").unwrap();
        assert!(kw.contains("synthetic"));
        assert!(kw.contains("test"));
    }

    #[test]
    fn embl_writes_a_valid_record() {
        let rec = parse(SAMPLE_REF_EMBL).unwrap();
        let text = write(&rec);
        // Must contain the required header lines.
        assert!(text.starts_with("ID"));
        assert!(text.contains("RN   [1]"));
        assert!(text.contains("RN   [2]"));
        assert!(text.contains("RT   \"A synthetic test article\""));
        assert!(text.contains("RX   PUBMED; 12345678."));
        assert!(text.contains("SQ"));
        assert!(text.trim_end().ends_with("//"));
    }

    #[test]
    fn embl_round_trip_preserves_references() {
        let rec = parse(SAMPLE_REF_EMBL).unwrap();
        let text = write(&rec);
        let reparsed = parse(&text).unwrap();
        assert_eq!(reparsed.references.len(), rec.references.len());
        assert_eq!(reparsed.references[0].pubmed, "12345678");
        assert_eq!(reparsed.references[0].title, "A synthetic test article");
        assert_eq!(
            reparsed.references[1].consortium,
            "Synthetic Sequencing Consortium"
        );
    }

    #[test]
    fn embl_round_trip_preserves_features_and_sequence() {
        let rec = parse(SAMPLE_REF_EMBL).unwrap();
        let text = write(&rec);
        let reparsed = parse(&text).unwrap();
        assert_eq!(rec.seq.as_str(), reparsed.seq.as_str());
        assert_eq!(rec.features.len(), reparsed.features.len());
    }
}
