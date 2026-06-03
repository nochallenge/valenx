//! GFF3 and GTF feature-annotation parser.
//!
//! GFF3 (Generic Feature Format v3) and GTF (Gene Transfer Format,
//! a.k.a. GFF 2.5) share the same nine tab-delimited columns —
//! `seqid`, `source`, `type`, `start`, `end`, `score`, `strand`,
//! `phase` and a ninth **attributes** column. The only difference is
//! that column-9 syntax:
//!
//! - **GFF3** — `key=value;key=value`, values percent-encoded.
//! - **GTF** — `key "value"; key "value";`, values double-quoted.
//!
//! Coordinates are **1-based, inclusive** in both formats (unlike BED).
//! This module implements [`GffRecord`], a [`GffDialect`] selector, and
//! [`GffFile`]; the parser auto-detects the dialect from the attribute
//! syntax when asked.

use crate::error::{GenomicsError, Result};

/// Which column-9 attribute syntax a record uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum GffDialect {
    /// GFF3 — `key=value;` attributes.
    #[default]
    Gff3,
    /// GTF / GFF2.5 — `key "value";` attributes.
    Gtf,
}

/// Strand of a GFF / GTF feature.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum GffStrand {
    /// `+` — forward strand.
    Plus,
    /// `-` — reverse strand.
    Minus,
    /// `.` — strand not applicable.
    #[default]
    Unstranded,
    /// `?` — strand relevant but unknown.
    Unknown,
}

impl GffStrand {
    /// The single-character code.
    pub fn code(self) -> char {
        match self {
            GffStrand::Plus => '+',
            GffStrand::Minus => '-',
            GffStrand::Unstranded => '.',
            GffStrand::Unknown => '?',
        }
    }

    /// Parses a strand character.
    pub fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "+" => GffStrand::Plus,
            "-" => GffStrand::Minus,
            "." => GffStrand::Unstranded,
            "?" => GffStrand::Unknown,
            other => {
                return Err(GenomicsError::parse(
                    "gff",
                    format!("strand must be +/-/./?, got `{other}`"),
                ))
            }
        })
    }
}

/// One GFF3 / GTF feature record.
#[derive(Clone, Debug, PartialEq)]
pub struct GffRecord {
    /// `seqid` — the landmark / contig name.
    pub seqid: String,
    /// `source` — the annotation source / algorithm.
    pub source: String,
    /// `type` — the feature type (`gene`, `exon`, `CDS`, …).
    pub feature_type: String,
    /// `start` — 1-based inclusive start.
    pub start: u64,
    /// `end` — 1-based inclusive end (`>= start`).
    pub end: u64,
    /// `score` — a numeric score (`None` for `.`).
    pub score: Option<f64>,
    /// `strand`.
    pub strand: GffStrand,
    /// `phase` — `0`, `1`, `2` for CDS features (`None` for `.`).
    pub phase: Option<u8>,
    /// The ordered column-9 attributes as key→value pairs.
    pub attributes: Vec<(String, String)>,
    /// Which dialect the record was parsed as.
    pub dialect: GffDialect,
}

impl GffRecord {
    /// Length of the feature in bases — inclusive coordinates, so
    /// `end - start + 1`.
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start) + 1
    }

    /// `true` for a degenerate `end < start` record.
    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }

    /// Looks up the first value of a column-9 attribute by key.
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// The GFF3 `ID` attribute, or the GTF `gene_id`, whichever the
    /// dialect uses for the primary identifier.
    pub fn primary_id(&self) -> Option<&str> {
        match self.dialect {
            GffDialect::Gff3 => self.attr("ID"),
            GffDialect::Gtf => self.attr("gene_id"),
        }
    }

    /// Parses one tab-delimited line with the given `dialect`.
    pub fn parse_line(line: &str, dialect: GffDialect) -> Result<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            return Err(GenomicsError::parse(
                "gff",
                format!("line has {} columns, need 9", f.len()),
            ));
        }
        let start = f[3]
            .parse::<u64>()
            .map_err(|_| GenomicsError::parse("gff", "bad start"))?;
        let end = f[4]
            .parse::<u64>()
            .map_err(|_| GenomicsError::parse("gff", "bad end"))?;
        if end < start {
            return Err(GenomicsError::invalid_record(
                "gff",
                format!("end ({end}) < start ({start})"),
            ));
        }
        let score = if f[5] == "." {
            None
        } else {
            Some(
                f[5].parse::<f64>()
                    .map_err(|_| GenomicsError::parse("gff", "bad score"))?,
            )
        };
        let strand = GffStrand::parse(f[6])?;
        let phase = if f[7] == "." {
            None
        } else {
            let p = f[7]
                .parse::<u8>()
                .map_err(|_| GenomicsError::parse("gff", "bad phase"))?;
            if p > 2 {
                return Err(GenomicsError::parse("gff", "phase must be 0, 1 or 2"));
            }
            Some(p)
        };
        let attributes = match dialect {
            GffDialect::Gff3 => parse_gff3_attrs(f[8]),
            GffDialect::Gtf => parse_gtf_attrs(f[8]),
        };
        Ok(GffRecord {
            seqid: f[0].to_string(),
            source: f[1].to_string(),
            feature_type: f[2].to_string(),
            start,
            end,
            score,
            strand,
            phase,
            attributes,
            dialect,
        })
    }

    /// Renders the record back to a tab-delimited line in its dialect.
    pub fn to_line(&self) -> String {
        let score = self
            .score
            .map(|s| {
                if s.fract() == 0.0 {
                    format!("{}", s as i64)
                } else {
                    format!("{s}")
                }
            })
            .unwrap_or_else(|| ".".to_string());
        let phase = self
            .phase
            .map(|p| p.to_string())
            .unwrap_or_else(|| ".".to_string());
        let attrs = match self.dialect {
            GffDialect::Gff3 => self
                .attributes
                .iter()
                .map(|(k, v)| format!("{}={}", k, gff3_encode(v)))
                .collect::<Vec<_>>()
                .join(";"),
            GffDialect::Gtf => {
                let mut s = self
                    .attributes
                    .iter()
                    .map(|(k, v)| format!("{k} \"{v}\""))
                    .collect::<Vec<_>>()
                    .join("; ");
                if !s.is_empty() {
                    s.push(';');
                }
                s
            }
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seqid,
            self.source,
            self.feature_type,
            self.start,
            self.end,
            score,
            self.strand.code(),
            phase,
            attrs,
        )
    }
}

/// Parses a GFF3 `key=value;key=value` attribute column.
fn parse_gff3_attrs(col: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    if col == "." || col.is_empty() {
        return out;
    }
    for kv in col.split(';') {
        let kv = kv.trim();
        if kv.is_empty() {
            continue;
        }
        if let Some((k, v)) = kv.split_once('=') {
            out.push((k.to_string(), gff3_decode(v)));
        } else {
            out.push((kv.to_string(), String::new()));
        }
    }
    out
}

/// Parses a GTF `key "value"; key "value";` attribute column.
fn parse_gtf_attrs(col: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    if col == "." || col.is_empty() {
        return out;
    }
    for kv in col.split(';') {
        let kv = kv.trim();
        if kv.is_empty() {
            continue;
        }
        // `key "value"` — split at the first space.
        if let Some((k, v)) = kv.split_once(char::is_whitespace) {
            let v = v.trim().trim_matches('"');
            out.push((k.to_string(), v.to_string()));
        } else {
            out.push((kv.to_string(), String::new()));
        }
    }
    out
}

/// Percent-decodes the GFF3-reserved characters in an attribute value.
///
/// `%XX` sequences decode to a single byte; the surrounding text is
/// copied codepoint-by-codepoint so multi-byte UTF-8 (a gene symbol
/// like "α-tubulin" carrying U+03B1 GREEK SMALL LETTER ALPHA) survives
/// the round-trip intact. The previous byte-by-byte `bytes[i] as char`
/// fast path would splatter mojibake into the output (every byte above
/// 0x7F became its own Latin-1 char).
fn gff3_decode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            // Hex-digit lookahead is single-byte ASCII either way —
            // `as char` is safe here.
            let hi = (bytes[i + 1] as char).to_digit(16);
            let lo = (bytes[i + 2] as char).to_digit(16);
            if let (Some(h), Some(l)) = (hi, lo) {
                out.push((h * 16 + l) as u8 as char);
                i += 3;
                continue;
            }
        }
        // Copy the next full UTF-8 codepoint verbatim.
        let ch = s[i..]
            .chars()
            .next()
            .expect("i is in-bounds inside a &str");
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

/// Percent-encodes the GFF3-reserved characters (`;`, `=`, `&`, `,`,
/// tab, newline) in an attribute value.
fn gff3_encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            ';' | '=' | '&' | ',' | '\t' | '\n' | '\r' | '%' => {
                out.push_str(&format!("%{:02X}", c as u32));
            }
            _ => out.push(c),
        }
    }
    out
}

/// A parsed GFF3 / GTF file — directive lines plus feature records.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GffFile {
    /// Leading `##`/`#!` directive and `#` comment lines, verbatim.
    pub directives: Vec<String>,
    /// Every feature record, in file order.
    pub records: Vec<GffRecord>,
    /// The dialect used for the whole file.
    pub dialect: GffDialect,
}

impl GffFile {
    /// Parses a document with an explicit `dialect`.
    pub fn parse(text: &str, dialect: GffDialect) -> Result<Self> {
        let mut directives = Vec::new();
        let mut records = Vec::new();
        for line in text.lines() {
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('#') {
                directives.push(trimmed.to_string());
                continue;
            }
            records.push(GffRecord::parse_line(trimmed, dialect)?);
        }
        Ok(GffFile {
            directives,
            records,
            dialect,
        })
    }

    /// Parses a document, auto-detecting the dialect: a column-9 token
    /// containing `="` GTF-style quoting or no `=` at all but a quoted
    /// value is treated as GTF; otherwise GFF3. The heuristic inspects
    /// the first data line.
    pub fn parse_autodetect(text: &str) -> Result<Self> {
        let dialect = detect_dialect(text);
        Self::parse(text, dialect)
    }

    /// Renders the file back to text (with a trailing newline).
    pub fn to_string_lines(&self) -> String {
        let mut s = String::new();
        for d in &self.directives {
            s.push_str(d);
            s.push('\n');
        }
        for r in &self.records {
            s.push_str(&r.to_line());
            s.push('\n');
        }
        s
    }

    /// All records whose `type` equals `feature_type`.
    pub fn by_type<'a>(
        &'a self,
        feature_type: &'a str,
    ) -> impl Iterator<Item = &'a GffRecord> {
        self.records
            .iter()
            .filter(move |r| r.feature_type == feature_type)
    }
}

/// Inspects the first data line of `text` and guesses the dialect.
fn detect_dialect(text: &str) -> GffDialect {
    if text.lines().any(|l| l.starts_with("##gff-version 3")) {
        return GffDialect::Gff3;
    }
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if let Some(col9) = t.split('\t').nth(8) {
            // GTF tokens look like `gene_id "x";`; GFF3 like `ID=x;`.
            let gtf_like = col9.contains('"')
                && col9
                    .split(';')
                    .filter(|s| !s.trim().is_empty())
                    .any(|s| {
                        let s = s.trim();
                        s.split_once(char::is_whitespace)
                            .map(|(_, v)| v.trim().starts_with('"'))
                            .unwrap_or(false)
                    });
            return if gtf_like {
                GffDialect::Gtf
            } else {
                GffDialect::Gff3
            };
        }
        break;
    }
    GffDialect::Gff3
}

#[cfg(test)]
mod tests {
    use super::*;

    const GFF3: &str = "##gff-version 3\n\
chr1\tphytozome\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC\n\
chr1\tphytozome\tmRNA\t1000\t2000\t.\t+\t.\tID=mrna1;Parent=gene1\n\
chr1\tphytozome\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=mrna1\n\
chr1\tphytozome\tCDS\t1050\t1200\t.\t+\t0\tID=cds1;Parent=mrna1\n";

    const GTF: &str = "chr1\tHAVANA\tgene\t100\t900\t.\t+\t.\tgene_id \"G1\"; gene_name \"ABC\";\n\
chr1\tHAVANA\texon\t100\t300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n";

    #[test]
    fn parse_gff3() {
        let f = GffFile::parse(GFF3, GffDialect::Gff3).unwrap();
        assert_eq!(f.records.len(), 4);
        assert_eq!(f.records[0].feature_type, "gene");
        assert_eq!(f.records[0].attr("ID"), Some("gene1"));
        assert_eq!(f.records[0].attr("Name"), Some("ABC"));
        assert_eq!(f.records[0].primary_id(), Some("gene1"));
        assert_eq!(f.records[3].phase, Some(0));
        assert_eq!(f.records[0].len(), 1001); // inclusive
    }

    #[test]
    fn parse_gtf() {
        let f = GffFile::parse(GTF, GffDialect::Gtf).unwrap();
        assert_eq!(f.records.len(), 2);
        assert_eq!(f.records[0].attr("gene_id"), Some("G1"));
        assert_eq!(f.records[0].attr("gene_name"), Some("ABC"));
        assert_eq!(f.records[0].primary_id(), Some("G1"));
    }

    #[test]
    fn autodetect_picks_dialect() {
        assert_eq!(
            GffFile::parse_autodetect(GFF3).unwrap().dialect,
            GffDialect::Gff3
        );
        assert_eq!(
            GffFile::parse_autodetect(GTF).unwrap().dialect,
            GffDialect::Gtf
        );
    }

    #[test]
    fn rejects_end_before_start() {
        let bad = "chr1\ts\tgene\t2000\t1000\t.\t+\t.\tID=x\n";
        assert!(GffFile::parse(bad, GffDialect::Gff3).is_err());
    }

    #[test]
    fn rejects_bad_phase() {
        let bad = "chr1\ts\tCDS\t10\t20\t.\t+\t5\tID=x\n";
        assert!(GffFile::parse(bad, GffDialect::Gff3).is_err());
    }

    #[test]
    fn by_type_filters() {
        let f = GffFile::parse(GFF3, GffDialect::Gff3).unwrap();
        assert_eq!(f.by_type("exon").count(), 1);
        assert_eq!(f.by_type("gene").count(), 1);
    }

    #[test]
    fn gff3_percent_encoding_round_trips() {
        let mut rec = GffRecord::parse_line(
            "chr1\ts\tgene\t1\t2\t.\t+\t.\tID=x;Note=a%3Bb",
            GffDialect::Gff3,
        )
        .unwrap();
        assert_eq!(rec.attr("Note"), Some("a;b"));
        // Re-render must re-encode the semicolon.
        let line = rec.to_line();
        assert!(line.contains("a%3Bb"));
        rec.attributes.clear();
        assert!(rec.attr("Note").is_none());
    }

    #[test]
    fn round_trip_gff3() {
        let f = GffFile::parse(GFF3, GffDialect::Gff3).unwrap();
        let text = f.to_string_lines();
        let f2 = GffFile::parse(&text, GffDialect::Gff3).unwrap();
        assert_eq!(f.records, f2.records);
    }
}
