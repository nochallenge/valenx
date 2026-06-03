//! VCF (Variant Call Format) parser and writer.
//!
//! VCF is the standard text format for genetic variants: a `##`-line
//! metadata block, one `#CHROM …` column-header line, and a body of
//! tab-delimited variant records each with eight mandatory columns
//! (`CHROM POS ID REF ALT QUAL FILTER INFO`) plus, when samples are
//! present, a `FORMAT` column and one genotype column per sample.
//!
//! This module implements [`VcfRecord`] (one variant site),
//! [`VcfHeader`] (the typed metadata) and [`VcfFile`] (the pair). It
//! covers VCF 4.x text; the BCF binary encoding is out of v1 scope
//! (the same BGZF gap noted in [`crate::format::sam`]).

use crate::error::{GenomicsError, Result};
use std::collections::BTreeMap;

// --------------------------------------------------------------------
// Header metadata
// --------------------------------------------------------------------

/// A typed `##INFO` / `##FORMAT` / `##FILTER` / `##contig` metadata
/// line. The fields are kept as an ordered key→value map so an unusual
/// attribute survives the round-trip.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetaLine {
    /// The line kind without the `##` (`"INFO"`, `"FORMAT"`, `"contig"`,
    /// `"FILTER"`, `"ALT"`, …).
    pub kind: String,
    /// Whether the value was an angle-bracket attribute block
    /// (`<ID=…,…>`) — `false` for a plain `##key=value` line.
    pub structured: bool,
    /// For a plain line: the raw value. For a structured line: empty.
    pub raw_value: String,
    /// For a structured line: the ordered `ID=…,Number=…` attributes.
    pub attrs: Vec<(String, String)>,
}

impl MetaLine {
    /// Looks up a structured attribute by key.
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// Parses one `##…` line.
    pub fn parse(line: &str) -> Result<Self> {
        let body = line
            .strip_prefix("##")
            .ok_or_else(|| GenomicsError::parse("vcf", "meta line must start with `##`"))?;
        let (kind, value) = body
            .split_once('=')
            .ok_or_else(|| GenomicsError::parse("vcf", "meta line missing `=`"))?;
        if let Some(inner) = value.strip_prefix('<').and_then(|v| v.strip_suffix('>')) {
            let attrs = parse_attr_block(inner)?;
            Ok(MetaLine {
                kind: kind.to_string(),
                structured: true,
                raw_value: String::new(),
                attrs,
            })
        } else {
            Ok(MetaLine {
                kind: kind.to_string(),
                structured: false,
                raw_value: value.to_string(),
                attrs: Vec::new(),
            })
        }
    }

    /// Renders the line back to `##…` text.
    pub fn to_vcf_string(&self) -> String {
        if self.structured {
            let inner: Vec<String> = self
                .attrs
                .iter()
                .map(|(k, v)| {
                    // Re-quote the Description value (it nearly always
                    // contains spaces).
                    if v.contains(' ') || v.contains(',') {
                        format!("{k}=\"{v}\"")
                    } else {
                        format!("{k}={v}")
                    }
                })
                .collect();
            format!("##{}=<{}>", self.kind, inner.join(","))
        } else {
            format!("##{}={}", self.kind, self.raw_value)
        }
    }
}

/// Splits a `<…>` attribute block respecting double-quoted values
/// (`Description="a, b"` must not split at the inner comma).
fn parse_attr_block(inner: &str) -> Result<Vec<(String, String)>> {
    let mut attrs = Vec::new();
    let mut chars = inner.chars().peekable();
    while chars.peek().is_some() {
        // key
        let mut key = String::new();
        for c in chars.by_ref() {
            if c == '=' {
                break;
            }
            key.push(c);
        }
        if key.is_empty() {
            break;
        }
        // value (quoted or bare)
        let mut value = String::new();
        if chars.peek() == Some(&'"') {
            chars.next(); // opening quote
            let mut escaped = false;
            for c in chars.by_ref() {
                if escaped {
                    value.push(c);
                    escaped = false;
                } else if c == '\\' {
                    escaped = true;
                } else if c == '"' {
                    break;
                } else {
                    value.push(c);
                }
            }
            // Consume the comma after the closing quote, if any.
            if chars.peek() == Some(&',') {
                chars.next();
            }
        } else {
            for c in chars.by_ref() {
                if c == ',' {
                    break;
                }
                value.push(c);
            }
        }
        attrs.push((key.trim().to_string(), value));
    }
    if attrs.is_empty() {
        return Err(GenomicsError::parse("vcf", "empty `<…>` attribute block"));
    }
    Ok(attrs)
}

/// The typed VCF header: the fileformat line, every `##` meta line, and
/// the sample-name list from the `#CHROM` column-header line.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct VcfHeader {
    /// The `##fileformat` value (e.g. `"VCFv4.2"`).
    pub fileformat: String,
    /// Every `##` meta line, in file order.
    pub meta: Vec<MetaLine>,
    /// The sample names from the `#CHROM` line (empty for a sites-only
    /// VCF with no `FORMAT` column).
    pub samples: Vec<String>,
}

impl VcfHeader {
    /// A minimal header for a VCF 4.2 sites-only file.
    pub fn sites_only() -> Self {
        VcfHeader {
            fileformat: "VCFv4.2".to_string(),
            meta: Vec::new(),
            samples: Vec::new(),
        }
    }

    /// All `##INFO` meta lines.
    pub fn info_lines(&self) -> impl Iterator<Item = &MetaLine> {
        self.meta.iter().filter(|m| m.kind == "INFO")
    }

    /// All `##FORMAT` meta lines.
    pub fn format_lines(&self) -> impl Iterator<Item = &MetaLine> {
        self.meta.iter().filter(|m| m.kind == "FORMAT")
    }

    /// Renders the header (all `##` lines plus the `#CHROM` line),
    /// newline-joined, no trailing newline.
    pub fn to_vcf_string(&self) -> String {
        let mut out = vec![format!("##fileformat={}", self.fileformat)];
        for m in &self.meta {
            out.push(m.to_vcf_string());
        }
        let mut chrom =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
        if !self.samples.is_empty() {
            chrom.push_str("\tFORMAT");
            for s in &self.samples {
                chrom.push('\t');
                chrom.push_str(s);
            }
        }
        out.push(chrom);
        out.join("\n")
    }
}

// --------------------------------------------------------------------
// Variant record
// --------------------------------------------------------------------

/// One VCF variant record.
#[derive(Clone, Debug, PartialEq)]
pub struct VcfRecord {
    /// `CHROM` — the reference contig name.
    pub chrom: String,
    /// `POS` — the 1-based reference position.
    pub pos: i64,
    /// `ID` — semicolon-joined variant identifiers (empty for `.`).
    pub id: Vec<String>,
    /// `REF` — the reference allele.
    pub reference: String,
    /// `ALT` — the alternate alleles (empty for `.`).
    pub alt: Vec<String>,
    /// `QUAL` — the phred-scaled site quality (`None` for `.`).
    pub qual: Option<f64>,
    /// `FILTER` — the filter labels (`["PASS"]` or specific failures;
    /// empty for `.`).
    pub filter: Vec<String>,
    /// `INFO` — the parsed key→value attributes. A flag attribute maps
    /// to an empty string.
    pub info: BTreeMap<String, String>,
    /// `FORMAT` — the per-sample field keys (empty for a sites-only
    /// record).
    pub format: Vec<String>,
    /// One genotype column per sample; each is the `:`-split field
    /// values aligned to [`format`](Self::format).
    pub samples: Vec<Vec<String>>,
}

impl VcfRecord {
    /// Builds a minimal SNV record at `chrom:pos` with one ALT allele.
    pub fn snv(
        chrom: impl Into<String>,
        pos: i64,
        reference: impl Into<String>,
        alt: impl Into<String>,
    ) -> Self {
        VcfRecord {
            chrom: chrom.into(),
            pos,
            id: Vec::new(),
            reference: reference.into(),
            alt: vec![alt.into()],
            qual: None,
            filter: Vec::new(),
            info: BTreeMap::new(),
            format: Vec::new(),
            samples: Vec::new(),
        }
    }

    /// `true` when every REF and ALT allele is a single base — a pure
    /// SNV / MNV site with no indel.
    pub fn is_snv(&self) -> bool {
        self.reference.len() == 1 && self.alt.iter().all(|a| a.len() == 1)
    }

    /// `true` when any ALT allele differs in length from REF — an
    /// insertion or deletion is present.
    pub fn is_indel(&self) -> bool {
        self.alt.iter().any(|a| a.len() != self.reference.len())
    }

    /// `true` when more than one ALT allele is recorded.
    pub fn is_multiallelic(&self) -> bool {
        self.alt.len() > 1
    }

    /// Reads an INFO attribute by key. A present flag returns
    /// `Some("")`.
    pub fn info_get(&self, key: &str) -> Option<&str> {
        self.info.get(key).map(|s| s.as_str())
    }

    /// Reads one sample's value for a `FORMAT` field, e.g. `"GT"`.
    pub fn sample_field(&self, sample_idx: usize, field: &str) -> Option<&str> {
        let fi = self.format.iter().position(|f| f == field)?;
        self.samples.get(sample_idx)?.get(fi).map(|s| s.as_str())
    }

    /// Validates that every genotype column has at most as many fields
    /// as `FORMAT` declares — a corrupt-record guard.
    pub fn validate(&self) -> Result<()> {
        for (i, s) in self.samples.iter().enumerate() {
            if s.len() > self.format.len() {
                return Err(GenomicsError::invalid_record(
                    "vcf",
                    format!(
                        "sample {i} has {} fields, FORMAT declares {}",
                        s.len(),
                        self.format.len()
                    ),
                ));
            }
        }
        Ok(())
    }

    /// Parses one tab-delimited VCF body line. `n_samples` is the count
    /// declared by the header — used to validate the column count.
    pub fn parse_line(line: &str, n_samples: usize) -> Result<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        let need = if n_samples > 0 { 9 + n_samples } else { 8 };
        if f.len() < need {
            return Err(GenomicsError::parse(
                "vcf",
                format!("record has {} columns, need {need}", f.len()),
            ));
        }
        let pos = f[1]
            .parse::<i64>()
            .map_err(|_| GenomicsError::parse("vcf", "bad POS"))?;
        let id = split_dot(f[2], ';');
        let reference = f[3].to_string();
        let alt = split_dot(f[4], ',');
        let qual = if f[5] == "." {
            None
        } else {
            Some(
                f[5].parse::<f64>()
                    .map_err(|_| GenomicsError::parse("vcf", "bad QUAL"))?,
            )
        };
        let filter = split_dot(f[6], ';');
        let info = parse_info(f[7]);
        let (format, samples) = if n_samples > 0 {
            let fmt: Vec<String> = f[8].split(':').map(|s| s.to_string()).collect();
            let mut samples = Vec::with_capacity(n_samples);
            for col in &f[9..9 + n_samples] {
                samples.push(col.split(':').map(|s| s.to_string()).collect());
            }
            (fmt, samples)
        } else {
            (Vec::new(), Vec::new())
        };
        let rec = VcfRecord {
            chrom: f[0].to_string(),
            pos,
            id,
            reference,
            alt,
            qual,
            filter,
            info,
            format,
            samples,
        };
        rec.validate()?;
        Ok(rec)
    }

    /// Renders the record back to one tab-delimited VCF line.
    pub fn to_vcf_line(&self) -> String {
        let id = join_or_dot(&self.id, ";");
        let alt = join_or_dot(&self.alt, ",");
        let qual = match self.qual {
            Some(q) => format_qual(q),
            None => ".".to_string(),
        };
        let filter = join_or_dot(&self.filter, ";");
        let info = render_info(&self.info);
        let mut s = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom, self.pos, id, self.reference, alt, qual, filter, info
        );
        if !self.format.is_empty() {
            s.push('\t');
            s.push_str(&self.format.join(":"));
            for col in &self.samples {
                s.push('\t');
                s.push_str(&col.join(":"));
            }
        }
        s
    }
}

fn split_dot(field: &str, sep: char) -> Vec<String> {
    if field == "." || field.is_empty() {
        Vec::new()
    } else {
        field.split(sep).map(|s| s.to_string()).collect()
    }
}

fn join_or_dot(parts: &[String], sep: &str) -> String {
    if parts.is_empty() {
        ".".to_string()
    } else {
        parts.join(sep)
    }
}

fn parse_info(field: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    if field == "." || field.is_empty() {
        return map;
    }
    for kv in field.split(';') {
        match kv.split_once('=') {
            Some((k, v)) => {
                map.insert(k.to_string(), v.to_string());
            }
            None => {
                map.insert(kv.to_string(), String::new());
            }
        }
    }
    map
}

fn render_info(map: &BTreeMap<String, String>) -> String {
    if map.is_empty() {
        return ".".to_string();
    }
    let parts: Vec<String> = map
        .iter()
        .map(|(k, v)| {
            if v.is_empty() {
                k.clone()
            } else {
                format!("{k}={v}")
            }
        })
        .collect();
    parts.join(";")
}

/// Renders a QUAL value: integral values lose the trailing `.0`.
fn format_qual(q: f64) -> String {
    if q.fract() == 0.0 && q.abs() < 1e15 {
        format!("{}", q as i64)
    } else {
        format!("{q:.2}")
    }
}

// --------------------------------------------------------------------
// VCF file
// --------------------------------------------------------------------

/// A parsed VCF file — its header and every variant record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VcfFile {
    /// The typed header.
    pub header: VcfHeader,
    /// Every variant record, in file order.
    pub records: Vec<VcfRecord>,
}

impl VcfFile {
    /// An empty sites-only VCF.
    pub fn new() -> Self {
        VcfFile {
            header: VcfHeader::sites_only(),
            records: Vec::new(),
        }
    }

    /// Parses an entire VCF document.
    pub fn parse(text: &str) -> Result<Self> {
        let mut header = VcfHeader::default();
        let mut records = Vec::new();
        let mut seen_chrom_line = false;
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("##") {
                if let Some(ff) = rest.strip_prefix("fileformat=") {
                    header.fileformat = ff.to_string();
                } else {
                    header.meta.push(MetaLine::parse(line)?);
                }
            } else if let Some(rest) = line.strip_prefix("#CHROM") {
                // Sample columns follow the 9 fixed columns.
                let cols: Vec<&str> = rest.split('\t').collect();
                // cols[0] is empty (the text after "#CHROM"); fixed
                // columns are POS..FORMAT; samples start at index 9.
                if cols.len() > 9 {
                    header.samples =
                        cols[9..].iter().map(|s| s.to_string()).collect();
                }
                seen_chrom_line = true;
            } else if line.starts_with('#') {
                // An unrecognised comment line — tolerate it.
                continue;
            } else {
                if !seen_chrom_line {
                    return Err(GenomicsError::parse(
                        "vcf",
                        "record before the `#CHROM` header line",
                    ));
                }
                records.push(VcfRecord::parse_line(line, header.samples.len())?);
            }
        }
        if header.fileformat.is_empty() {
            header.fileformat = "VCFv4.2".to_string();
        }
        Ok(VcfFile { header, records })
    }

    /// Renders the whole file back to VCF text (with a trailing
    /// newline).
    pub fn to_vcf_string(&self) -> String {
        let mut s = self.header.to_vcf_string();
        s.push('\n');
        for rec in &self.records {
            s.push_str(&rec.to_vcf_line());
            s.push('\n');
        }
        s
    }

    /// Count of records that are pure SNV sites.
    pub fn snv_count(&self) -> usize {
        self.records.iter().filter(|r| r.is_snv()).count()
    }

    /// Count of records carrying an indel.
    pub fn indel_count(&self) -> usize {
        self.records.iter().filter(|r| r.is_indel()).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, see notes\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\n\
chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30;AF=0.5\tGT:DP\t0/1:30\n\
chr1\t200\t.\tAT\tA\t99\tPASS\tDP=40\tGT:DP\t1/1:40\n\
chr1\t300\t.\tC\tT,G\t.\tq10\tDP=12\tGT:DP\t1/2:12\n";

    #[test]
    fn parse_full_file() {
        let f = VcfFile::parse(SAMPLE).unwrap();
        assert_eq!(f.header.fileformat, "VCFv4.2");
        assert_eq!(f.header.samples, vec!["sampleA"]);
        assert_eq!(f.header.info_lines().count(), 2);
        assert_eq!(f.header.format_lines().count(), 2);
        assert_eq!(f.records.len(), 3);
        assert_eq!(f.snv_count(), 2);
        assert_eq!(f.indel_count(), 1);
    }

    #[test]
    fn record_fields_typed_correctly() {
        let f = VcfFile::parse(SAMPLE).unwrap();
        let r0 = &f.records[0];
        assert_eq!(r0.pos, 100);
        assert_eq!(r0.id, vec!["rs1"]);
        assert_eq!(r0.qual, Some(50.0));
        assert_eq!(r0.filter, vec!["PASS"]);
        assert_eq!(r0.info_get("DP"), Some("30"));
        assert_eq!(r0.sample_field(0, "GT"), Some("0/1"));
        assert_eq!(r0.sample_field(0, "DP"), Some("30"));
        assert!(r0.is_snv());

        let r1 = &f.records[1];
        assert!(r1.is_indel());
        assert!(!r1.is_snv());

        let r2 = &f.records[2];
        assert!(r2.is_multiallelic());
        assert_eq!(r2.qual, None);
    }

    #[test]
    fn meta_line_quoted_description_round_trips() {
        let ml = MetaLine::parse(
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, see notes\">",
        )
        .unwrap();
        assert!(ml.structured);
        assert_eq!(ml.attr("ID"), Some("AF"));
        assert_eq!(ml.attr("Description"), Some("Allele frequency, see notes"));
        let s = ml.to_vcf_string();
        let ml2 = MetaLine::parse(&s).unwrap();
        assert_eq!(ml, ml2);
    }

    #[test]
    fn round_trip_preserves_records() {
        let f = VcfFile::parse(SAMPLE).unwrap();
        let text = f.to_vcf_string();
        let f2 = VcfFile::parse(&text).unwrap();
        assert_eq!(f.records, f2.records);
    }

    #[test]
    fn record_before_chrom_line_is_an_error() {
        let bad = "##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t.\t.\t.\n";
        assert!(VcfFile::parse(bad).is_err());
    }

    #[test]
    fn snv_builder_and_info() {
        let mut r = VcfRecord::snv("chr2", 5, "A", "T");
        r.info.insert("DP".to_string(), "10".to_string());
        let line = r.to_vcf_line();
        assert!(line.contains("DP=10"));
        assert!(line.starts_with("chr2\t5\t."));
    }
}
