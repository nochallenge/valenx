//! BED interval-format parser and writer.
//!
//! BED stores genomic intervals as tab-delimited lines with 3 to 12
//! columns. The first three (`chrom`, `chromStart`, `chromEnd`) are
//! mandatory; the rest (`name`, `score`, `strand`, `thickStart`,
//! `thickEnd`, `itemRgb`, `blockCount`, `blockSizes`, `blockStarts`)
//! are optional and positional. BED coordinates are **0-based,
//! half-open** — `[start, end)` — unlike VCF / SAM (1-based).
//!
//! This module implements [`BedRecord`] and [`BedFile`] with the full
//! BED12 column set, a strand helper, and length / overlap utilities.

use crate::error::{GenomicsError, Result};

/// Strand of a BED feature.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum BedStrand {
    /// `+` — forward / plus strand.
    Plus,
    /// `-` — reverse / minus strand.
    Minus,
    /// `.` — strand unknown or not applicable.
    #[default]
    Unknown,
}

impl BedStrand {
    /// The single-character BED code.
    pub fn code(self) -> char {
        match self {
            BedStrand::Plus => '+',
            BedStrand::Minus => '-',
            BedStrand::Unknown => '.',
        }
    }

    /// Parses a BED strand character.
    pub fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "+" => BedStrand::Plus,
            "-" => BedStrand::Minus,
            "." => BedStrand::Unknown,
            other => {
                return Err(GenomicsError::parse(
                    "bed",
                    format!("strand must be +/-/., got `{other}`"),
                ))
            }
        })
    }
}

/// One BED record. Columns 4-12 are optional; absent columns are
/// `None`. Coordinates are 0-based half-open.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BedRecord {
    /// `chrom` — the reference contig name.
    pub chrom: String,
    /// `chromStart` — 0-based inclusive start.
    pub start: u64,
    /// `chromEnd` — 0-based exclusive end (`> start`).
    pub end: u64,
    /// `name` — the feature name (col 4).
    pub name: Option<String>,
    /// `score` — `0..=1000` display score (col 5).
    pub score: Option<i64>,
    /// `strand` (col 6).
    pub strand: Option<BedStrand>,
    /// `thickStart` — 0-based start of thick drawing (col 7).
    pub thick_start: Option<u64>,
    /// `thickEnd` — 0-based end of thick drawing (col 8).
    pub thick_end: Option<u64>,
    /// `itemRgb` — `r,g,b` display colour (col 9), kept verbatim.
    pub item_rgb: Option<String>,
    /// `blockCount` — number of sub-blocks (col 10).
    pub block_count: Option<u64>,
    /// `blockSizes` — comma-list of block sizes (col 11), verbatim.
    pub block_sizes: Option<String>,
    /// `blockStarts` — comma-list of block starts (col 12), verbatim.
    pub block_starts: Option<String>,
}

impl BedRecord {
    /// Builds a minimal BED3 interval.
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Result<Self> {
        if end <= start {
            return Err(GenomicsError::invalid_record(
                "bed",
                format!("end ({end}) must be > start ({start})"),
            ));
        }
        Ok(BedRecord {
            chrom: chrom.into(),
            start,
            end,
            name: None,
            score: None,
            strand: None,
            thick_start: None,
            thick_end: None,
            item_rgb: None,
            block_count: None,
            block_sizes: None,
            block_starts: None,
        })
    }

    /// Length of the interval in bases (`end - start`).
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// `true` for a zero-length interval (never produced by the
    /// validating constructors, but possible after manual edits).
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    /// `true` when this record overlaps `other` on the same contig.
    /// Half-open intervals overlap when `a.start < b.end &&
    /// b.start < a.end`.
    pub fn overlaps(&self, other: &BedRecord) -> bool {
        self.chrom == other.chrom
            && self.start < other.end
            && other.start < self.end
    }

    /// Parses one tab-delimited BED line.
    pub fn parse_line(line: &str) -> Result<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            return Err(GenomicsError::parse(
                "bed",
                format!("line has {} columns, need >= 3", f.len()),
            ));
        }
        let start = f[1]
            .parse::<u64>()
            .map_err(|_| GenomicsError::parse("bed", "bad chromStart"))?;
        let end = f[2]
            .parse::<u64>()
            .map_err(|_| GenomicsError::parse("bed", "bad chromEnd"))?;
        let mut rec = BedRecord::new(f[0], start, end)?;
        if let Some(v) = f.get(3) {
            rec.name = Some(v.to_string());
        }
        if let Some(v) = f.get(4) {
            if *v != "." {
                rec.score = Some(
                    v.parse::<i64>()
                        .map_err(|_| GenomicsError::parse("bed", "bad score"))?,
                );
            }
        }
        if let Some(v) = f.get(5) {
            rec.strand = Some(BedStrand::parse(v)?);
        }
        if let Some(v) = f.get(6) {
            rec.thick_start = v.parse::<u64>().ok();
        }
        if let Some(v) = f.get(7) {
            rec.thick_end = v.parse::<u64>().ok();
        }
        if let Some(v) = f.get(8) {
            rec.item_rgb = Some(v.to_string());
        }
        if let Some(v) = f.get(9) {
            rec.block_count = v.parse::<u64>().ok();
        }
        if let Some(v) = f.get(10) {
            rec.block_sizes = Some(v.to_string());
        }
        if let Some(v) = f.get(11) {
            rec.block_starts = Some(v.to_string());
        }
        Ok(rec)
    }

    /// Renders the record to a tab-delimited BED line. Trailing absent
    /// optional columns are omitted; a present later column forces
    /// placeholders for the earlier absent ones.
    pub fn to_bed_line(&self) -> String {
        let mut cols: Vec<String> = vec![
            self.chrom.clone(),
            self.start.to_string(),
            self.end.to_string(),
        ];
        // Build the optional tail, tracking the last set index.
        let opt: Vec<Option<String>> = vec![
            self.name.clone(),
            self.score.map(|s| s.to_string()),
            self.strand.map(|s| s.code().to_string()),
            self.thick_start.map(|v| v.to_string()),
            self.thick_end.map(|v| v.to_string()),
            self.item_rgb.clone(),
            self.block_count.map(|v| v.to_string()),
            self.block_sizes.clone(),
            self.block_starts.clone(),
        ];
        let last_set = opt.iter().rposition(|o| o.is_some());
        if let Some(last) = last_set {
            for o in opt.iter().take(last + 1) {
                cols.push(o.clone().unwrap_or_else(|| ".".to_string()));
            }
        }
        cols.join("\t")
    }
}

/// A parsed BED file — a flat list of records (BED has no header,
/// though `track`/`browser` lines may appear; they are kept verbatim).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BedFile {
    /// Leading `track` / `browser` directive lines, verbatim.
    pub directives: Vec<String>,
    /// Every interval record, in file order.
    pub records: Vec<BedRecord>,
}

impl BedFile {
    /// An empty BED file.
    pub fn new() -> Self {
        BedFile::default()
    }

    /// Parses an entire BED document.
    pub fn parse(text: &str) -> Result<Self> {
        let mut directives = Vec::new();
        let mut records = Vec::new();
        for line in text.lines() {
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with("track")
                || trimmed.starts_with("browser")
                || trimmed.starts_with('#')
            {
                directives.push(trimmed.to_string());
                continue;
            }
            records.push(BedRecord::parse_line(trimmed)?);
        }
        Ok(BedFile {
            directives,
            records,
        })
    }

    /// Renders the file back to BED text (with a trailing newline).
    pub fn to_bed_string(&self) -> String {
        let mut s = String::new();
        for d in &self.directives {
            s.push_str(d);
            s.push('\n');
        }
        for r in &self.records {
            s.push_str(&r.to_bed_line());
            s.push('\n');
        }
        s
    }

    /// Total bases covered by all intervals (overlaps double-counted).
    pub fn total_span(&self) -> u64 {
        self.records.iter().map(|r| r.len()).sum()
    }

    /// Sorts records by `(chrom, start, end)` in place.
    pub fn sort(&mut self) {
        self.records
            .sort_by(|a, b| (&a.chrom, a.start, a.end).cmp(&(&b.chrom, b.start, b.end)));
    }

    /// Merges overlapping / book-ended same-contig intervals into a new
    /// `BedFile` of BED3 records. The receiver is sorted first.
    pub fn merge(&self) -> BedFile {
        let mut sorted = self.clone();
        sorted.sort();
        let mut out: Vec<BedRecord> = Vec::new();
        for r in sorted.records {
            match out.last_mut() {
                Some(last)
                    if last.chrom == r.chrom && r.start <= last.end =>
                {
                    last.end = last.end.max(r.end);
                }
                _ => {
                    // Re-emit as a clean BED3 record.
                    if let Ok(b) = BedRecord::new(r.chrom, r.start, r.end) {
                        out.push(b);
                    }
                }
            }
        }
        BedFile {
            directives: Vec::new(),
            records: out,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "track name=demo\n\
chr1\t100\t200\tfeatA\t500\t+\n\
chr1\t150\t250\tfeatB\t600\t-\n\
chr1\t400\t500\tfeatC\t100\t.\n\
chr2\t0\t1000\n";

    #[test]
    fn parse_full_file() {
        let f = BedFile::parse(SAMPLE).unwrap();
        assert_eq!(f.directives.len(), 1);
        assert_eq!(f.records.len(), 4);
        assert_eq!(f.records[0].name.as_deref(), Some("featA"));
        assert_eq!(f.records[0].strand, Some(BedStrand::Plus));
        assert_eq!(f.records[3].chrom, "chr2");
        assert_eq!(f.records[3].len(), 1000);
    }

    #[test]
    fn coordinates_are_half_open() {
        let r = BedRecord::new("chr1", 10, 20).unwrap();
        assert_eq!(r.len(), 10);
    }

    #[test]
    fn rejects_end_le_start() {
        assert!(BedRecord::new("chr1", 20, 10).is_err());
        assert!(BedRecord::new("chr1", 10, 10).is_err());
    }

    #[test]
    fn overlap_detection() {
        let f = BedFile::parse(SAMPLE).unwrap();
        assert!(f.records[0].overlaps(&f.records[1])); // 100-200 vs 150-250
        assert!(!f.records[0].overlaps(&f.records[2])); // 100-200 vs 400-500
        assert!(!f.records[0].overlaps(&f.records[3])); // different chrom
    }

    #[test]
    fn merge_collapses_overlaps() {
        let f = BedFile::parse(SAMPLE).unwrap();
        let m = f.merge();
        // chr1: 100-200 + 150-250 -> 100-250; 400-500 separate. chr2 one.
        assert_eq!(m.records.len(), 3);
        assert_eq!((m.records[0].start, m.records[0].end), (100, 250));
        assert_eq!((m.records[1].start, m.records[1].end), (400, 500));
    }

    #[test]
    fn round_trip_preserves_records() {
        let f = BedFile::parse(SAMPLE).unwrap();
        let text = f.to_bed_string();
        let f2 = BedFile::parse(&text).unwrap();
        assert_eq!(f.records, f2.records);
    }

    #[test]
    fn total_span_sums_lengths() {
        let f = BedFile::parse(SAMPLE).unwrap();
        assert_eq!(f.total_span(), 100 + 100 + 100 + 1000);
    }
}
