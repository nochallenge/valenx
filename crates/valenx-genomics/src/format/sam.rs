//! SAM (Sequence Alignment / Map) parser and writer.
//!
//! SAM is the de-facto text format for storing read alignments — the
//! human-readable sibling of BAM. A SAM file is a header section of
//! `@`-prefixed lines followed by a body of tab-delimited alignment
//! records, each with the eleven mandatory fields and any number of
//! optional `TAG:TYPE:VALUE` fields.
//!
//! This module implements [`SamRecord`] (one alignment), [`SamHeader`]
//! (the typed header lines) and [`SamFile`] (the pair), plus a
//! [`Cigar`] type with full SAM-spec op coverage and a flag-bit helper.
//!
//! ## v1 scope
//!
//! Text SAM only. Binary BAM (BGZF-compressed) and the CRAM
//! reference-compressed format are out of scope for the v1 — they need
//! a BGZF / DEFLATE codec that the Round-6 budget keeps off the
//! dependency tree; [`crate::format`] documents the gap. Optional tags
//! are parsed into a typed [`SamTag`] value list; the rare `B`
//! (numeric-array) tag type is stored as its raw string.

use crate::error::{GenomicsError, Result};
use std::fmt;

// --------------------------------------------------------------------
// CIGAR
// --------------------------------------------------------------------

/// A single CIGAR operation kind (the SAM spec's nine op codes).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum CigarKind {
    /// `M` — alignment match (sequence match *or* mismatch).
    Match,
    /// `I` — insertion to the reference.
    Ins,
    /// `D` — deletion from the reference.
    Del,
    /// `N` — skipped region from the reference (e.g. an mRNA intron).
    Skip,
    /// `S` — soft clip: clipped bases present in `SEQ`.
    SoftClip,
    /// `H` — hard clip: clipped bases absent from `SEQ`.
    HardClip,
    /// `P` — padding (silent deletion from a padded reference).
    Pad,
    /// `=` — sequence match.
    Equal,
    /// `X` — sequence mismatch.
    Diff,
}

impl CigarKind {
    /// The single-character SAM op code.
    pub fn code(self) -> char {
        match self {
            CigarKind::Match => 'M',
            CigarKind::Ins => 'I',
            CigarKind::Del => 'D',
            CigarKind::Skip => 'N',
            CigarKind::SoftClip => 'S',
            CigarKind::HardClip => 'H',
            CigarKind::Pad => 'P',
            CigarKind::Equal => '=',
            CigarKind::Diff => 'X',
        }
    }

    /// Parses a SAM op code; `None` for an unrecognised character.
    pub fn from_code(c: char) -> Option<Self> {
        Some(match c {
            'M' => CigarKind::Match,
            'I' => CigarKind::Ins,
            'D' => CigarKind::Del,
            'N' => CigarKind::Skip,
            'S' => CigarKind::SoftClip,
            'H' => CigarKind::HardClip,
            'P' => CigarKind::Pad,
            '=' => CigarKind::Equal,
            'X' => CigarKind::Diff,
            _ => return None,
        })
    }

    /// `true` if the op consumes a base from the query (`SEQ`).
    pub fn consumes_query(self) -> bool {
        matches!(
            self,
            CigarKind::Match
                | CigarKind::Ins
                | CigarKind::SoftClip
                | CigarKind::Equal
                | CigarKind::Diff
        )
    }

    /// `true` if the op consumes a base from the reference.
    pub fn consumes_ref(self) -> bool {
        matches!(
            self,
            CigarKind::Match
                | CigarKind::Del
                | CigarKind::Skip
                | CigarKind::Equal
                | CigarKind::Diff
        )
    }
}

/// A run-length CIGAR operation: a count and an op kind.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CigarOp {
    /// The run length (always positive in a well-formed CIGAR).
    pub len: u32,
    /// The operation kind.
    pub kind: CigarKind,
}

/// Per-CIGAR-op run-length cap. Round-6 fix: a 100-byte SAM with
/// `CIGAR=4294967295M` would otherwise drive `build_pileup` to
/// allocate ~4.3 billion `BTreeMap` entries — a parse-time DoS that
/// outsizes the input by ten orders of magnitude. 1 M bases is more
/// than enough for any real alignment op (the longest credible run
/// is an `N` skip across a multi-megabase intron, but even those
/// don't exceed ~5 Mb; we choose 1 M as a defensible upper bound
/// that catches the attack while leaving a generous safety margin
/// for legitimate spliced reads).
pub const MAX_CIGAR_OP_LEN: u32 = 1_000_000;

/// Sister cap to [`MAX_CIGAR_OP_LEN`]: the maximum number of *ops* a
/// single CIGAR string may contain. Without this cap, a tiny SAM
/// input like `"1I1I1I…"` (millions of single-op tokens) parses in
/// linear time but produces a `Vec<CigarOp>` with millions of entries,
/// blowing memory and feeding the same amplification-DoS class the
/// per-op cap defuses for length. 1 million ops is well past any
/// real read (whole-genome long reads cap at ~50k ops); we choose 1M
/// as a defensible safety margin that catches the attack.
pub const MAX_CIGAR_OPS: usize = 1_000_000;

/// A parsed CIGAR string — an ordered list of run-length operations.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Cigar {
    /// The operations in 5′→3′ query order.
    pub ops: Vec<CigarOp>,
}

impl Cigar {
    /// An empty CIGAR (the SAM `*` placeholder).
    pub fn new() -> Self {
        Cigar::default()
    }

    /// `true` if the CIGAR has no operations (an unmapped read).
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// Parses a CIGAR string such as `"8M1I4M1D3M"`. The single
    /// character `*` parses to an empty CIGAR.
    ///
    /// Each run length is bounded by [`MAX_CIGAR_OP_LEN`]; an op
    /// past the cap is rejected at parse time to defuse the
    /// amplification-DoS class (a tiny SAM input that drives a
    /// downstream pileup builder to allocate billions of entries).
    pub fn parse(s: &str) -> Result<Self> {
        if s == "*" || s.is_empty() {
            return Ok(Cigar::new());
        }
        let mut ops = Vec::new();
        let mut num = 0u32;
        let mut saw_digit = false;
        for c in s.chars() {
            if let Some(d) = c.to_digit(10) {
                num = num
                    .checked_mul(10)
                    .and_then(|n| n.checked_add(d))
                    .ok_or_else(|| GenomicsError::parse("sam", "CIGAR run length overflow"))?;
                saw_digit = true;
            } else {
                let kind = CigarKind::from_code(c).ok_or_else(|| {
                    GenomicsError::parse("sam", format!("unknown CIGAR op `{c}`"))
                })?;
                if !saw_digit {
                    return Err(GenomicsError::parse("sam", "CIGAR op without a length"));
                }
                if num > MAX_CIGAR_OP_LEN {
                    return Err(GenomicsError::parse(
                        "sam",
                        format!(
                            "CIGAR op length {num} exceeds the {MAX_CIGAR_OP_LEN}-base cap \
                             (amplification-DoS guard)"
                        ),
                    ));
                }
                if ops.len() >= MAX_CIGAR_OPS {
                    return Err(GenomicsError::parse(
                        "sam",
                        format!(
                            "CIGAR contains more than {MAX_CIGAR_OPS} ops \
                             (amplification-DoS guard — sister cap to MAX_CIGAR_OP_LEN)"
                        ),
                    ));
                }
                ops.push(CigarOp { len: num, kind });
                num = 0;
                saw_digit = false;
            }
        }
        if saw_digit {
            return Err(GenomicsError::parse("sam", "CIGAR ends with a trailing length"));
        }
        Ok(Cigar { ops })
    }

    /// Total query (`SEQ`) length the CIGAR consumes — soft-clips
    /// included, hard-clips excluded.
    pub fn query_len(&self) -> usize {
        self.ops
            .iter()
            .filter(|o| o.kind.consumes_query())
            .map(|o| o.len as usize)
            .sum()
    }

    /// Total reference length the alignment spans.
    pub fn ref_len(&self) -> usize {
        self.ops
            .iter()
            .filter(|o| o.kind.consumes_ref())
            .map(|o| o.len as usize)
            .sum()
    }

    /// Renders the CIGAR back to its SAM string (`*` when empty).
    pub fn to_sam_string(&self) -> String {
        if self.ops.is_empty() {
            return "*".to_string();
        }
        let mut s = String::new();
        for op in &self.ops {
            s.push_str(&op.len.to_string());
            s.push(op.kind.code());
        }
        s
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_sam_string())
    }
}

// --------------------------------------------------------------------
// SAM flag bits
// --------------------------------------------------------------------

/// The twelve SAM `FLAG` bits, exposed as a typed wrapper around the
/// raw `u16`.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct SamFlags(pub u16);

impl SamFlags {
    /// `0x1` — the read is paired in sequencing.
    pub const PAIRED: u16 = 0x1;
    /// `0x2` — the read is mapped in a proper pair.
    pub const PROPER_PAIR: u16 = 0x2;
    /// `0x4` — the read itself is unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// `0x8` — the mate is unmapped.
    pub const MATE_UNMAPPED: u16 = 0x8;
    /// `0x10` — the read maps to the reverse strand.
    pub const REVERSE: u16 = 0x10;
    /// `0x20` — the mate maps to the reverse strand.
    pub const MATE_REVERSE: u16 = 0x20;
    /// `0x40` — the read is the first in the pair.
    pub const FIRST: u16 = 0x40;
    /// `0x80` — the read is the second in the pair.
    pub const LAST: u16 = 0x80;
    /// `0x100` — a secondary (multi-mapping) alignment.
    pub const SECONDARY: u16 = 0x100;
    /// `0x200` — the read failed QC.
    pub const QC_FAIL: u16 = 0x200;
    /// `0x400` — the read is a PCR / optical duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// `0x800` — a supplementary (chimeric) alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;

    /// `true` when the given bit mask is set.
    pub fn has(self, bit: u16) -> bool {
        self.0 & bit != 0
    }

    /// Sets or clears the given bit mask.
    pub fn set(&mut self, bit: u16, on: bool) {
        if on {
            self.0 |= bit;
        } else {
            self.0 &= !bit;
        }
    }

    /// Shorthand: `true` when the [`UNMAPPED`](Self::UNMAPPED) bit is set.
    pub fn is_unmapped(self) -> bool {
        self.has(Self::UNMAPPED)
    }

    /// Shorthand: `true` when the [`REVERSE`](Self::REVERSE) bit is set.
    pub fn is_reverse(self) -> bool {
        self.has(Self::REVERSE)
    }

    /// Shorthand: `true` when the [`PAIRED`](Self::PAIRED) bit is set.
    pub fn is_paired(self) -> bool {
        self.has(Self::PAIRED)
    }

    /// Shorthand: `true` when the [`DUPLICATE`](Self::DUPLICATE) bit is set.
    pub fn is_duplicate(self) -> bool {
        self.has(Self::DUPLICATE)
    }

    /// `true` for a secondary or supplementary alignment — the records
    /// most callers want to skip when counting "primary" alignments.
    pub fn is_secondary_or_supplementary(self) -> bool {
        self.has(Self::SECONDARY) || self.has(Self::SUPPLEMENTARY)
    }
}

// --------------------------------------------------------------------
// Optional tags
// --------------------------------------------------------------------

/// The typed value of a SAM optional `TAG:TYPE:VALUE` field.
#[derive(Clone, Debug, PartialEq)]
pub enum SamTagValue {
    /// `A` — a printable character.
    Char(char),
    /// `i` — a signed 32-bit integer.
    Int(i64),
    /// `f` — a single-precision float.
    Float(f64),
    /// `Z` — a printable string.
    Str(String),
    /// `H` — a byte array rendered as a hex string (kept verbatim).
    Hex(String),
    /// `B` — a numeric array; stored as its raw `subtype,v1,v2,…` text.
    Array(String),
}

impl SamTagValue {
    /// The single-character SAM type code for this value.
    pub fn type_code(&self) -> char {
        match self {
            SamTagValue::Char(_) => 'A',
            SamTagValue::Int(_) => 'i',
            SamTagValue::Float(_) => 'f',
            SamTagValue::Str(_) => 'Z',
            SamTagValue::Hex(_) => 'H',
            SamTagValue::Array(_) => 'B',
        }
    }

    /// Renders the value to its SAM `VALUE` text (no `TAG:TYPE:` prefix).
    pub fn to_value_string(&self) -> String {
        match self {
            SamTagValue::Char(c) => c.to_string(),
            SamTagValue::Int(i) => i.to_string(),
            SamTagValue::Float(x) => format!("{x}"),
            SamTagValue::Str(s) => s.clone(),
            SamTagValue::Hex(s) => s.clone(),
            SamTagValue::Array(s) => s.clone(),
        }
    }
}

/// One optional SAM field: a two-character tag plus a typed value.
#[derive(Clone, Debug, PartialEq)]
pub struct SamTag {
    /// The two-character tag name (`NM`, `MD`, `AS`, `RG`, …).
    pub tag: String,
    /// The typed value.
    pub value: SamTagValue,
}

impl SamTag {
    /// Parses one `TAG:TYPE:VALUE` token.
    pub fn parse(token: &str) -> Result<Self> {
        let mut parts = token.splitn(3, ':');
        let tag = parts
            .next()
            .filter(|t| t.len() == 2)
            .ok_or_else(|| GenomicsError::parse("sam", "optional tag must be `XX:T:V`"))?
            .to_string();
        let ty = parts
            .next()
            .ok_or_else(|| GenomicsError::parse("sam", "optional tag missing type"))?;
        let raw = parts
            .next()
            .ok_or_else(|| GenomicsError::parse("sam", "optional tag missing value"))?;
        let value = match ty {
            "A" => SamTagValue::Char(raw.chars().next().ok_or_else(|| {
                GenomicsError::parse("sam", "empty `A`-type tag value")
            })?),
            "i" => SamTagValue::Int(
                raw.parse::<i64>()
                    .map_err(|_| GenomicsError::parse("sam", "bad `i`-type tag value"))?,
            ),
            "f" => SamTagValue::Float(
                raw.parse::<f64>()
                    .map_err(|_| GenomicsError::parse("sam", "bad `f`-type tag value"))?,
            ),
            "Z" => SamTagValue::Str(raw.to_string()),
            "H" => SamTagValue::Hex(raw.to_string()),
            "B" => SamTagValue::Array(raw.to_string()),
            other => {
                return Err(GenomicsError::parse(
                    "sam",
                    format!("unknown optional-tag type `{other}`"),
                ))
            }
        };
        Ok(SamTag { tag, value })
    }

    /// Renders the tag back to its `TAG:TYPE:VALUE` text.
    pub fn to_sam_string(&self) -> String {
        format!(
            "{}:{}:{}",
            self.tag,
            self.value.type_code(),
            self.value.to_value_string()
        )
    }
}

// --------------------------------------------------------------------
// SAM record
// --------------------------------------------------------------------

/// One SAM alignment record — the eleven mandatory fields plus any
/// optional tags.
#[derive(Clone, Debug, PartialEq)]
pub struct SamRecord {
    /// `QNAME` — the read / template name.
    pub qname: String,
    /// `FLAG` — the bitwise flags.
    pub flags: SamFlags,
    /// `RNAME` — the reference name (`*` for unmapped).
    pub rname: String,
    /// `POS` — 1-based leftmost mapping position (`0` for unmapped).
    pub pos: i64,
    /// `MAPQ` — mapping quality (`255` = unavailable).
    pub mapq: u8,
    /// `CIGAR` — the parsed alignment CIGAR.
    pub cigar: Cigar,
    /// `RNEXT` — the mate's reference name (`*` / `=`).
    pub rnext: String,
    /// `PNEXT` — the mate's 1-based position (`0` if unavailable).
    pub pnext: i64,
    /// `TLEN` — observed template length (signed; `0` if unavailable).
    pub tlen: i64,
    /// `SEQ` — the read bases (empty for the SAM `*` placeholder).
    pub seq: String,
    /// `QUAL` — the per-base ASCII Phred+33 qualities (empty for `*`).
    pub qual: String,
    /// The optional `TAG:TYPE:VALUE` fields, in file order.
    pub tags: Vec<SamTag>,
}

impl SamRecord {
    /// An empty unmapped record carrying only a name.
    pub fn unmapped(qname: impl Into<String>) -> Self {
        SamRecord {
            qname: qname.into(),
            flags: SamFlags(SamFlags::UNMAPPED),
            rname: "*".to_string(),
            pos: 0,
            mapq: 0,
            cigar: Cigar::new(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq: String::new(),
            qual: String::new(),
            tags: Vec::new(),
        }
    }

    /// `true` when the record's [`UNMAPPED`](SamFlags::UNMAPPED) bit is set.
    pub fn is_unmapped(&self) -> bool {
        self.flags.is_unmapped()
    }

    /// 1-based inclusive reference end position of the alignment, or
    /// `None` for an unmapped record. `end = pos + ref_len - 1`.
    pub fn ref_end(&self) -> Option<i64> {
        if self.is_unmapped() || self.pos == 0 {
            return None;
        }
        let span = self.cigar.ref_len() as i64;
        Some(self.pos + span.max(1) - 1)
    }

    /// Looks up an optional tag by its two-character name.
    pub fn tag(&self, name: &str) -> Option<&SamTag> {
        self.tags.iter().find(|t| t.tag == name)
    }

    /// Validates that `SEQ`, `QUAL` and the CIGAR query length agree —
    /// the single most common source of a corrupt SAM record.
    pub fn validate(&self) -> Result<()> {
        let seq_set = self.seq != "*" && !self.seq.is_empty();
        let qual_set = self.qual != "*" && !self.qual.is_empty();
        if seq_set && qual_set && self.seq.len() != self.qual.len() {
            return Err(GenomicsError::invalid_record(
                "sam",
                format!(
                    "SEQ length {} != QUAL length {}",
                    self.seq.len(),
                    self.qual.len()
                ),
            ));
        }
        if seq_set && !self.cigar.is_empty() {
            let qlen = self.cigar.query_len();
            if qlen != self.seq.len() {
                return Err(GenomicsError::invalid_record(
                    "sam",
                    format!("CIGAR query length {qlen} != SEQ length {}", self.seq.len()),
                ));
            }
        }
        Ok(())
    }

    /// Parses one tab-delimited SAM body line.
    pub fn parse_line(line: &str) -> Result<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            return Err(GenomicsError::parse(
                "sam",
                format!("alignment line has {} fields, need >= 11", f.len()),
            ));
        }
        let flags = SamFlags(
            f[1]
                .parse::<u16>()
                .map_err(|_| GenomicsError::parse("sam", "bad FLAG"))?,
        );
        let pos = f[3]
            .parse::<i64>()
            .map_err(|_| GenomicsError::parse("sam", "bad POS"))?;
        let mapq = f[4]
            .parse::<u8>()
            .map_err(|_| GenomicsError::parse("sam", "bad MAPQ"))?;
        let cigar = Cigar::parse(f[5])?;
        let pnext = f[7]
            .parse::<i64>()
            .map_err(|_| GenomicsError::parse("sam", "bad PNEXT"))?;
        let tlen = f[8]
            .parse::<i64>()
            .map_err(|_| GenomicsError::parse("sam", "bad TLEN"))?;
        let seq = if f[9] == "*" { String::new() } else { f[9].to_string() };
        let qual = if f[10] == "*" { String::new() } else { f[10].to_string() };
        let mut tags = Vec::new();
        for tok in &f[11..] {
            if tok.is_empty() {
                continue;
            }
            tags.push(SamTag::parse(tok)?);
        }
        let rec = SamRecord {
            qname: f[0].to_string(),
            flags,
            rname: f[2].to_string(),
            pos,
            mapq,
            cigar,
            rnext: f[6].to_string(),
            pnext,
            tlen,
            seq,
            qual,
            tags,
        };
        rec.validate()?;
        Ok(rec)
    }

    /// Renders the record back to one tab-delimited SAM line.
    pub fn to_sam_line(&self) -> String {
        let seq = if self.seq.is_empty() { "*" } else { &self.seq };
        let qual = if self.qual.is_empty() { "*" } else { &self.qual };
        let mut s = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname,
            self.flags.0,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar.to_sam_string(),
            self.rnext,
            self.pnext,
            self.tlen,
            seq,
            qual,
        );
        for tag in &self.tags {
            s.push('\t');
            s.push_str(&tag.to_sam_string());
        }
        s
    }
}

// --------------------------------------------------------------------
// SAM header
// --------------------------------------------------------------------

/// A `@SQ` reference-sequence header entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RefSeq {
    /// `SN` — the reference name.
    pub name: String,
    /// `LN` — the reference length in bases.
    pub length: u64,
}

/// The typed SAM header: the `@HD` / `@SQ` / `@RG` / `@PG` lines.
///
/// `@SQ` records are typed into [`RefSeq`]; the remaining lines are
/// kept as their raw text so nothing is lost on a round-trip.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SamHeader {
    /// The `@HD` version string, if any (the `VN` value).
    pub version: Option<String>,
    /// The `@HD` sort order, if any (the `SO` value).
    pub sort_order: Option<String>,
    /// The `@SQ` reference-sequence dictionary, in file order.
    pub references: Vec<RefSeq>,
    /// Raw `@RG` read-group lines, verbatim.
    pub read_groups: Vec<String>,
    /// Raw `@PG` program lines, verbatim.
    pub programs: Vec<String>,
    /// Raw `@CO` free-text comment lines, verbatim.
    pub comments: Vec<String>,
}

impl SamHeader {
    /// An empty header.
    pub fn new() -> Self {
        SamHeader::default()
    }

    /// Looks up a reference by name, returning its 0-based index.
    pub fn ref_index(&self, name: &str) -> Option<usize> {
        self.references.iter().position(|r| r.name == name)
    }

    /// Parses one `@`-prefixed header line into the accumulating header.
    fn parse_into(&mut self, line: &str) -> Result<()> {
        let mut fields = line.split('\t');
        let kind = fields
            .next()
            .ok_or_else(|| GenomicsError::parse("sam", "empty header line"))?;
        // Collect `KEY:VALUE` tag fields.
        let mut kv: Vec<(&str, &str)> = Vec::new();
        for f in fields {
            if let Some((k, v)) = f.split_once(':') {
                kv.push((k, v));
            }
        }
        match kind {
            "@HD" => {
                for (k, v) in kv {
                    match k {
                        "VN" => self.version = Some(v.to_string()),
                        "SO" => self.sort_order = Some(v.to_string()),
                        _ => {}
                    }
                }
            }
            "@SQ" => {
                let mut name = None;
                let mut length = None;
                for (k, v) in kv {
                    match k {
                        "SN" => name = Some(v.to_string()),
                        "LN" => {
                            length = Some(v.parse::<u64>().map_err(|_| {
                                GenomicsError::parse("sam", "bad @SQ LN")
                            })?)
                        }
                        _ => {}
                    }
                }
                let name = name
                    .ok_or_else(|| GenomicsError::parse("sam", "@SQ missing SN"))?;
                let length = length
                    .ok_or_else(|| GenomicsError::parse("sam", "@SQ missing LN"))?;
                self.references.push(RefSeq { name, length });
            }
            "@RG" => self.read_groups.push(line.to_string()),
            "@PG" => self.programs.push(line.to_string()),
            "@CO" => self.comments.push(line.to_string()),
            other => {
                return Err(GenomicsError::parse(
                    "sam",
                    format!("unknown header line `{other}`"),
                ))
            }
        }
        Ok(())
    }

    /// Renders the header back to its `@`-prefixed lines (newline-joined,
    /// no trailing newline).
    pub fn to_sam_string(&self) -> String {
        let mut out: Vec<String> = Vec::new();
        if self.version.is_some() || self.sort_order.is_some() {
            let mut hd = "@HD".to_string();
            if let Some(v) = &self.version {
                hd.push_str(&format!("\tVN:{v}"));
            }
            if let Some(so) = &self.sort_order {
                hd.push_str(&format!("\tSO:{so}"));
            }
            out.push(hd);
        }
        for r in &self.references {
            out.push(format!("@SQ\tSN:{}\tLN:{}", r.name, r.length));
        }
        out.extend(self.read_groups.iter().cloned());
        out.extend(self.programs.iter().cloned());
        out.extend(self.comments.iter().cloned());
        out.join("\n")
    }
}

// --------------------------------------------------------------------
// SAM file
// --------------------------------------------------------------------

/// A parsed SAM file — its header and every alignment record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SamFile {
    /// The typed header.
    pub header: SamHeader,
    /// Every alignment record, in file order.
    pub records: Vec<SamRecord>,
}

impl SamFile {
    /// An empty SAM file.
    pub fn new() -> Self {
        SamFile::default()
    }

    /// Parses an entire SAM document.
    pub fn parse(text: &str) -> Result<Self> {
        let mut header = SamHeader::new();
        let mut records = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            if let Some(stripped) = line.strip_prefix('@') {
                // `strip_prefix` removed the `@`; rebuild for the parser.
                let _ = stripped;
                header.parse_into(line)?;
            } else {
                records.push(SamRecord::parse_line(line)?);
            }
        }
        Ok(SamFile { header, records })
    }

    /// Renders the whole file back to SAM text (with a trailing
    /// newline).
    pub fn to_sam_string(&self) -> String {
        let mut s = self.header.to_sam_string();
        if !s.is_empty() {
            s.push('\n');
        }
        for rec in &self.records {
            s.push_str(&rec.to_sam_line());
            s.push('\n');
        }
        s
    }

    /// Count of mapped (non-`UNMAPPED`) records.
    pub fn mapped_count(&self) -> usize {
        self.records.iter().filter(|r| !r.is_unmapped()).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
@SQ\tSN:chr2\tLN:500\n\
@PG\tID:bwa\tPN:bwa\n\
r1\t0\tchr1\t10\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\tAS:i:8\n\
r2\t16\tchr1\t20\t60\t4M1I3M\t*\t0\t0\tACGTACGT\t########\n\
r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n";

    #[test]
    fn cigar_parses_and_round_trips() {
        let c = Cigar::parse("8M1I4M1D3M").unwrap();
        assert_eq!(c.ops.len(), 5);
        assert_eq!(c.query_len(), 8 + 1 + 4 + 3);
        assert_eq!(c.ref_len(), 8 + 4 + 1 + 3);
        assert_eq!(c.to_sam_string(), "8M1I4M1D3M");
        assert!(Cigar::parse("*").unwrap().is_empty());
    }

    #[test]
    fn cigar_rejects_garbage() {
        assert!(Cigar::parse("8Z").is_err());
        assert!(Cigar::parse("M").is_err());
        assert!(Cigar::parse("8M3").is_err());
    }

    #[test]
    fn cigar_rejects_op_length_past_max_cap() {
        // Round-6 RED→GREEN: a single CIGAR op claiming u32::MAX bases
        // would let `build_pileup` insert ~4.3 B BTreeMap entries from
        // a 100-byte SAM input. Reject at parse time.
        let err = Cigar::parse("4294967295M").unwrap_err();
        match err {
            GenomicsError::Parse { format, detail } => {
                assert_eq!(format, "sam");
                assert!(detail.contains("exceeds"), "msg: {detail}");
                assert!(
                    detail.contains(&MAX_CIGAR_OP_LEN.to_string()),
                    "msg: {detail}"
                );
            }
            other => panic!("expected Parse, got {other:?}"),
        }
        // The cap edge: exactly MAX accepts, MAX+1 rejects.
        assert!(Cigar::parse(&format!("{MAX_CIGAR_OP_LEN}M")).is_ok());
        let over = MAX_CIGAR_OP_LEN + 1;
        assert!(Cigar::parse(&format!("{over}M")).is_err());
    }

    #[test]
    fn cigar_rejects_op_count_past_max_cap() {
        // Round-8 RED→GREEN: a CIGAR with millions of tiny ops
        // (`"1I1I1I…"`) parses in linear time but produces a Vec with
        // millions of entries — amplification-DoS sister to the
        // per-op-length cap. The MAX_CIGAR_OPS cap defuses it.
        //
        // We construct a CIGAR with MAX_CIGAR_OPS+1 ops to confirm the
        // cap edge. Building the string is O(MAX_CIGAR_OPS) work
        // (~5 MiB), but the test only runs once per CI.
        let mut s = String::with_capacity(MAX_CIGAR_OPS * 2 + 4);
        for _ in 0..=MAX_CIGAR_OPS {
            s.push_str("1I");
        }
        let err = Cigar::parse(&s).unwrap_err();
        match err {
            GenomicsError::Parse { format, detail } => {
                assert_eq!(format, "sam");
                assert!(
                    detail.contains("more than") && detail.contains(&MAX_CIGAR_OPS.to_string()),
                    "msg: {detail}"
                );
            }
            other => panic!("expected Parse, got {other:?}"),
        }
    }

    #[test]
    fn flags_round_trip() {
        let mut f = SamFlags::default();
        f.set(SamFlags::PAIRED, true);
        f.set(SamFlags::REVERSE, true);
        assert!(f.is_paired());
        assert!(f.is_reverse());
        f.set(SamFlags::REVERSE, false);
        assert!(!f.is_reverse());
    }

    #[test]
    fn tag_parses_each_type() {
        assert_eq!(
            SamTag::parse("NM:i:3").unwrap().value,
            SamTagValue::Int(3)
        );
        assert_eq!(
            SamTag::parse("RG:Z:group1").unwrap().value,
            SamTagValue::Str("group1".to_string())
        );
        assert!(matches!(
            SamTag::parse("XF:f:1.5").unwrap().value,
            SamTagValue::Float(_)
        ));
        assert!(SamTag::parse("BAD").is_err());
    }

    #[test]
    fn parse_full_file() {
        let f = SamFile::parse(SAMPLE).unwrap();
        assert_eq!(f.header.version.as_deref(), Some("1.6"));
        assert_eq!(f.header.references.len(), 2);
        assert_eq!(f.header.ref_index("chr2"), Some(1));
        assert_eq!(f.records.len(), 3);
        assert_eq!(f.mapped_count(), 2);
        assert!(f.records[2].is_unmapped());
        assert_eq!(f.records[0].tag("NM").unwrap().value, SamTagValue::Int(0));
    }

    #[test]
    fn ref_end_uses_cigar_span() {
        let f = SamFile::parse(SAMPLE).unwrap();
        // r1: pos 10, 8M -> end 17.
        assert_eq!(f.records[0].ref_end(), Some(17));
        // r2: pos 20, 4M1I3M -> ref span 7 -> end 26.
        assert_eq!(f.records[1].ref_end(), Some(26));
        // r3 unmapped.
        assert_eq!(f.records[2].ref_end(), None);
    }

    #[test]
    fn record_validates_seq_qual_cigar() {
        let mut rec = SamRecord::unmapped("x");
        rec.seq = "ACGT".to_string();
        rec.qual = "III".to_string();
        assert!(rec.validate().is_err());
        rec.qual = "IIII".to_string();
        rec.cigar = Cigar::parse("3M").unwrap();
        assert!(rec.validate().is_err());
        rec.cigar = Cigar::parse("4M").unwrap();
        assert!(rec.validate().is_ok());
    }

    #[test]
    fn round_trip_preserves_records() {
        let f = SamFile::parse(SAMPLE).unwrap();
        let text = f.to_sam_string();
        let f2 = SamFile::parse(&text).unwrap();
        assert_eq!(f, f2);
    }
}
