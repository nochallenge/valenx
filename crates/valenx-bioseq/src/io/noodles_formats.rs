//! Native VCF (variant-call) and SAM (alignment) parsing.
//!
//! This module rounds out the crate's file-format coverage. The
//! in-house readers in [`fasta`](super::fasta),
//! [`fastq`](super::fastq), [`genbank`](super::genbank) and
//! [`embl`](super::embl) cover the *sequence* formats; the two formats
//! parsed here — VCF and SAM — describe *variants against a reference*
//! and *read alignments*, and are delegated to the pure-Rust
//! [`noodles`](https://crates.io/crates/noodles) family
//! (`noodles-vcf`, `noodles-sam`). That keeps the whole stack free of
//! the C `htslib` dependency that the reference tools (`bcftools`,
//! `samtools`) build on, with no external process and no FFI.
//!
//! The noodles record models are rich; this module deliberately
//! projects each record down to a small, owned, `serde`-friendly
//! struct holding the fields most analyses actually use. Parsing the
//! full INFO/FORMAT/genotype and CIGAR/tag detail is left to callers
//! that reach for noodles directly.
//!
//! Only the text encodings are handled here. The bgzf-compressed
//! binary siblings (BAM, BCF) are out of scope for this entry point.
//!
//! # Examples
//!
//! ```
//! use valenx_bioseq::io::noodles_formats::read_vcf_str;
//!
//! // Columns are tab-separated (built with `join` to keep the tabs
//! // unambiguous in this doc comment).
//! let header_line = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"].join("\t");
//! let data_line = ["chr1", "100", "rs1", "A", "G", ".", ".", "."].join("\t");
//! let vcf = format!(
//!     "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n{header_line}\n{data_line}\n"
//! );
//! let variants = read_vcf_str(&vcf).unwrap();
//! assert_eq!(variants.len(), 1);
//! assert_eq!(variants[0].chrom, "chr1");
//! assert_eq!(variants[0].pos, 100);
//! assert_eq!(variants[0].alt, vec!["G".to_string()]);
//! ```

use std::io::BufRead;
use std::path::Path;

use crate::error::{BioseqError, Result};

// --- VCF --------------------------------------------------------------

/// One variant record projected from a VCF data line.
///
/// Positions are 1-based, matching the VCF specification (and the
/// reference tools). Coordinates are stored as `usize`; the optional
/// fields mirror the VCF `.` "missing" convention with [`Option`] /
/// empty collections.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct VariantRecord {
    /// `CHROM` — the reference-sequence (chromosome/contig) name.
    pub chrom: String,
    /// `POS` — 1-based reference position of the variant start.
    pub pos: usize,
    /// `ID` — the list of variant identifiers (e.g. dbSNP `rs` ids);
    /// empty when the column was `.`.
    pub ids: Vec<String>,
    /// `REF` — the reference allele bases.
    pub reference: String,
    /// `ALT` — the alternate alleles; empty when the column was `.`.
    pub alt: Vec<String>,
}

/// Parses an entire VCF document held in memory into [`VariantRecord`]s.
///
/// The `##`/`#CHROM` header is parsed (and required) but only the data
/// lines are returned. Returns [`BioseqError::Parse`] if the text is
/// not valid VCF.
pub fn read_vcf_str(text: &str) -> Result<Vec<VariantRecord>> {
    read_vcf(text.as_bytes())
}

/// Parses VCF from any buffered byte source into [`VariantRecord`]s.
///
/// Accepts anything implementing [`BufRead`] (a file, a `&[u8]`, a
/// network stream). For a path, see [`read_vcf_path`].
pub fn read_vcf<R: BufRead>(reader: R) -> Result<Vec<VariantRecord>> {
    use noodles_vcf as vcf;

    let mut rdr = vcf::io::Reader::new(reader);
    let header = rdr
        .read_header()
        .map_err(|e| BioseqError::parse("vcf", format!("invalid VCF header: {e}")))?;

    let mut out = Vec::new();
    for result in rdr.record_bufs(&header) {
        let rec =
            result.map_err(|e| BioseqError::parse("vcf", format!("invalid VCF record: {e}")))?;

        let pos = rec.variant_start().map(|p| p.get()).unwrap_or(0);
        let ids = rec.ids().as_ref().iter().cloned().collect();
        let alt = rec.alternate_bases().as_ref().to_vec();

        out.push(VariantRecord {
            chrom: rec.reference_sequence_name().to_string(),
            pos,
            ids,
            reference: rec.reference_bases().to_string(),
            alt,
        });
    }
    Ok(out)
}

/// Reads and parses a VCF file from `path` into [`VariantRecord`]s.
pub fn read_vcf_path<P: AsRef<Path>>(path: P) -> Result<Vec<VariantRecord>> {
    let file = std::fs::File::open(path.as_ref())
        .map_err(|e| BioseqError::parse("vcf", format!("cannot open VCF: {e}")))?;
    read_vcf(std::io::BufReader::new(file))
}

// --- SAM --------------------------------------------------------------

/// One read alignment projected from a SAM data line.
///
/// The reference name is resolved from the alignment's reference id
/// against the SAM header's `@SQ` dictionary; unmapped reads (no
/// reference) leave it [`None`]. The start position is 1-based.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct AlignmentRecord {
    /// `QNAME` — the read (query template) name; [`None`] when `*`.
    pub name: Option<String>,
    /// `FLAG` — the raw SAM bitwise flag.
    pub flags: u16,
    /// Reference-sequence name the read aligns to (`RNAME`), resolved
    /// via the header; [`None`] for unmapped reads.
    pub reference_name: Option<String>,
    /// `POS` — 1-based leftmost alignment position; [`None`] when
    /// unmapped.
    pub pos: Option<usize>,
    /// `SEQ` — the read bases (`*` yields an empty string).
    pub sequence: String,
}

impl AlignmentRecord {
    /// `true` if the `UNMAPPED` (0x4) flag bit is set.
    pub fn is_unmapped(&self) -> bool {
        self.flags & 0x4 != 0
    }
}

/// Parses an entire SAM document held in memory into
/// [`AlignmentRecord`]s.
pub fn read_sam_str(text: &str) -> Result<Vec<AlignmentRecord>> {
    read_sam(text.as_bytes())
}

/// Parses SAM from any buffered byte source into [`AlignmentRecord`]s.
pub fn read_sam<R: BufRead>(reader: R) -> Result<Vec<AlignmentRecord>> {
    use noodles_sam as sam;

    let mut rdr = sam::io::Reader::new(reader);
    let header = rdr
        .read_header()
        .map_err(|e| BioseqError::parse("sam", format!("invalid SAM header: {e}")))?;
    let refs = header.reference_sequences();

    let mut out = Vec::new();
    for result in rdr.record_bufs(&header) {
        let rec =
            result.map_err(|e| BioseqError::parse("sam", format!("invalid SAM record: {e}")))?;

        let name = rec.name().map(|n| n.to_string());
        let reference_name = rec
            .reference_sequence_id()
            .and_then(|id| refs.get_index(id).map(|(name, _meta)| name.to_string()));
        let pos = rec.alignment_start().map(|p| p.get());
        let sequence = String::from_utf8_lossy(rec.sequence().as_ref()).into_owned();

        out.push(AlignmentRecord {
            name,
            flags: u16::from(rec.flags()),
            reference_name,
            pos,
            sequence,
        });
    }
    Ok(out)
}

/// Reads and parses a SAM file from `path` into [`AlignmentRecord`]s.
pub fn read_sam_path<P: AsRef<Path>>(path: P) -> Result<Vec<AlignmentRecord>> {
    let file = std::fs::File::open(path.as_ref())
        .map_err(|e| BioseqError::parse("sam", format!("cannot open SAM: {e}")))?;
    read_sam(std::io::BufReader::new(file))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Writes `content` to a uniquely-named temp file and returns its
    /// path. Exercises the path-based readers end to end.
    fn temp_file(stem: &str, ext: &str, content: &str) -> std::path::PathBuf {
        let pid = std::process::id();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let mut path = std::env::temp_dir();
        path.push(format!("valenx_bioseq_{stem}_{pid}_{nanos}.{ext}"));
        let mut f = std::fs::File::create(&path).expect("create temp file");
        f.write_all(content.as_bytes()).expect("write temp file");
        path
    }

    const SAMPLE_VCF: &str = "\
##fileformat=VCFv4.3
##contig=<ID=chr1>
##contig=<ID=chr2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs1\tA\tG\t50\tPASS\t.
chr1\t250\t.\tC\tT,A\t.\t.\t.
chr2\t17\trs99\tGAT\tG\t30\tPASS\t.
";

    #[test]
    fn vcf_parses_records_ids_and_alleles() {
        let variants = read_vcf_str(SAMPLE_VCF).expect("parse VCF");
        assert_eq!(variants.len(), 3, "record count");

        // First record: a simple SNV with a dbSNP id.
        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].pos, 100);
        assert_eq!(variants[0].ids, vec!["rs1".to_string()]);
        assert_eq!(variants[0].reference, "A");
        assert_eq!(variants[0].alt, vec!["G".to_string()]);

        // Second record: missing id (`.`) and two alternate alleles.
        assert!(variants[1].ids.is_empty(), "missing ID -> empty");
        assert_eq!(
            variants[1].alt,
            vec!["T".to_string(), "A".to_string()],
            "multi-allelic ALT"
        );

        // Third record: a deletion on a second contig.
        assert_eq!(variants[2].chrom, "chr2");
        assert_eq!(variants[2].pos, 17);
        assert_eq!(variants[2].reference, "GAT");
        assert_eq!(variants[2].alt, vec!["G".to_string()]);
    }

    #[test]
    fn vcf_round_trips_via_temp_file() {
        let path = temp_file("variants", "vcf", SAMPLE_VCF);
        let variants = read_vcf_path(&path).expect("parse VCF from file");
        let _ = std::fs::remove_file(&path);
        assert_eq!(variants.len(), 3);
        assert_eq!(variants[0].chrom, "chr1");
    }

    const SAMPLE_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:500
r001\t0\tchr1\t7\t30\t8M\t*\t0\t0\tTTAGATAA\t*
r002\t16\tchr2\t9\t30\t5M\t*\t0\t0\tAACGT\t*
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\t*
";

    #[test]
    fn sam_parses_names_positions_and_sequences() {
        let alns = read_sam_str(SAMPLE_SAM).expect("parse SAM");
        assert_eq!(alns.len(), 3, "record count");

        assert_eq!(alns[0].name.as_deref(), Some("r001"));
        assert_eq!(alns[0].flags, 0);
        assert_eq!(alns[0].reference_name.as_deref(), Some("chr1"));
        assert_eq!(alns[0].pos, Some(7));
        assert_eq!(alns[0].sequence, "TTAGATAA");
        assert!(!alns[0].is_unmapped());

        assert_eq!(alns[1].reference_name.as_deref(), Some("chr2"));
        assert_eq!(alns[1].pos, Some(9));
        assert_eq!(alns[1].sequence, "AACGT");

        // The unmapped read: FLAG 0x4, no reference, no position.
        assert!(alns[2].is_unmapped(), "0x4 flag");
        assert_eq!(alns[2].reference_name, None);
        assert_eq!(alns[2].pos, None);
        assert_eq!(alns[2].sequence, "ACGTACGT");
    }

    #[test]
    fn sam_round_trips_via_temp_file() {
        let path = temp_file("aln", "sam", SAMPLE_SAM);
        let alns = read_sam_path(&path).expect("parse SAM from file");
        let _ = std::fs::remove_file(&path);
        assert_eq!(alns.len(), 3);
        assert_eq!(alns[1].name.as_deref(), Some("r002"));
    }
}
