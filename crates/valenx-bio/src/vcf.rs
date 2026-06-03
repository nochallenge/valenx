//! VCF (Variant Call Format) canonical type.
//!
//! Holds variant records produced by callers like `bcftools call`,
//! `gatk HaplotypeCaller`, or `run_deepvariant`. Mirrors the shape
//! of [`crate::Sequence`] / [`crate::Alignment`] / [`crate::FastqRecord`]:
//! a small struct that schema-validates only what the format
//! mandates, with downstream interpretation (INFO key parsing, genotype
//! decoding, variant typing) left to callers that care about specific
//! keys.
//!
//! BCF (binary VCF) and bgzf-compressed VCF are intentionally out of
//! scope — convert with `bcftools view` first. See
//! [`crate::format::vcf`] for the text reader.
//!
//! # Missing-value sentinels
//!
//! VCF uses `"."` as the missing-value sentinel for ID / QUAL / ALT /
//! FILTER. INFO is special: `"."` there means "no annotations" and
//! rides through as the literal `"."` string.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A parsed VCF document.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Vcf {
    /// `##` meta-information lines, verbatim.
    pub header: Vec<String>,
    /// Sample IDs parsed from the `#CHROM` line. Empty when no
    /// per-sample columns are present.
    pub samples: Vec<String>,
    pub records: Vec<VcfRecord>,
}

/// One VCF data row.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    /// Variant ID. `None` when the file had `"."`.
    pub id: Option<String>,
    /// Reference allele. Always populated.
    pub ref_allele: String,
    /// Alternate alleles. May be empty (the file had `"."`).
    pub alt: Vec<String>,
    /// Phred-scaled variant call quality. `None` when the file had
    /// `"."`.
    pub qual: Option<f64>,
    /// Filter status. Empty when the file had `"."` (= unfiltered).
    /// `["PASS"]` is the conventional "passed all filters" marker.
    pub filter: Vec<String>,
    /// INFO field, raw. Parsed downstream by callers that care about
    /// specific keys (DP, AF, MQ, etc.) — keeping it as a string here
    /// keeps the canonical type schema-agnostic.
    pub info: String,
    /// FORMAT field. `Some(...)` when column 9 was present in the
    /// row, regardless of whether per-sample columns followed —
    /// preserving the FORMAT info is informationally richer than
    /// suppressing it for the (rare, technically malformed) case
    /// where a row carries FORMAT but zero samples.
    pub format: Option<String>,
    /// Per-sample column strings, in the same order as `Vcf.samples`.
    pub samples: Vec<String>,
}

/// VCF parse errors.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum VcfError {
    /// Generic malformed-line catch-all.
    #[error("line {line}: {msg}")]
    Bad {
        /// 1-based line number of the offending line.
        line: usize,
        /// Short human-readable explanation.
        msg: String,
    },
    /// A data row's per-sample column count did not match the count
    /// promised by the `#CHROM` header line. A `#CHROM` line with
    /// fewer than 10 columns (no sample names) followed by data rows
    /// carrying sample columns is the classic header/data drift bug.
    /// Without this check the per-sample arrays would silently shift
    /// to wrong subjects.
    #[error("line {line}: sample count {got} differs from header sample count {expected}")]
    SampleCountMismatch {
        /// 1-based line number of the offending data row.
        line: usize,
        /// Number of sample columns the `#CHROM` header advertised
        /// (0 if no `#CHROM` was seen but a data row still carries
        /// sample columns — see [`Self::Bad`] guidance there).
        expected: usize,
        /// Number of sample columns present on the data row.
        got: usize,
    },
}

impl VcfRecord {
    /// True if the record has at least one alternate allele.
    pub fn has_alt(&self) -> bool {
        !self.alt.is_empty()
    }
    /// True if FILTER is `["PASS"]` or empty (= unfiltered).
    pub fn is_pass(&self) -> bool {
        self.filter.is_empty() || self.filter == ["PASS".to_string()]
    }
}
