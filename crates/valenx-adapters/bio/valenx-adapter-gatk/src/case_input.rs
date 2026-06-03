//! `[bio.gatk]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "gatk.haplotype-caller"
//!
//! [bio.gatk]
//! reference   = "ref.fa"          # FASTA, must be indexed (.fai + .dict)
//! input_bam   = "aligned.bam"     # sorted + indexed BAM
//! output_vcf  = "calls.vcf"       # relative to workdir
//! intervals   = "regions.bed"     # optional BED for region restriction
//! java_heap   = "8g"              # required, e.g. "8g" / "16G" / "4096m"
//! extra_args  = ["--annotation", "Coverage"]   # optional, defaults to []
//! ```
//!
//! `java_heap` matches `^\d+[gmGM]$` — the bare-number-suffix form
//! GATK's `--java-options "-Xmx..."` accepts. We validate it here so a
//! typo doesn't surface as a Java startup failure mid-run.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct GatkInput {
    /// Path to the reference FASTA. GATK additionally requires the
    /// `.fai` index and `.dict` sequence dictionary to sit next to it
    /// — those are produced via `samtools faidx` and
    /// `gatk CreateSequenceDictionary` respectively, but the case
    /// input only carries the FASTA path itself.
    pub reference: PathBuf,
    /// Sorted + indexed BAM with read group(s) populated. Producing
    /// the BAM is upstream of the adapter (typically BWA → samtools
    /// sort → samtools index).
    pub input_bam: PathBuf,
    /// Output VCF path, relative to the workdir.
    pub output_vcf: PathBuf,
    /// Optional BED for restricting calling to specific regions.
    /// `None` runs whole-genome (or whole-BAM) calling.
    pub intervals: Option<PathBuf>,
    /// JVM heap size, e.g. "8g", "16G", "4096m". Format validated by
    /// [`is_valid_heap`].
    pub java_heap: String,
    pub extra_args: Vec<String>,
}

/// Returns true if `s` matches the GATK `--java-options "-Xmx..."`
/// heap-size shape: one or more digits followed by a single `g`, `m`,
/// `G`, or `M` suffix. Hand-rolled so we don't pay the regex-crate
/// dependency cost for one validator.
pub fn is_valid_heap(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }
    // Last byte must be the suffix; everything before must be digits.
    // Heap strings are always ASCII so byte indexing is fine.
    let bytes = s.as_bytes();
    let suffix = bytes[bytes.len() - 1];
    if !matches!(suffix, b'g' | b'G' | b'm' | b'M') {
        return false;
    }
    let prefix = &bytes[..bytes.len() - 1];
    if prefix.is_empty() {
        return false;
    }
    prefix.iter().all(|b| b.is_ascii_digit())
}

impl GatkInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display()))
        })?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("gatk"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.gatk] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.gatk].reference required (path to FASTA + .fai + .dict)"
                ))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.gatk].reference must not be empty"
            )));
        }

        let input_bam_str = block
            .get("input_bam")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.gatk].input_bam required (path to sorted+indexed BAM)"
                ))
            })?;
        if input_bam_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.gatk].input_bam must not be empty"
            )));
        }

        let output_vcf_str = block
            .get("output_vcf")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.gatk].output_vcf required (path to output VCF)"
                ))
            })?;
        if output_vcf_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.gatk].output_vcf must not be empty"
            )));
        }

        let intervals = block
            .get("intervals")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        let java_heap_str = block
            .get("java_heap")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.gatk].java_heap required (e.g. \"8g\" or \"4096m\")"
                ))
            })?;
        if !is_valid_heap(java_heap_str) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.gatk].java_heap `{java_heap_str}` invalid — \
                 expected digits followed by `g`, `G`, `m`, or `M` \
                 (e.g. \"8g\", \"4096m\")"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.gatk].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.gatk].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            reference: PathBuf::from(reference_str),
            input_bam: PathBuf::from(input_bam_str),
            output_vcf: PathBuf::from(output_vcf_str),
            intervals,
            java_heap: java_heap_str.to_string(),
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("gatk");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "gatk.haplotype-caller"

[bio.gatk]
reference  = "ref.fa"
input_bam  = "aligned.bam"
output_vcf = "calls.vcf"
java_heap  = "8g"
"#,
        )
        .unwrap();
        let input = GatkInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.input_bam, PathBuf::from("aligned.bam"));
        assert_eq!(input.output_vcf, PathBuf::from("calls.vcf"));
        assert_eq!(input.java_heap, "8g");
        assert_eq!(input.intervals, None);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_intervals() {
        // `intervals` and `extra_args` round-trip through, with a
        // capital-M heap to confirm the validator accepts the
        // case-insensitive suffix.
        let d = tempdir("gatk");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "gatk.haplotype-caller"

[bio.gatk]
reference   = "ref.fa"
input_bam   = "aligned.bam"
output_vcf  = "calls.vcf"
intervals   = "regions.bed"
java_heap   = "4096M"
extra_args  = ["--annotation", "Coverage"]
"#,
        )
        .unwrap();
        let input = GatkInput::from_case_dir(&d).unwrap();
        assert_eq!(input.intervals, Some(PathBuf::from("regions.bed")));
        assert_eq!(input.java_heap, "4096M");
        assert_eq!(
            input.extra_args,
            vec!["--annotation".to_string(), "Coverage".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_bad_heap() {
        // `8gb` is a common typo (the JVM accepts only `g`/`G`/`m`/`M`,
        // not the multi-letter `gb`/`mb` units). Catch it at parse
        // time so the user sees a clear error rather than a Java
        // startup failure.
        let d = tempdir("gatk");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "gatk.haplotype-caller"

[bio.gatk]
reference  = "ref.fa"
input_bam  = "aligned.bam"
output_vcf = "calls.vcf"
java_heap  = "8gb"
"#,
        )
        .unwrap();
        let err = GatkInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("java_heap"), "msg: {msg}");
        assert!(msg.contains("8gb"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("gatk");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = GatkInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.gatk]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn is_valid_heap_helper_unit_tests() {
        // Belt-and-braces direct unit coverage of the validator —
        // makes the rules explicit and protects against future
        // refactors that might widen the accepted shape.
        assert!(is_valid_heap("8g"));
        assert!(is_valid_heap("4096m"));
        assert!(is_valid_heap("16G"));
        assert!(is_valid_heap("1024M"));
        assert!(!is_valid_heap(""));
        assert!(!is_valid_heap("g"));
        assert!(!is_valid_heap("8gb"));
        assert!(!is_valid_heap("8 g"));
        assert!(!is_valid_heap("8.5g"));
        assert!(!is_valid_heap("8k"));
    }
}
