//! `[bio.deepvariant]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "deepvariant.call"
//!
//! [bio.deepvariant]
//! reference   = "ref.fa"
//! input_bam   = "aligned.bam"
//! output_vcf  = "calls.vcf"
//! model_type  = "WGS"           # one of WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA
//! num_shards  = 8               # ≥ 1
//! extra_args  = ["--vcf_stats_report=true"]   # optional, defaults to []
//! ```
//!
//! `model_type` selects the trained DeepVariant model. The five values
//! map directly to the Google-published model bundles:
//!
//! - `WGS` — whole-genome sequencing, Illumina short reads (the
//!   default).
//! - `WES` — whole-exome sequencing, Illumina short reads.
//! - `PACBIO` — PacBio HiFi long reads.
//! - `ONT_R104` — Oxford Nanopore Technologies R10.4 long reads.
//! - `HYBRID_PACBIO_ILLUMINA` — the experimental hybrid model that
//!   takes Illumina + PacBio reads together.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// The five DeepVariant model bundles the Broad and Google publish.
/// Module-public so the UI can surface the supported values without
/// redefining them here.
pub const SUPPORTED_MODEL_TYPES: &[&str] =
    &["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"];

#[derive(Clone, Debug, PartialEq)]
pub struct DeepVariantInput {
    pub reference: PathBuf,
    pub input_bam: PathBuf,
    pub output_vcf: PathBuf,
    /// One of `WGS`, `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA`.
    pub model_type: String,
    /// Number of parallel shards DeepVariant uses for the
    /// make_examples step. ≥ 1; the user typically picks
    /// `nproc - 1`.
    pub num_shards: u32,
    pub extra_args: Vec<String>,
}

impl DeepVariantInput {
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
            .and_then(|v| v.get("deepvariant"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.deepvariant] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.deepvariant].reference required (path to FASTA)"
                ))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.deepvariant].reference must not be empty"
            )));
        }

        let input_bam_str = block
            .get("input_bam")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.deepvariant].input_bam required (path to sorted+indexed BAM)"
                ))
            })?;
        if input_bam_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.deepvariant].input_bam must not be empty"
            )));
        }

        let output_vcf_str = block
            .get("output_vcf")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.deepvariant].output_vcf required (path to output VCF)"
                ))
            })?;
        if output_vcf_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.deepvariant].output_vcf must not be empty"
            )));
        }

        let model_type = block
            .get("model_type")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.deepvariant].model_type required (one of {SUPPORTED_MODEL_TYPES:?})"
                ))
            })?;
        if !SUPPORTED_MODEL_TYPES.contains(&model_type) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.deepvariant].model_type `{model_type}` not recognised — \
                 expected one of {SUPPORTED_MODEL_TYPES:?}"
            )));
        }

        let num_shards = match block.get("num_shards") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.deepvariant].num_shards must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.deepvariant].num_shards must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.deepvariant].num_shards required (>= 1)"
                )));
            }
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.deepvariant].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.deepvariant].extra_args entries must be strings"
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
            model_type: model_type.to_string(),
            num_shards,
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
        let d = tempdir("deepvariant");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "deepvariant.call"

[bio.deepvariant]
reference   = "ref.fa"
input_bam   = "aligned.bam"
output_vcf  = "calls.vcf"
model_type  = "WGS"
num_shards  = 8
"#,
        )
        .unwrap();
        let input = DeepVariantInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.input_bam, PathBuf::from("aligned.bam"));
        assert_eq!(input.output_vcf, PathBuf::from("calls.vcf"));
        assert_eq!(input.model_type, "WGS");
        assert_eq!(input.num_shards, 8);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_model_type() {
        // `WGS_v2` isn't a published bundle — must be rejected up
        // front rather than failing inside DeepVariant's Python entry
        // point.
        let d = tempdir("deepvariant");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "deepvariant.call"

[bio.deepvariant]
reference   = "ref.fa"
input_bam   = "aligned.bam"
output_vcf  = "calls.vcf"
model_type  = "WGS_v2"
num_shards  = 8
"#,
        )
        .unwrap();
        let err = DeepVariantInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("PACBIO"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_shards() {
        // num_shards = 0 would deadlock the make_examples step (no
        // workers spawned); reject up front.
        let d = tempdir("deepvariant");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "deepvariant.call"

[bio.deepvariant]
reference   = "ref.fa"
input_bam   = "aligned.bam"
output_vcf  = "calls.vcf"
model_type  = "WGS"
num_shards  = 0
"#,
        )
        .unwrap();
        let err = DeepVariantInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("num_shards"), "msg: {msg}");
        assert!(msg.contains(">= 1"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
