//! `[bio.bcftools]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "bcftools.call"
//!
//! [bio.bcftools]
//! action     = "call"             # one of "view", "call", "filter", "concat"
//! input      = "aligned.bam"      # required for view/filter (VCF) and call (BAM)
//! inputs     = ["a.vcf", "b.vcf"] # required for concat (>= 2 entries)
//! output     = "out.vcf"          # required for all actions
//! reference  = "ref.fa"           # required only for `call`
//! threads    = 4                  # optional, defaults to 1
//! extra_args = ["-v"]             # optional, defaults to []
//! ```
//!
//! `action` selects which bcftools subcommand the adapter wraps:
//!
//! - `view`   — VCF/BCF↔VCF/BCF conversion or subsetting. `input`
//!   required.
//! - `call`   — variant calling from a BAM. `input` (BAM) and
//!   `reference` (FASTA) both required.
//! - `filter` — apply soft / hard filters to a VCF. `input` required.
//! - `concat` — concatenate multiple VCFs/BCFs from the same sample.
//!   `inputs` must hold ≥ 2 entries; `input` is forbidden.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical bcftools action list. Module-public so the UI can surface
/// the supported values without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["view", "call", "filter", "concat"];

#[derive(Clone, Debug, PartialEq)]
pub struct BcftoolsInput {
    /// One of: "view" | "call" | "filter" | "concat".
    pub action: String,
    /// For view/filter: input VCF/BCF path.
    /// For call: input BAM path.
    /// For concat: leave empty (reads `inputs`).
    pub input: Option<PathBuf>,
    /// For concat: inputs to concatenate.
    pub inputs: Vec<PathBuf>,
    /// Output filename (relative to workdir). Required for all actions.
    pub output: PathBuf,
    /// For call: reference FASTA. Otherwise None.
    pub reference: Option<PathBuf>,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl BcftoolsInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display())))?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("bcftools"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.bcftools] section",
                    case_toml.display()
                ))
            })?;

        let action = block
            .get("action")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.bcftools].action required (one of {SUPPORTED_ACTIONS:?})"
                ))
            })?;
        if !SUPPORTED_ACTIONS.contains(&action) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bcftools].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        // `output` is required for every action.
        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.bcftools].output required (path to output file)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bcftools].output must not be empty"
            )));
        }
        let output = PathBuf::from(output_str);

        // `input` (single-path scalar). Optional in the TOML but
        // validated per-action below.
        let input = block
            .get("input")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        // `inputs` (array of paths) — only meaningful for `concat`.
        let inputs: Vec<PathBuf> = match block.get("inputs") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].inputs must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.bcftools].inputs entries must be strings"
                        ))
                    })?;
                    if s.is_empty() {
                        return Err(AdapterError::Other(anyhow::anyhow!(
                            "[bio.bcftools].inputs entries must not be empty"
                        )));
                    }
                    out.push(PathBuf::from(s));
                }
                out
            }
            None => Vec::new(),
        };

        // `reference` is required only for `call`.
        let reference = block
            .get("reference")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        // Per-action validation.
        match action {
            "view" | "filter" => {
                if input.is_none() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].input required for action `{action}`"
                    )));
                }
                if !inputs.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].inputs is only valid for action `concat`; \
                         action `{action}` uses `input`"
                    )));
                }
            }
            "call" => {
                if input.is_none() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].input required for action `call` \
                         (path to BAM)"
                    )));
                }
                if reference.is_none() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].reference required for action `call` \
                         (path to FASTA)"
                    )));
                }
                if !inputs.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].inputs is only valid for action `concat`; \
                         action `call` uses `input`"
                    )));
                }
            }
            "concat" => {
                if input.is_some() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].input is not valid for action `concat`; \
                         use `inputs` (array of >= 2 paths) instead"
                    )));
                }
                if inputs.len() < 2 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].inputs must contain >= 2 entries for \
                         action `concat`, got {}",
                        inputs.len()
                    )));
                }
            }
            _ => unreachable!("action validated against SUPPORTED_ACTIONS above"),
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].threads must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.bcftools].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.bcftools].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            action: action.to_string(),
            input,
            inputs,
            output,
            reference,
            threads,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_view_minimal() {
        // view with input + output. Defaults: 1 thread, no extras, no
        // reference, no inputs[].
        let d = tempdir("bcftools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bcftools.view"

[bio.bcftools]
action = "view"
input  = "in.vcf"
output = "out.vcf"
"#,
        )
        .unwrap();
        let input = BcftoolsInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "view");
        assert_eq!(input.input, Some(PathBuf::from("in.vcf")));
        assert_eq!(input.output, PathBuf::from("out.vcf"));
        assert!(input.inputs.is_empty());
        assert_eq!(input.reference, None);
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_call_with_reference() {
        // call requires both `input` (BAM) and `reference` (FASTA).
        let d = tempdir("bcftools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bcftools.call"

[bio.bcftools]
action     = "call"
input      = "aligned.bam"
output     = "calls.vcf"
reference  = "ref.fa"
threads    = 8
extra_args = ["-v"]
"#,
        )
        .unwrap();
        let input = BcftoolsInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "call");
        assert_eq!(input.input, Some(PathBuf::from("aligned.bam")));
        assert_eq!(input.reference, Some(PathBuf::from("ref.fa")));
        assert_eq!(input.output, PathBuf::from("calls.vcf"));
        assert_eq!(input.threads, 8);
        assert_eq!(input.extra_args, vec!["-v".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_concat_with_inputs() {
        // concat takes `inputs = [...]` (>= 2) and rejects scalar
        // `input`. Output is required.
        let d = tempdir("bcftools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bcftools.concat"

[bio.bcftools]
action = "concat"
inputs = ["chunk1.vcf", "chunk2.vcf", "chunk3.vcf"]
output = "merged.vcf"
"#,
        )
        .unwrap();
        let input = BcftoolsInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "concat");
        assert_eq!(input.input, None);
        assert_eq!(input.inputs.len(), 3);
        assert_eq!(input.inputs[0], PathBuf::from("chunk1.vcf"));
        assert_eq!(input.inputs[2], PathBuf::from("chunk3.vcf"));
        assert_eq!(input.output, PathBuf::from("merged.vcf"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_call_without_reference() {
        // `call` is the only action that needs a reference; the
        // adapter must reject the case at parse time rather than
        // letting bcftools crash partway through with a half-written
        // output.
        let d = tempdir("bcftools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bcftools.call"

[bio.bcftools]
action = "call"
input  = "aligned.bam"
output = "calls.vcf"
"#,
        )
        .unwrap();
        let err = BcftoolsInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("reference"), "msg: {msg}");
        assert!(msg.contains("call"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `bcftools mpileup` is a real subcommand but isn't one this
        // adapter wraps — must be rejected up front.
        let d = tempdir("bcftools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bcftools.mpileup"

[bio.bcftools]
action = "mpileup"
input  = "aligned.bam"
output = "out.vcf"
"#,
        )
        .unwrap();
        let err = BcftoolsInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("concat"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
