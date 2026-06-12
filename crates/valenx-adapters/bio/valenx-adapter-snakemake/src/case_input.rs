//! `[bio.snakemake]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "snakemake.run"
//!
//! [bio.snakemake]
//! snakefile   = "Snakefile"            # required; path to the rulefile
//! targets     = ["all"]                # optional; rule / file targets
//! cores       = 4                      # optional, defaults to 1
//! use_conda   = true                   # optional; --use-conda
//! dry_run     = false                  # optional; -n
//! config_file = "config.yaml"          # optional; --configfile <path>
//! extra_args  = ["--rerun-incomplete"] # optional; pass-through
//! ```
//!
//! Snakemake is the rule-based / data-flow alternative to Nextflow.
//! The adapter just composes a `snakemake -s <Snakefile>` invocation
//! and lets the Snakemake CLI itself report back through stderr.
//!
//! Targets are arbitrary rule names *or* output file paths — Snakemake
//! resolves both — so we keep them as untyped strings rather than
//! forcing them through `PathBuf`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SnakemakeInput {
    pub snakefile: PathBuf,
    pub targets: Vec<String>,
    pub cores: u32,
    pub use_conda: bool,
    pub dry_run: bool,
    pub config_file: Option<PathBuf>,
    pub extra_args: Vec<String>,
}

impl SnakemakeInput {
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
            .and_then(|v| v.get("snakemake"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.snakemake] section",
                    case_toml.display()
                ))
            })?;

        let snakefile_str = block
            .get("snakefile")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.snakemake].snakefile required"))
            })?;
        if snakefile_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.snakemake].snakefile must not be empty"
            )));
        }

        let targets = match block.get("targets") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.snakemake].targets must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.snakemake].targets entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        let cores = match block.get("cores") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.snakemake].cores must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.snakemake].cores must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let use_conda = block
            .get("use_conda")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let dry_run = block
            .get("dry_run")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let config_file = match block.get("config_file") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.snakemake].config_file must be a string path"
                    ))
                })?;
                if s.is_empty() {
                    None
                } else {
                    Some(PathBuf::from(s))
                }
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.snakemake].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.snakemake].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            snakefile: PathBuf::from(snakefile_str),
            targets,
            cores,
            use_conda,
            dry_run,
            config_file,
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
        // Just a Snakefile — every other field defaults.
        let d = tempdir("snakemake");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "snakemake.run"

[bio.snakemake]
snakefile = "Snakefile"
"#,
        )
        .unwrap();
        let input = SnakemakeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.snakefile, PathBuf::from("Snakefile"));
        assert!(input.targets.is_empty());
        assert_eq!(input.cores, 1);
        assert!(!input.use_conda);
        assert!(!input.dry_run);
        assert!(input.config_file.is_none());
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_targets_and_conda() {
        // Conda-isolated rule run with explicit targets and an
        // outboard config file. The 8-core hint is the typical
        // local-laptop default for a multi-rule pipeline.
        let d = tempdir("snakemake");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "snakemake.run"

[bio.snakemake]
snakefile   = "rules/Snakefile"
targets     = ["align", "results/all.bam"]
cores       = 8
use_conda   = true
config_file = "config.yaml"
"#,
        )
        .unwrap();
        let input = SnakemakeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.snakefile, PathBuf::from("rules/Snakefile"));
        assert_eq!(
            input.targets,
            vec!["align".to_string(), "results/all.bam".to_string()]
        );
        assert_eq!(input.cores, 8);
        assert!(input.use_conda);
        assert_eq!(input.config_file, Some(PathBuf::from("config.yaml")));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_dry_run() {
        // Dry-run mode plus an extra-args pass-through. `-n` is the
        // canonical Snakemake "show me the plan" flag.
        let d = tempdir("snakemake");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "snakemake.run"

[bio.snakemake]
snakefile  = "Snakefile"
dry_run    = true
extra_args = ["--rerun-incomplete"]
"#,
        )
        .unwrap();
        let input = SnakemakeInput::from_case_dir(&d).unwrap();
        assert!(input.dry_run);
        assert_eq!(input.extra_args, vec!["--rerun-incomplete".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_cores() {
        // `--cores 0` would error far downstream in Snakemake itself
        // ("at least 1 core required"); reject upfront for a clearer
        // error path.
        let d = tempdir("snakemake");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "snakemake.run"

[bio.snakemake]
snakefile = "Snakefile"
cores     = 0
"#,
        )
        .unwrap();
        let err = SnakemakeInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("cores must be >= 1"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
