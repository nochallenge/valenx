//! `[bio.mrbayes]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mrbayes.mcmc"
//!
//! [bio.mrbayes]
//! nexus      = "data.nex"
//! batch      = false                  # optional, defaults to false
//! extra_args = ["--no-autoclose"]     # optional, defaults to []
//! ```
//!
//! MrBayes' primary input is a NEXUS file with at least a DATA block
//! and (typically) a MRBAYES block embedding the model / MCMC
//! parameters and `mcmc` command. The adapter doesn't generate the
//! NEXUS; the user authors it and references it from `nexus`.
//!
//! `batch = true` adds `-i` (interactive-off) so MrBayes runs the
//! embedded commands and exits cleanly rather than waiting on stdin
//! at the prompt — the right default for non-interactive automation.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct MrBayesInput {
    /// Path to the NEXUS data file (relative to the case directory,
    /// or absolute). Must include a MrBayes block driving the run.
    pub nexus: PathBuf,
    /// When true, append `-i` so MrBayes runs the embedded
    /// commands non-interactively and exits cleanly. Defaults to
    /// false (mirrors MrBayes's own default).
    pub batch: bool,
    /// Additional CLI arguments appended to the invocation.
    pub extra_args: Vec<String>,
}

impl MrBayesInput {
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
            .and_then(|v| v.get("mrbayes"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mrbayes] section",
                    case_toml.display()
                ))
            })?;

        let nexus = block
            .get("nexus")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.mrbayes].nexus required")))?;
        if nexus.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mrbayes].nexus must not be empty"
            )));
        }

        let batch = block
            .get("batch")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mrbayes].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.mrbayes].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            nexus: PathBuf::from(nexus),
            batch,
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
        let d = tempdir("mrbayes-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mrbayes.mcmc"

[bio.mrbayes]
nexus = "data.nex"
"#,
        )
        .unwrap();
        let input = MrBayesInput::from_case_dir(&d).unwrap();
        assert_eq!(input.nexus, PathBuf::from("data.nex"));
        // Defaults: interactive (no -i), no extras.
        assert!(!input.batch);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_batch_mode() {
        // Non-interactive automation run. `batch = true` lifts the
        // `-i` flag in the adapter so MrBayes doesn't wait on stdin
        // at the prompt after running the embedded commands.
        let d = tempdir("mrbayes-batch");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mrbayes.mcmc"

[bio.mrbayes]
nexus      = "primates.nex"
batch      = true
extra_args = ["--no-autoclose"]
"#,
        )
        .unwrap();
        let input = MrBayesInput::from_case_dir(&d).unwrap();
        assert_eq!(input.nexus, PathBuf::from("primates.nex"));
        assert!(input.batch);
        assert_eq!(input.extra_args, vec!["--no-autoclose".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_nexus() {
        // The NEXUS file is the entire input — empty string means
        // MrBayes has nothing to run. Reject up front.
        let d = tempdir("mrbayes-nonex");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mrbayes.mcmc"

[bio.mrbayes]
nexus = ""
"#,
        )
        .unwrap();
        let err = MrBayesInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("nexus"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
