//! `[bio.nextflow]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "nextflow.run"
//!
//! [bio.nextflow]
//! pipeline   = "nf-core/rnaseq"             # required; registry slug or local .nf path
//! profile    = "docker"                      # optional; -profile <name>
//! resume     = true                          # optional; defaults to false
//! config     = "nextflow.config"             # optional; -c <file>
//! params     = { input = "samplesheet.csv" } # optional; flattens to `--key value`
//! extra_args = ["-with-report", "report.html"] # optional; pass-through
//! ```
//!
//! Nextflow accepts a registry slug (`nf-core/rnaseq`) or a path to a
//! local `.nf` file. The adapter resolves relative paths against the
//! case directory; bare slugs go through verbatim so the underlying
//! Nextflow CLI can fetch them from its registry.
//!
//! `params` lands as `--<key> <value>` pairs in the final invocation;
//! we keep them in a `BTreeMap` so iteration order is deterministic
//! (alphabetical by key) — important for reproducibility and for
//! diffing two `case.toml` files that drive the same pipeline.

use std::collections::BTreeMap;
use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct NextflowInput {
    pub pipeline: String,
    pub params: BTreeMap<String, String>,
    pub profile: Option<String>,
    pub resume: bool,
    pub config: Option<PathBuf>,
    pub extra_args: Vec<String>,
}

impl NextflowInput {
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
            .and_then(|v| v.get("nextflow"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.nextflow] section",
                    case_toml.display()
                ))
            })?;

        let pipeline_raw = block
            .get("pipeline")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.nextflow].pipeline required"))
            })?;
        let pipeline = pipeline_raw.trim().to_string();
        if pipeline.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.nextflow].pipeline must not be empty"
            )));
        }

        // params: optional inline TOML table → BTreeMap<String, String>.
        // We coerce values to strings so the downstream CLI flags are
        // a uniform shape; integers / floats / bools stringify the
        // obvious way (Nextflow accepts string-form values for every
        // primitive parameter).
        let params = match block.get("params") {
            Some(v) => {
                let table = v.as_table().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.nextflow].params must be a table of key=value entries"
                    ))
                })?;
                let mut out: BTreeMap<String, String> = BTreeMap::new();
                for (k, v) in table {
                    let s = match v {
                        toml::Value::String(s) => s.clone(),
                        toml::Value::Integer(i) => i.to_string(),
                        toml::Value::Float(f) => f.to_string(),
                        toml::Value::Boolean(b) => b.to_string(),
                        _ => {
                            return Err(AdapterError::Other(anyhow::anyhow!(
                                "[bio.nextflow].params.{k} must be a scalar (string, integer, float, or bool)"
                            )));
                        }
                    };
                    out.insert(k.clone(), s);
                }
                out
            }
            None => BTreeMap::new(),
        };

        let profile = match block.get("profile") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.nextflow].profile must be a string"))
                })?;
                if s.is_empty() {
                    None
                } else {
                    Some(s.to_string())
                }
            }
            None => None,
        };

        let resume = block
            .get("resume")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let config = match block.get("config") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.nextflow].config must be a string path"
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
                        "[bio.nextflow].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.nextflow].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            pipeline,
            params,
            profile,
            resume,
            config,
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
        // Just a pipeline reference — every other field defaults.
        let d = tempdir("nextflow");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nextflow.run"

[bio.nextflow]
pipeline = "nf-core/rnaseq"
"#,
        )
        .unwrap();
        let input = NextflowInput::from_case_dir(&d).unwrap();
        assert_eq!(input.pipeline, "nf-core/rnaseq");
        assert!(input.params.is_empty());
        assert!(input.profile.is_none());
        assert!(!input.resume);
        assert!(input.config.is_none());
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_params_and_profile() {
        // Typical nf-core invocation: pipeline + a docker profile +
        // `--input` / `--outdir` style param map. The BTreeMap
        // ensures the params come back in alphabetical-by-key order
        // regardless of TOML source order.
        let d = tempdir("nextflow");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nextflow.run"

[bio.nextflow]
pipeline = "nf-core/rnaseq"
profile  = "docker"
config   = "nextflow.config"

[bio.nextflow.params]
outdir = "results"
input  = "samplesheet.csv"
"#,
        )
        .unwrap();
        let input = NextflowInput::from_case_dir(&d).unwrap();
        assert_eq!(input.profile.as_deref(), Some("docker"));
        assert_eq!(input.config, Some(PathBuf::from("nextflow.config")));
        // BTreeMap iteration is alphabetical; `input` < `outdir`.
        let keys: Vec<&String> = input.params.keys().collect();
        assert_eq!(keys, vec!["input", "outdir"]);
        assert_eq!(
            input.params.get("input").map(String::as_str),
            Some("samplesheet.csv")
        );
        assert_eq!(
            input.params.get("outdir").map(String::as_str),
            Some("results")
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_resume() {
        let d = tempdir("nextflow");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nextflow.run"

[bio.nextflow]
pipeline   = "nf-core/sarek"
resume     = true
extra_args = ["-with-report", "report.html"]
"#,
        )
        .unwrap();
        let input = NextflowInput::from_case_dir(&d).unwrap();
        assert!(input.resume);
        assert_eq!(
            input.extra_args,
            vec!["-with-report".to_string(), "report.html".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_pipeline() {
        // Whitespace-only pipeline trims to empty and must error —
        // running `nextflow run ""` would otherwise produce a
        // confusing CLI error far downstream.
        let d = tempdir("nextflow");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nextflow.run"

[bio.nextflow]
pipeline = "   "
"#,
        )
        .unwrap();
        let err = NextflowInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("pipeline must not be empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
