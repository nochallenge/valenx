//! `[bio.cwltool]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "cwltool.run"
//!
//! [bio.cwltool]
//! workflow   = "workflow.cwl"  # required; .cwl tool / workflow document
//! # inputs   = "inputs.json"   # optional; CWL input-object JSON / YAML
//! output_dir = "results"       # required; --outdir target subdir
//! extra_args = []              # optional; extra cwltool args
//! ```
//!
//! cwltool is the reference implementation of the
//! [Common Workflow Language](https://www.commonwl.org/) — the
//! cross-tool standard for describing analytical workflows in YAML /
//! JSON. The adapter composes a
//! `cwltool --outdir <output_dir> [extras...] <workflow> [inputs]`
//! invocation; cwltool itself stages tools (in-process or via Docker /
//! Singularity / podman per the workflow's `DockerRequirement`) and
//! writes the workflow's declared `outputs` into `<output_dir>/`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CwltoolInput {
    pub workflow: PathBuf,
    pub inputs: Option<PathBuf>,
    pub output_dir: String,
    pub extra_args: Vec<String>,
}

impl CwltoolInput {
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
            .and_then(|v| v.get("cwltool"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.cwltool] section",
                    case_toml.display()
                ))
            })?;

        let workflow_str = block
            .get("workflow")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cwltool].workflow required"))
            })?;
        if workflow_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cwltool].workflow must not be empty"
            )));
        }

        // `inputs` is the optional CWL input-object document
        // (`.json` / `.yaml`). When omitted, cwltool will still run a
        // workflow that takes no inputs (or only defaulted inputs).
        let inputs = block
            .get("inputs")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        let output_dir = block
            .get("output_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cwltool].output_dir required"))
            })?;
        if output_dir.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cwltool].output_dir must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.cwltool].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.cwltool].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            workflow: PathBuf::from(workflow_str),
            inputs,
            output_dir: output_dir.to_string(),
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
        // Just workflow + output_dir — inputs is None, extras empty.
        // Mirrors the canonical "run a self-contained CWL workflow"
        // setup where the workflow has no required external inputs.
        let d = tempdir("cwltool");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cwltool.run"

[bio.cwltool]
workflow   = "workflow.cwl"
output_dir = "results"
"#,
        )
        .unwrap();
        let input = CwltoolInput::from_case_dir(&d).unwrap();
        assert_eq!(input.workflow, PathBuf::from("workflow.cwl"));
        assert_eq!(input.inputs, None);
        assert_eq!(input.output_dir, "results");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_inputs_and_extras() {
        // Workflow + inputs JSON + pass-through cwltool flags.
        // `--parallel` and `--cachedir` are common for batch reruns
        // of a multi-step CWL pipeline.
        let d = tempdir("cwltool");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cwltool.run"

[bio.cwltool]
workflow   = "pipelines/main.cwl"
inputs     = "inputs.json"
output_dir = "out"
extra_args = ["--parallel", "--cachedir", ".cwl-cache"]
"#,
        )
        .unwrap();
        let input = CwltoolInput::from_case_dir(&d).unwrap();
        assert_eq!(input.workflow, PathBuf::from("pipelines/main.cwl"));
        assert_eq!(input.inputs, Some(PathBuf::from("inputs.json")));
        assert_eq!(input.output_dir, "out");
        assert_eq!(
            input.extra_args,
            vec![
                "--parallel".to_string(),
                "--cachedir".to_string(),
                ".cwl-cache".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_workflow() {
        // `[bio.cwltool]` exists but `workflow` is missing — surface
        // the missing-field error up front rather than letting
        // prepare() emit a confusing path error.
        let d = tempdir("cwltool");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cwltool.run"

[bio.cwltool]
output_dir = "results"
"#,
        )
        .unwrap();
        let err = CwltoolInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("workflow required"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
