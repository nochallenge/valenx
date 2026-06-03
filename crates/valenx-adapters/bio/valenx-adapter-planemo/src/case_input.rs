//! `[bio.planemo]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "planemo.run"
//!
//! [bio.planemo]
//! workflow        = "workflow.ga"   # required; .ga / .gxwf.yml Galaxy workflow
//! # inputs        = "inputs.json"   # optional; tool / workflow inputs JSON
//! output_basename = "report"        # required; collect() filter for outputs
//! action          = "run"           # optional; one of "run" | "test" | "lint"
//! extra_args      = []              # optional; extra Planemo args
//! ```
//!
//! Planemo is the Galaxy project's CLI for developing tools and running
//! workflows without standing up a Galaxy server. The same binary
//! handles tool linting, test execution, and one-off workflow runs;
//! `action` selects the subcommand the adapter wraps.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical Planemo action list. Module-public so the UI can surface
/// the supported values without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["run", "test", "lint"];

#[derive(Clone, Debug, PartialEq)]
pub struct PlanemoInput {
    pub workflow: PathBuf,
    pub inputs: Option<PathBuf>,
    pub output_basename: String,
    /// One of: "run" | "test" | "lint". Defaults to "run".
    pub action: String,
    pub extra_args: Vec<String>,
}

impl PlanemoInput {
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
            .and_then(|v| v.get("planemo"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.planemo] section",
                    case_toml.display()
                ))
            })?;

        let workflow_str = block
            .get("workflow")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.planemo].workflow required"))
            })?;
        if workflow_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.planemo].workflow must not be empty"
            )));
        }

        let inputs = block
            .get("inputs")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.planemo].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.planemo].output_basename must not be empty"
            )));
        }

        let action = match block.get("action") {
            Some(v) => v
                .as_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.planemo].action must be a string"))
                })?
                .to_string(),
            None => "run".to_string(),
        };
        if !SUPPORTED_ACTIONS.contains(&action.as_str()) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.planemo].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.planemo].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.planemo].extra_args entries must be strings"
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
            output_basename: output_basename.to_string(),
            action,
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
        // Just workflow + output_basename — action defaults to "run",
        // inputs is None, extras empty.
        let d = tempdir("planemo");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "planemo.run"

[bio.planemo]
workflow        = "workflow.ga"
output_basename = "report"
"#,
        )
        .unwrap();
        let input = PlanemoInput::from_case_dir(&d).unwrap();
        assert_eq!(input.workflow, PathBuf::from("workflow.ga"));
        assert_eq!(input.inputs, None);
        assert_eq!(input.output_basename, "report");
        assert_eq!(input.action, "run");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_inputs_and_action() {
        // `test` action with an inputs JSON and pass-through extras.
        let d = tempdir("planemo");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "planemo.test"

[bio.planemo]
workflow        = "workflow.gxwf.yml"
inputs          = "inputs.json"
output_basename = "results"
action          = "test"
extra_args      = ["--no_cleanup"]
"#,
        )
        .unwrap();
        let input = PlanemoInput::from_case_dir(&d).unwrap();
        assert_eq!(input.workflow, PathBuf::from("workflow.gxwf.yml"));
        assert_eq!(input.inputs, Some(PathBuf::from("inputs.json")));
        assert_eq!(input.output_basename, "results");
        assert_eq!(input.action, "test");
        assert_eq!(input.extra_args, vec!["--no_cleanup".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `serve` is a real Planemo subcommand but isn't one this
        // adapter wraps — must be rejected up front rather than forwarded
        // to the binary blindly.
        let d = tempdir("planemo");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "planemo.serve"

[bio.planemo]
workflow        = "workflow.ga"
output_basename = "report"
action          = "serve"
"#,
        )
        .unwrap();
        let err = PlanemoInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("lint"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
