//! `[bio.cromwell]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "cromwell.run"
//!
//! [bio.cromwell]
//! jar             = "/path/to/cromwell-86.jar"  # required; Broad's distribution JAR
//! workflow        = "workflow.wdl"              # required; WDL source
//! # inputs        = "inputs.json"               # optional; workflow inputs JSON
//! output_basename = "metadata"                  # required; collect() filter for metadata JSON
//! action          = "run"                       # optional; one of "run" | "submit" | "validate"
//! extra_args      = []                          # optional; extra Cromwell args
//! ```
//!
//! Cromwell ships as a Java JAR (`cromwell-<version>.jar`) — there's
//! no `cromwell` launcher binary on PATH. The user supplies the jar
//! path via `jar`; the adapter probes that `java` itself is installed
//! and composes `java -jar <jar> <action> <workflow> [-i <inputs>]
//! [extras...]` from `prepare()`.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical Cromwell action list. Module-public so the UI can surface
/// the supported subcommands without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["run", "submit", "validate"];

#[derive(Clone, Debug, PartialEq)]
pub struct CromwellInput {
    /// Absolute (or case-relative) path to the `cromwell-<version>.jar`
    /// distributed by the Broad Institute.
    pub jar: PathBuf,
    /// Path to the WDL workflow source.
    pub workflow: PathBuf,
    /// Optional path to a workflow inputs JSON. When present, surfaces
    /// as the `-i <inputs>` flag pair.
    pub inputs: Option<PathBuf>,
    /// Filename stem `collect()` uses to filter metadata JSON outputs.
    pub output_basename: String,
    /// One of: "run" | "submit" | "validate". Defaults to "run".
    pub action: String,
    /// Additional CLI arguments appended to the `java -jar` call —
    /// useful for `--metadata-output <path>` or `-Dconfig.file=<file>`.
    pub extra_args: Vec<String>,
}

impl CromwellInput {
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
            .and_then(|v| v.get("cromwell"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.cromwell] section",
                    case_toml.display()
                ))
            })?;

        let jar_str = block.get("jar").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.cromwell].jar required (path to cromwell-<version>.jar)"
            ))
        })?;
        if jar_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cromwell].jar must not be empty"
            )));
        }

        let workflow_str = block
            .get("workflow")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cromwell].workflow required"))
            })?;
        if workflow_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cromwell].workflow must not be empty"
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
                AdapterError::Other(anyhow::anyhow!("[bio.cromwell].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cromwell].output_basename must not be empty"
            )));
        }

        let action = match block.get("action") {
            Some(v) => v
                .as_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.cromwell].action must be a string"))
                })?
                .to_string(),
            None => "run".to_string(),
        };
        if !SUPPORTED_ACTIONS.contains(&action.as_str()) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cromwell].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.cromwell].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.cromwell].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            jar: PathBuf::from(jar_str),
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
    fn parses_minimal_with_default_action() {
        // Just jar + workflow + output_basename — action defaults to
        // "run", inputs is None, extras empty.
        let d = tempdir("cromwell-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cromwell.run"

[bio.cromwell]
jar             = "/opt/cromwell/cromwell-86.jar"
workflow        = "workflow.wdl"
output_basename = "metadata"
"#,
        )
        .unwrap();
        let input = CromwellInput::from_case_dir(&d).unwrap();
        assert_eq!(input.jar, PathBuf::from("/opt/cromwell/cromwell-86.jar"));
        assert_eq!(input.workflow, PathBuf::from("workflow.wdl"));
        assert_eq!(input.inputs, None);
        assert_eq!(input.output_basename, "metadata");
        assert_eq!(input.action, "run");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_inputs_and_action() {
        // `validate` action with an inputs JSON and pass-through extras.
        let d = tempdir("cromwell-full");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cromwell.validate"

[bio.cromwell]
jar             = "/opt/cromwell/cromwell-86.jar"
workflow        = "workflow.wdl"
inputs          = "inputs.json"
output_basename = "metadata"
action          = "validate"
extra_args      = ["-Dconfig.file=cromwell.conf"]
"#,
        )
        .unwrap();
        let input = CromwellInput::from_case_dir(&d).unwrap();
        assert_eq!(input.jar, PathBuf::from("/opt/cromwell/cromwell-86.jar"));
        assert_eq!(input.workflow, PathBuf::from("workflow.wdl"));
        assert_eq!(input.inputs, Some(PathBuf::from("inputs.json")));
        assert_eq!(input.output_basename, "metadata");
        assert_eq!(input.action, "validate");
        assert_eq!(
            input.extra_args,
            vec!["-Dconfig.file=cromwell.conf".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `server` is a real Cromwell subcommand but isn't one this
        // adapter wraps — it spins up an HTTP daemon, which doesn't
        // fit the prepare/run/collect cycle. Reject up front rather
        // than forwarding blindly.
        let d = tempdir("cromwell-badaction");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cromwell.server"

[bio.cromwell]
jar             = "/opt/cromwell/cromwell-86.jar"
workflow        = "workflow.wdl"
output_basename = "metadata"
action          = "server"
"#,
        )
        .unwrap();
        let err = CromwellInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("validate"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
