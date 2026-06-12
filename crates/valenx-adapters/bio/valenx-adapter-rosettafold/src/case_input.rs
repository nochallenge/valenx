//! `[bio.rosettafold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "rosettafold.predict"
//!
//! [bio.rosettafold]
//! script           = "predict.py"
//! python           = "python3"        # optional, default python3
//! fasta            = "query.fasta"
//! output_basename  = "predicted"
//! ```
//!
//! RoseTTAFold has no canonical pip-installable distribution — the
//! Baker lab's repo is a research codebase that the user clones,
//! installs into a conda env, and drives via per-site predict scripts.
//! The schema mirrors that: the user supplies their own `predict.py`,
//! a FASTA query, and an `output_basename` prefix the script writes
//! its outputs under (`<basename>.pdb`, `<basename>.npz`, etc.).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RoseTTAFoldInput {
    pub script: PathBuf,
    pub python: String,
    pub fasta: PathBuf,
    pub output_basename: String,
}

impl RoseTTAFoldInput {
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
            .and_then(|v| v.get("rosettafold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.rosettafold] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rosettafold].script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let fasta = block.get("fasta").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!("[bio.rosettafold].fasta required"))
        })?;
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.rosettafold].output_basename required"
                ))
            })?
            .to_string();
        Ok(Self {
            script: PathBuf::from(script),
            python,
            fasta: PathBuf::from(fasta),
            output_basename,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_defaults() {
        let d = tempdir("rosettafold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosettafold.predict"

[bio.rosettafold]
script          = "predict.py"
fasta           = "query.fasta"
output_basename = "predicted"
"#,
        )
        .unwrap();
        let input = RoseTTAFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("predict.py"));
        assert_eq!(input.fasta, PathBuf::from("query.fasta"));
        assert_eq!(input.output_basename, "predicted");
        // Default python interpreter when omitted.
        assert_eq!(input.python, "python3");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("rosettafold-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = RoseTTAFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.rosettafold]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_python_override() {
        // RoseTTAFold is typically pinned to a conda env distinct from
        // the system Python — the user supplies the env's interpreter
        // path. Round-trip cleanly through the case-input parser.
        let d = tempdir("rosettafold-py");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosettafold.predict"

[bio.rosettafold]
script          = "predict.py"
python          = "/opt/conda/envs/RoseTTAFold/bin/python"
fasta           = "query.fasta"
output_basename = "predicted"
"#,
        )
        .unwrap();
        let input = RoseTTAFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/RoseTTAFold/bin/python");
        let _ = std::fs::remove_dir_all(&d);
    }
}
