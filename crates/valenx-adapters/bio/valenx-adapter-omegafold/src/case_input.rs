//! `[bio.omegafold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "omegafold.predict"
//!
//! [bio.omegafold]
//! fasta            = "query.fasta"
//! output_basename  = "predicted"
//! python           = "python3"        # optional, default python3
//! # model_dir      = "/path/to/omegafold/checkpoints"
//! ```
//!
//! Unlike ESMFold / RoseTTAFold, OmegaFold ships as an installed Python
//! package with its own CLI (`omegafold <fasta> <output_dir>`), so the
//! user doesn't need to provide a predict script. The `python` knob is
//! kept around for the fallback path (`python -m omegafold ...`) when
//! the standalone `omegafold` binary isn't on PATH but the package was
//! installed into a Python environment. `model_dir` lets the user point
//! at a pre-downloaded checkpoint directory so OmegaFold doesn't try to
//! re-fetch the weights every run.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct OmegaFoldInput {
    pub fasta: PathBuf,
    pub output_basename: String,
    pub python: String,
    pub model_dir: Option<PathBuf>,
}

impl OmegaFoldInput {
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
            .and_then(|v| v.get("omegafold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.omegafold] section",
                    case_toml.display()
                ))
            })?;
        let fasta = block.get("fasta").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!("[bio.omegafold].fasta required"))
        })?;
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.omegafold].output_basename required"))
            })?
            .to_string();
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let model_dir = block
            .get("model_dir")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        Ok(Self {
            fasta: PathBuf::from(fasta),
            output_basename,
            python,
            model_dir,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_defaults() {
        let d = tempdir("omegafold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "omegafold.predict"

[bio.omegafold]
fasta           = "query.fasta"
output_basename = "predicted"
"#,
        )
        .unwrap();
        let input = OmegaFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.fasta, PathBuf::from("query.fasta"));
        assert_eq!(input.output_basename, "predicted");
        // Default python interpreter, no model_dir override.
        assert_eq!(input.python, "python3");
        assert_eq!(input.model_dir, None);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("omegafold-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = OmegaFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.omegafold]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_model_dir_and_python_override() {
        // Pre-downloaded weights live outside of the case directory, so
        // `model_dir` is typically absolute. The python interpreter
        // gets overridden when OmegaFold lives in its own conda env.
        let d = tempdir("omegafold-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "omegafold.predict"

[bio.omegafold]
fasta           = "query.fasta"
output_basename = "predicted"
python          = "/opt/conda/envs/omegafold/bin/python"
model_dir       = "/data/omegafold/checkpoints"
"#,
        )
        .unwrap();
        let input = OmegaFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/omegafold/bin/python");
        assert_eq!(
            input.model_dir,
            Some(PathBuf::from("/data/omegafold/checkpoints"))
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
