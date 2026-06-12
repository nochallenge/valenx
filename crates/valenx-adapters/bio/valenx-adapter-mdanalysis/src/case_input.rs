//! `[bio.mdanalysis]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mdanalysis.script"
//!
//! [bio.mdanalysis]
//! script  = "analyse_traj.py"
//! python  = "python3"          # optional, defaults to python3
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct MdAnalysisInput {
    pub script: PathBuf,
    pub python: String,
}

impl MdAnalysisInput {
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
            .and_then(|v| v.get("mdanalysis"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mdanalysis] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.mdanalysis].script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        Ok(Self {
            script: PathBuf::from(script),
            python,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case() {
        let d = tempdir("mdanalysis");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mdanalysis.script"

[bio.mdanalysis]
script = "rmsd.py"
"#,
        )
        .unwrap();
        let input = MdAnalysisInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("rmsd.py"));
        assert_eq!(input.python, "python3");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("mdanalysis");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = MdAnalysisInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.mdanalysis]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_python_override() {
        // MDAnalysis is sometimes pinned to a conda env distinct
        // from the system Python — verify a custom interpreter path
        // round-trips cleanly through the case-input parser.
        let d = tempdir("mdanalysis");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mdanalysis.script"

[bio.mdanalysis]
script = "analyse.py"
python = "/opt/conda/envs/mda/bin/python"
"#,
        )
        .unwrap();
        let input = MdAnalysisInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/mda/bin/python");
        let _ = std::fs::remove_dir_all(&d);
    }
}
