//! `[bio.pyrosetta]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "pyrosetta.script"
//!
//! [bio.pyrosetta]
//! script          = "design.py"
//! python          = "python3"          # optional, defaults to python3
//! input_pdb       = "input.pdb"        # optional — passed via params.json if present
//! output_basename = "design"
//! ```
//!
//! PyRosetta exposes the entire Rosetta C++ core through Python
//! bindings, letting users drive Rosetta from Python scripts rather
//! than XML protocols. The adapter stages the user's `.py` driver
//! (and optional input PDB) into the workdir and writes a flat
//! `valenx_params.json` so the script can read the parsed knobs
//! without re-parsing case.toml.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PyRosettaInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute).
    pub script: PathBuf,
    /// Python interpreter to invoke. Defaults to `python3`; pass an
    /// absolute path or environment-pinned name to override.
    pub python: String,
    /// Optional input PDB the script will operate on. Surfaced in
    /// `valenx_params.json` so the script can read it back without
    /// re-parsing `case.toml`. None when the script generates its
    /// own structures de novo.
    pub input_pdb: Option<PathBuf>,
    /// Output filename stem. Surfaced in `valenx_params.json`;
    /// `collect()` walks the workdir for `<basename>*.pdb` decoys
    /// and `*.sc` scorefiles.
    pub output_basename: String,
}

impl PyRosettaInput {
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
            .and_then(|v| v.get("pyrosetta"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.pyrosetta] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.pyrosetta].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pyrosetta].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.pyrosetta].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pyrosetta].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_pdb,
            output_basename: output_basename.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("pyrosetta-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pyrosetta.script"

[bio.pyrosetta]
script          = "design.py"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = PyRosettaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_pdb, None);
        assert_eq!(input.output_basename, "design");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_pdb() {
        // Script that operates on an existing structure (e.g. a
        // FastDesign / FastRelax driver). Pinning a non-default
        // Python interpreter exercises the override path too —
        // common in HPC installs where PyRosetta lives in a
        // dedicated conda env.
        let d = tempdir("pyrosetta-inpdb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pyrosetta.script"

[bio.pyrosetta]
script          = "fastrelax.py"
python          = "/opt/conda/envs/pyrosetta/bin/python"
input_pdb       = "scaffold.pdb"
output_basename = "relaxed"
"#,
        )
        .unwrap();
        let input = PyRosettaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input_pdb, Some(PathBuf::from("scaffold.pdb")));
        assert_eq!(input.python, "/opt/conda/envs/pyrosetta/bin/python");
        assert_eq!(input.output_basename, "relaxed");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename drives the `collect()` filter — empty
        // string would surface every PDB in the workdir, including
        // the staged input. Reject up front.
        let d = tempdir("pyrosetta-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pyrosetta.script"

[bio.pyrosetta]
script          = "design.py"
output_basename = ""
"#,
        )
        .unwrap();
        let err = PyRosettaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
