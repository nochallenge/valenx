//! `[bio.openmm]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "openmm.script"
//!
//! [bio.openmm]
//! script      = "minimise.py"
//! python      = "python3"          # optional, defaults to python3
//! output_pdb  = "minimised.pdb"    # optional, default minimised.pdb
//! output_dcd  = "trajectory.dcd"   # optional, default trajectory.dcd
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct OpenMmInput {
    pub script: PathBuf,
    pub python: String,
    /// Filename (relative to the workdir) the user's script writes
    /// the post-minimisation / post-equilibration coordinates to.
    /// `collect()` parses this PDB into a typed `valenx_bio`
    /// `Structure` when it exists.
    pub output_pdb: String,
    /// Filename (relative to the workdir) the user's script writes
    /// the simulation trajectory to. `collect()` lists this as a
    /// `Native` artifact; full DCD parsing lands in Task 11.
    pub output_dcd: String,
}

impl OpenMmInput {
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
            .and_then(|v| v.get("openmm"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.openmm] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.openmm].script required")))?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let output_pdb = block
            .get("output_pdb")
            .and_then(|v| v.as_str())
            .unwrap_or("minimised.pdb")
            .to_string();
        let output_dcd = block
            .get("output_dcd")
            .and_then(|v| v.as_str())
            .unwrap_or("trajectory.dcd")
            .to_string();
        Ok(Self {
            script: PathBuf::from(script),
            python,
            output_pdb,
            output_dcd,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_defaults() {
        let d = tempdir("openmm");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openmm.script"

[bio.openmm]
script = "minimise.py"
"#,
        )
        .unwrap();
        let input = OpenMmInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("minimise.py"));
        assert_eq!(input.python, "python3");
        // Defaults match the docstring schema.
        assert_eq!(input.output_pdb, "minimised.pdb");
        assert_eq!(input.output_dcd, "trajectory.dcd");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("openmm");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = OpenMmInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.openmm]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_output_filename_overrides() {
        // A simulation that writes "equilibrated.pdb" + "run.dcd"
        // — verify both overrides round-trip.
        let d = tempdir("openmm");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openmm.script"

[bio.openmm]
script     = "run.py"
python     = "/opt/conda/envs/openmm/bin/python"
output_pdb = "equilibrated.pdb"
output_dcd = "run.dcd"
"#,
        )
        .unwrap();
        let input = OpenMmInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/openmm/bin/python");
        assert_eq!(input.output_pdb, "equilibrated.pdb");
        assert_eq!(input.output_dcd, "run.dcd");
        let _ = std::fs::remove_dir_all(&d);
    }
}
