//! `[bio.mdtraj]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mdtraj.analyze"
//!
//! [bio.mdtraj]
//! script           = "analyze.py"
//! python           = "python3"           # optional, defaults to "python3"
//! trajectory       = "md.xtc"
//! topology         = "md.pdb"
//! output_basename  = "analysis"
//! ```
//!
//! MDTraj is the Pande / VanderSpoel / Beauchamp lab's Python library
//! for molecular-dynamics trajectory analysis — sister to MDAnalysis,
//! with wider format support (XTC / TRR / DCD / NetCDF / HDF5 / LAMMPS
//! / TNG / ...) and tighter OpenMM integration. The user authors an
//! `analyze.py` driver that does `import mdtraj`, loads `trajectory`
//! against `topology`, runs the analysis (RMSD / RMSF / SASA / radius
//! of gyration / hydrogen bonds / dihedrals / contact maps / ...), and
//! writes results into the workdir. We just stage the three inputs and
//! invoke `python <script>` after dropping a `valenx_params.json` next
//! to it so the script can read the parsed knobs without re-parsing
//! case.toml.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct MdtrajInput {
    /// Path to the user-authored Python driver script (relative to the
    /// case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3` so
    /// the adapter walks PATH; can be pinned to an absolute path for
    /// users with multiple Python installs / venvs.
    pub python: String,
    /// Path to the input trajectory file (relative to the case
    /// directory, or absolute). MDTraj sniffs the format from the
    /// extension — `.xtc`, `.trr`, `.dcd`, `.nc`, `.h5`, `.lammpstrj`,
    /// `.tng`, `.pdb` (multi-model), and others all work.
    pub trajectory: PathBuf,
    /// Path to the topology file (relative to the case directory, or
    /// absolute). Required because most trajectory formats (XTC / TRR /
    /// DCD / NetCDF) carry only coordinates — bonds, residues, and
    /// chains live in the separate topology (`.pdb` / `.psf` /
    /// `.prmtop` / `.gro` / ...). MDTraj uses it to populate the
    /// `Trajectory.topology` attribute.
    pub topology: PathBuf,
    /// Filename stem for outputs. The script writes `<basename>*.csv`
    /// (analysis tables), `<basename>*.npz` (NumPy archives),
    /// `<basename>*.h5` (MDTraj HDF5), and `<basename>*.png` (plots)
    /// into the workdir.
    pub output_basename: String,
}

impl MdtrajInput {
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
            .and_then(|v| v.get("mdtraj"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mdtraj] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.mdtraj].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mdtraj].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python
        // tolerates other extensions but `import mdtraj` workflows
        // are conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mdtraj].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let trajectory = block
            .get("trajectory")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.mdtraj].trajectory required (path to trajectory file)"
                ))
            })?;
        if trajectory.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mdtraj].trajectory must not be empty"
            )));
        }

        let topology = block
            .get("topology")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.mdtraj].topology required (path to topology file)"
                ))
            })?;
        if topology.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mdtraj].topology must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.mdtraj].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mdtraj].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            trajectory: PathBuf::from(trajectory),
            topology: PathBuf::from(topology),
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
        let d = tempdir("mdtraj-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mdtraj.analyze"

[bio.mdtraj]
script          = "analyze.py"
trajectory      = "md.xtc"
topology        = "md.pdb"
output_basename = "analysis"
"#,
        )
        .unwrap();
        let input = MdtrajInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("analyze.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.trajectory, PathBuf::from("md.xtc"));
        assert_eq!(input.topology, PathBuf::from("md.pdb"));
        assert_eq!(input.output_basename, "analysis");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_python_override() {
        // MDTraj is sometimes pinned to a conda env distinct from the
        // system Python (the conda-forge build is the easiest install
        // on macOS / Windows) — verify a custom interpreter path
        // round-trips cleanly through the case-input parser.
        let d = tempdir("mdtraj-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mdtraj.analyze"

[bio.mdtraj]
script          = "rmsd.py"
python          = "/opt/conda/envs/mdtraj/bin/python"
trajectory      = "long.dcd"
topology        = "system.psf"
output_basename = "rmsd"
"#,
        )
        .unwrap();
        let input = MdtrajInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/mdtraj/bin/python");
        assert_eq!(input.trajectory, PathBuf::from("long.dcd"));
        assert_eq!(input.topology, PathBuf::from("system.psf"));
        assert_eq!(input.output_basename, "rmsd");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.R` from a
        // copy-paste off the Seurat sister, `.txt`); catch it at
        // parse time so the user gets a clear error before Python is
        // invoked.
        let d = tempdir("mdtraj-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mdtraj.analyze"

[bio.mdtraj]
script          = "analyze.R"
trajectory      = "md.xtc"
topology        = "md.pdb"
output_basename = "analysis"
"#,
        )
        .unwrap();
        let err = MdtrajInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
