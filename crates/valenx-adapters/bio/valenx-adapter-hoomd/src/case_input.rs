//! `[bio.hoomd]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "hoomd.simulate"
//!
//! [bio.hoomd]
//! script           = "production.py"
//! python           = "python3"           # optional, defaults to "python3"
//! output_basename  = "production"
//! ```
//!
//! HOOMD-blue is the Glotzer lab's GPU-native particle molecular
//! dynamics engine. From v3 onward there is no native CLI — every
//! simulation is driven by a user-authored Python script that does
//! `import hoomd` and constructs / runs the simulation. The adapter
//! therefore doesn't generate Python; it stages the user-supplied
//! script into the workdir, drops a flat `valenx_params.json` next
//! to it, and spawns `python <script>`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct HoomdInput {
    /// Path to the user-authored HOOMD-blue Python driver script
    /// (relative to the case directory, or absolute). Must end in
    /// `.py` (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs (HOOMD-blue is
    /// often installed inside a dedicated conda env with the matching
    /// CUDA toolkit).
    pub python: String,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.gsd` (HOOMD trajectories) and `<basename>*.h5`
    /// (HDF5 outputs) into the workdir.
    pub output_basename: String,
}

impl HoomdInput {
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
            .and_then(|v| v.get("hoomd"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.hoomd] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.hoomd].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hoomd].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). HOOMD-blue v3+
        // is purely Python-driven; flagging the wrong extension up
        // front saves a confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hoomd].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.hoomd].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hoomd].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
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
        let d = tempdir("hoomd-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hoomd.simulate"

[bio.hoomd]
script          = "production.py"
output_basename = "production"
"#,
        )
        .unwrap();
        let input = HoomdInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("production.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.output_basename, "production");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_pinned_python() {
        // Pinned interpreter — typical when HOOMD-blue lives in a
        // dedicated conda env matched to a specific CUDA toolkit.
        let d = tempdir("hoomd-pinned");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hoomd.simulate"

[bio.hoomd]
script          = "lj_fluid.py"
python          = "/opt/conda/envs/hoomd/bin/python"
output_basename = "lj_fluid"
"#,
        )
        .unwrap();
        let input = HoomdInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("lj_fluid.py"));
        assert_eq!(input.python, "/opt/conda/envs/hoomd/bin/python");
        assert_eq!(input.output_basename, "lj_fluid");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.in` from a
        // copy-paste off the LAMMPS sister, `.txt`); catch it at
        // parse time so the user gets a clear error before Python is
        // invoked.
        let d = tempdir("hoomd-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hoomd.simulate"

[bio.hoomd]
script          = "production.in"
output_basename = "production"
"#,
        )
        .unwrap();
        let err = HoomdInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
