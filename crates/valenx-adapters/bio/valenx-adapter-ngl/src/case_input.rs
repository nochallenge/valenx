//! `[bio.ngl]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "ngl.view"
//!
//! [bio.ngl]
//! script             = "viewer.py"
//! python             = "python3"          # optional, defaults to python3
//! # input_structure  = "structure.pdb"     # optional, omit when scripts fetch by PDB ID
//! output_basename    = "viewer"
//! ```
//!
//! NGL Viewer is the Rose lab's WebGL framework for molecular
//! visualization (MIT licensed). The core library is JavaScript,
//! widely embedded in Jupyter notebooks via the `nglview` Python
//! binding. This adapter wraps `nglview` so a user-authored
//! `viewer.py` can build an NGL view, export a self-contained HTML
//! page, render a still PNG, or dump the underlying state JSON, all
//! from the standard subprocess pattern.
//!
//! `input_structure` is optional: omit it for scripts that pull
//! coordinates by PDB ID over the network, or supply a path to an
//! existing `.pdb` / `.cif` / `.gro` file the script loads from
//! disk.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct NglInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input structure file (`.pdb`, `.cif`,
    /// `.gro`, etc.) the script loads. `None` means the script
    /// fetches coordinates by PDB ID over the network.
    pub input_structure: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.html` (self-contained viewer pages),
    /// `<basename>*.png` (rendered images), and `<basename>*.json`
    /// (NGL state JSON) into the workdir.
    pub output_basename: String,
}

impl NglInput {
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
            .and_then(|v| v.get("ngl"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.ngl] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.ngl].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ngl].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python tolerates
        // other extensions but `import nglview` workflows are
        // conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ngl].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_structure = match block.get("input_structure").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.ngl].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ngl].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_structure,
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
        let d = tempdir("ngl-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ngl.view"

[bio.ngl]
script          = "viewer.py"
output_basename = "viewer"
"#,
        )
        .unwrap();
        let input = NglInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("viewer.py"));
        assert_eq!(input.python, "python3");
        // No input_structure — script fetches by PDB ID.
        assert_eq!(input.input_structure, None);
        assert_eq!(input.output_basename, "viewer");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_structure() {
        // Pinned conda interpreter + a local PDB file the script
        // colors / annotates and renders to HTML.
        let d = tempdir("ngl-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ngl.view"

[bio.ngl]
script          = "render.py"
python          = "/opt/conda/envs/ngl/bin/python"
input_structure = "structure.pdb"
output_basename = "rendered"
"#,
        )
        .unwrap();
        let input = NglInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/ngl/bin/python");
        assert_eq!(input.input_structure, Some(PathBuf::from("structure.pdb")));
        assert_eq!(input.output_basename, "rendered");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.html`, `.json`
        // from a copy-paste off the output field); catch it at parse
        // time so the user gets a clear error before Python is
        // invoked.
        let d = tempdir("ngl-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ngl.view"

[bio.ngl]
script          = "viewer.html"
output_basename = "viewer"
"#,
        )
        .unwrap();
        let err = NglInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
