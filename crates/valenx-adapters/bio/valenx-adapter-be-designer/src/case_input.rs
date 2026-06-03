//! `[bio.be-designer]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "be-designer.design"
//!
//! [bio.be-designer]
//! script           = "design.py"
//! python           = "python3"           # optional, defaults to python3
//! # input_fasta    = "target.fa"          # optional target sequence
//! output_basename  = "output"
//! ```
//!
//! BE-Designer is the rgenome.net base-editor guide design tool: given
//! a target DNA window it enumerates sgRNAs whose editing window
//! covers the desired position for the supplied base editor (CBE / ABE
//! family). The adapter itself doesn't talk to the rgenome web API;
//! the user authors a `design.py` that invokes the local `bedesigner`
//! Python package (or shell-outs to the upstream pipeline) and we
//! spawn `python <script>` after staging script + optional FASTA into
//! the workdir.
//!
//! `input_fasta` is optional: omit when the script generates / fetches
//! its own target sequences, or supply a path to a `.fa` / `.fasta`
//! the script reads as the target window.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct BeDesignerInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input FASTA the script reads as the target
    /// window. `None` means the script supplies its own target.
    pub input_fasta: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.csv` (guide tables), `<basename>*.fasta`
    /// (designed sequence exports), and any `*.log` files into the
    /// workdir.
    pub output_basename: String,
}

impl BeDesignerInput {
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
            .and_then(|v| v.get("be-designer"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.be-designer] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.be-designer].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.be-designer].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python tolerates
        // other extensions but `import bedesigner` workflows are
        // conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.be-designer].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_fasta = match block.get("input_fasta").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.be-designer].output_basename required"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.be-designer].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_fasta,
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
        let d = tempdir("be-designer-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "be-designer.design"

[bio.be-designer]
script          = "design.py"
output_basename = "output"
"#,
        )
        .unwrap();
        let input = BeDesignerInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.python, "python3");
        // No input_fasta — script supplies its own target.
        assert_eq!(input.input_fasta, None);
        assert_eq!(input.output_basename, "output");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_fasta() {
        // Pinned conda interpreter + an existing FASTA target window
        // the script enumerates sgRNAs across.
        let d = tempdir("be-designer-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "be-designer.design"

[bio.be-designer]
script          = "design.py"
python          = "/opt/conda/envs/bedesigner/bin/python"
input_fasta     = "target.fa"
output_basename = "guides"
"#,
        )
        .unwrap();
        let input = BeDesignerInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/bedesigner/bin/python");
        assert_eq!(input.input_fasta, Some(PathBuf::from("target.fa")));
        assert_eq!(input.output_basename, "guides");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.fasta` from a
        // copy-paste off the input field); catch it at parse time so
        // the user gets a clear error before Python is invoked.
        let d = tempdir("be-designer-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "be-designer.design"

[bio.be-designer]
script          = "design.fa"
output_basename = "output"
"#,
        )
        .unwrap();
        let err = BeDesignerInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
