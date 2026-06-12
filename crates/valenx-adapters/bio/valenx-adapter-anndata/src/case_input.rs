//! `[bio.anndata]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "anndata.process"
//!
//! [bio.anndata]
//! script           = "process.py"
//! python           = "python3"           # optional, defaults to "python3"
//! # input_h5ad     = "matrix.h5ad"       # optional, omit to skip
//! output_basename  = "processed"
//! ```
//!
//! AnnData is the canonical Python single-cell data container library
//! (HDF5-backed `.h5ad` format that scanpy / scvi / etc. all read &
//! write). The adapter itself doesn't generate Python; the user
//! supplies a `process.py` that does `import anndata` and the actual
//! data work. We just spawn `python <script>` after staging the
//! script (and any optional `.h5ad` input) into the workdir.
//!
//! `input_h5ad` is optional: omit it for scripts that fetch /
//! synthesise their own data, or supply a path to an existing
//! `.h5ad` / `.h5` file the script reads via the staged filename.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct AnnDataInput {
    /// Path to the user-authored Python driver script (relative to the
    /// case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input `.h5ad` / `.h5` the script reads as
    /// its starting point. `None` means the script fetches /
    /// synthesises its own data.
    pub input_h5ad: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.h5ad` (AnnData files), `<basename>*.csv`
    /// (tables), and `<basename>*.png` (plots) into the workdir.
    pub output_basename: String,
}

impl AnnDataInput {
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
            .and_then(|v| v.get("anndata"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.anndata] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.anndata].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.anndata].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python
        // tolerates other extensions but `import anndata` workflows
        // are conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.anndata].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_h5ad = match block.get("input_h5ad").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.anndata].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.anndata].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_h5ad,
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
        let d = tempdir("anndata-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "anndata.process"

[bio.anndata]
script          = "process.py"
output_basename = "processed"
"#,
        )
        .unwrap();
        let input = AnnDataInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("process.py"));
        assert_eq!(input.python, "python3");
        // No input_h5ad — script fetches / synthesises its own.
        assert_eq!(input.input_h5ad, None);
        assert_eq!(input.output_basename, "processed");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_h5ad() {
        // Pinned Python interpreter + an input .h5ad the script reads.
        let d = tempdir("anndata-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "anndata.process"

[bio.anndata]
script          = "qc.py"
python          = "/opt/conda/envs/sc/bin/python"
input_h5ad      = "matrix.h5ad"
output_basename = "qc"
"#,
        )
        .unwrap();
        let input = AnnDataInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/sc/bin/python");
        assert_eq!(input.input_h5ad, Some(PathBuf::from("matrix.h5ad")));
        assert_eq!(input.output_basename, "qc");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.R` from a
        // copy-paste off the Seurat sister, `.txt`); catch it at
        // parse time so the user gets a clear error before Python is
        // invoked.
        let d = tempdir("anndata-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "anndata.process"

[bio.anndata]
script          = "process.R"
output_basename = "processed"
"#,
        )
        .unwrap();
        let err = AnnDataInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
