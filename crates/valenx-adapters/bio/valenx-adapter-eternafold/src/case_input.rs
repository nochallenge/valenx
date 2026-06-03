//! `[bio.eternafold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "eternafold.fold"
//!
//! [bio.eternafold]
//! script           = "fold.py"
//! python           = "python3"           # optional, defaults to python3
//! # input_fasta    = "rna.fa"             # optional
//! output_basename  = "fold"
//! ```
//!
//! EternaFold is the Eterna Project's ML-aware RNA secondary-structure
//! folder, trained on community-generated structure data from the
//! Eterna game. The reference C++ binary lives inside the upstream
//! repo; in practice most users access EternaFold via the
//! [`arnie`](https://github.com/DasLab/arnie) Python wrapper, which
//! bundles EternaFold alongside ViennaRNA, NUPACK, and several other
//! folders behind a single `bp_matrix(...)` / `mfe(...)` API. The
//! adapter targets that workflow: the user authors a `fold.py` that
//! does `from arnie.mfe import mfe; mfe(seq, package='eternafold')`
//! and the actual folding logic. We just spawn `python <script>`
//! after staging the script (and any optional `.fa` input) into the
//! workdir.
//!
//! `input_fasta` is optional: omit it for scripts that hardcode the
//! sequence or pull it from another source, or supply a path to an
//! existing `.fa` / `.fasta` file the script reads as its template.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct EternaFoldInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input FASTA the script reads. `None` means
    /// the script hardcodes the sequence or fetches it from another
    /// source.
    pub input_fasta: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.ct` (connect-table), `<basename>*.dot` (dot-bracket
    /// notation), and `<basename>*.csv` (MEA / probabilities) into
    /// the workdir.
    pub output_basename: String,
}

impl EternaFoldInput {
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
            .and_then(|v| v.get("eternafold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.eternafold] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.eternafold].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eternafold].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python
        // tolerates other extensions but `import arnie` workflows are
        // conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eternafold].script `{script}` must end in `.py`"
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
                AdapterError::Other(anyhow::anyhow!("[bio.eternafold].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eternafold].output_basename must not be empty"
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
        let d = tempdir("eternafold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eternafold.fold"

[bio.eternafold]
script          = "fold.py"
output_basename = "fold"
"#,
        )
        .unwrap();
        let input = EternaFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("fold.py"));
        assert_eq!(input.python, "python3");
        // No input_fasta — script hardcodes the sequence.
        assert_eq!(input.input_fasta, None);
        assert_eq!(input.output_basename, "fold");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_fasta() {
        // Pinned conda interpreter + an existing FASTA the script
        // folds via arnie's eternafold backend.
        let d = tempdir("eternafold-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eternafold.fold"

[bio.eternafold]
script          = "batch_fold.py"
python          = "/opt/conda/envs/rna/bin/python"
input_fasta     = "rna.fa"
output_basename = "batch"
"#,
        )
        .unwrap();
        let input = EternaFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/rna/bin/python");
        assert_eq!(input.input_fasta, Some(PathBuf::from("rna.fa")));
        assert_eq!(input.output_basename, "batch");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.fa` from a
        // copy-paste off the input field); catch it at parse time so
        // the user gets a clear error before Python is invoked.
        let d = tempdir("eternafold-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eternafold.fold"

[bio.eternafold]
script          = "fold.fa"
output_basename = "fold"
"#,
        )
        .unwrap();
        let err = EternaFoldInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
