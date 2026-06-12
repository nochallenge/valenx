//! `[bio.crispritz]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "crispritz.search"
//!
//! [bio.crispritz]
//! script           = "search.py"
//! python           = "python3"           # optional, defaults to python3
//! # input_fasta    = "target.fa"          # optional guide / target sequence
//! output_basename  = "output"
//! ```
//!
//! CRISPRitz is the Pinello lab off-target search tool (Cancellieri et
//! al., Bioinformatics): given a guide RNA and a reference genome it
//! enumerates putative off-target sites supporting mismatches, DNA /
//! RNA bulges, and PAM variants. The adapter doesn't reimplement
//! search logic; the user authors a `search.py` that imports the
//! upstream `crispritz` Python package
//! (github.com/pinellolab/CRISPRitz) and we spawn `python <script>`
//! after staging script + optional FASTA into the workdir.
//!
//! `input_fasta` is optional: omit when the script supplies guides /
//! targets inline, or supply a path to a `.fa` / `.fasta` the script
//! reads as the input.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CrispritzInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input FASTA the script reads as the guide
    /// or target. `None` means the script supplies its own input.
    pub input_fasta: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.txt` (off-target tables), `<basename>*.bed`
    /// (BED files), and any `*.log` files into the workdir.
    pub output_basename: String,
}

impl CrispritzInput {
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
            .and_then(|v| v.get("crispritz"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.crispritz] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.crispritz].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispritz].script must not be empty"
            )));
        }
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispritz].script `{script}` must end in `.py`"
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
                AdapterError::Other(anyhow::anyhow!("[bio.crispritz].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispritz].output_basename must not be empty"
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
        let d = tempdir("crispritz-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispritz.search"

[bio.crispritz]
script          = "search.py"
output_basename = "output"
"#,
        )
        .unwrap();
        let input = CrispritzInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("search.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_fasta, None);
        assert_eq!(input.output_basename, "output");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_fasta() {
        // Pinned conda interpreter + a FASTA guide file the script
        // searches off-targets for.
        let d = tempdir("crispritz-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispritz.search"

[bio.crispritz]
script          = "search.py"
python          = "/opt/conda/envs/crispritz/bin/python"
input_fasta     = "guides.fa"
output_basename = "offtargets"
"#,
        )
        .unwrap();
        let input = CrispritzInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/crispritz/bin/python");
        assert_eq!(input.input_fasta, Some(PathBuf::from("guides.fa")));
        assert_eq!(input.output_basename, "offtargets");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        let d = tempdir("crispritz-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispritz.search"

[bio.crispritz]
script          = "search.fa"
output_basename = "output"
"#,
        )
        .unwrap();
        let err = CrispritzInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
