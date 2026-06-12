//! `[bio.forecast]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "forecast.predict"
//!
//! [bio.forecast]
//! script           = "predict.py"
//! python           = "python3"           # optional, defaults to python3
//! # input_fasta    = "target.fa"          # optional target sequence
//! output_basename  = "output"
//! ```
//!
//! FORECasT is the Sanger CRISPR/Cas9 indel pattern predictor (Allen
//! et al., Nature Biotechnology). FORECasT is published under the
//! SelfTarget name (github.com/felicityallen/SelfTarget); the Python
//! module is `selftarget`. The adapter doesn't reimplement prediction
//! logic; the user authors a `predict.py` that imports `selftarget`
//! and we spawn `python <script>` after staging script + optional
//! FASTA into the workdir.
//!
//! `input_fasta` is optional: omit when the script supplies sequences
//! inline, or supply a path to a `.fa` / `.fasta` the script reads as
//! the target window.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ForecastInput {
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
    /// `<basename>*.csv` (predicted indel tables), `<basename>*.txt`
    /// (summary), and any `*.log` files into the workdir.
    pub output_basename: String,
}

impl ForecastInput {
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
            .and_then(|v| v.get("forecast"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.forecast] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.forecast].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.forecast].script must not be empty"
            )));
        }
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.forecast].script `{script}` must end in `.py`"
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
                AdapterError::Other(anyhow::anyhow!("[bio.forecast].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.forecast].output_basename must not be empty"
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
        let d = tempdir("forecast-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "forecast.predict"

[bio.forecast]
script          = "predict.py"
output_basename = "output"
"#,
        )
        .unwrap();
        let input = ForecastInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("predict.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_fasta, None);
        assert_eq!(input.output_basename, "output");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_fasta() {
        // Pinned conda interpreter + a FASTA target window the script
        // predicts indel distributions over.
        let d = tempdir("forecast-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "forecast.predict"

[bio.forecast]
script          = "predict.py"
python          = "/opt/conda/envs/forecast/bin/python"
input_fasta     = "target.fa"
output_basename = "predictions"
"#,
        )
        .unwrap();
        let input = ForecastInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/forecast/bin/python");
        assert_eq!(input.input_fasta, Some(PathBuf::from("target.fa")));
        assert_eq!(input.output_basename, "predictions");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        let d = tempdir("forecast-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "forecast.predict"

[bio.forecast]
script          = "predict.fa"
output_basename = "output"
"#,
        )
        .unwrap();
        let err = ForecastInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
