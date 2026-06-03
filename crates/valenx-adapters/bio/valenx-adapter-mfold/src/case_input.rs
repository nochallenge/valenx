//! `[bio.mfold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mfold.fold"
//!
//! [bio.mfold]
//! sequence        = "rna.seq"
//! output_basename = "fold"
//! temperature     = 37.0          # optional, defaults to 37.0 (Celsius)
//! extra_args      = []            # optional, defaults to []
//! ```
//!
//! mfold (and its UNAFold successor) is Michael Zuker's classic RNA
//! secondary-structure folder. The CLI uses `KEY=VALUE` argument syntax
//! (e.g. `mfold SEQ=rna.seq NA=RNA T=37`) — every knob is a positional
//! `KEY=VALUE` token rather than a `--flag value` pair.
//!
//! `sequence` accepts FASTA (`.fa` / `.fasta`) and the legacy mfold
//! single-sequence `.seq` format (one sequence with optional metadata
//! line). The file is read in place — relative paths resolve against
//! the case directory but the file is NOT staged into the workdir.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct MfoldInput {
    /// Path to the single-sequence input. Accepts FASTA (`.fa` /
    /// `.fasta`) or the legacy mfold `.seq` format. Relative paths
    /// resolve against the case directory; the file is read in place.
    pub sequence: PathBuf,
    /// Filename stem for outputs. mfold writes
    /// `<basename>.ct` (connect-table), `<basename>*.ps` / `*.pdf`
    /// (structure plots), and `<basename>*.out` (run log) into the
    /// workdir.
    pub output_basename: String,
    /// Folding temperature in degrees Celsius. mfold's classic
    /// thermodynamic parameters target 37 °C; warmer / cooler values
    /// shift base-pair stabilities. Default 37.0.
    pub temperature: f64,
    /// Additional CLI arguments appended to the mfold invocation.
    /// Useful for `MAX=<n>` (max number of structures), `WIN=<n>`
    /// (window for hairpin search), `RUN_TYPE=html`, or any `KEY=VALUE`
    /// knob the upstream tool grows.
    pub extra_args: Vec<String>,
}

impl MfoldInput {
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
            .and_then(|v| v.get("mfold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mfold] section",
                    case_toml.display()
                ))
            })?;

        let sequence = block
            .get("sequence")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.mfold].sequence required (path to FASTA or .seq input)"
                ))
            })?;
        if sequence.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mfold].sequence must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.mfold].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mfold].output_basename must not be empty"
            )));
        }

        // Folding temperature in Celsius. Admit integers
        // (e.g. `temperature = 37`) and validate finite.
        let temperature = match block.get("temperature") {
            Some(v) => {
                let raw = if let Some(f) = v.as_float() {
                    f
                } else if let Some(i) = v.as_integer() {
                    i as f64
                } else {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mfold].temperature must be a number"
                    )));
                };
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mfold].temperature must be finite, got {raw}"
                    )));
                }
                raw
            }
            None => 37.0,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mfold].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.mfold].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            sequence: PathBuf::from(sequence),
            output_basename: output_basename.to_string(),
            temperature,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        // Only the two required keys: sequence and output_basename.
        // Defaults: temperature 37.0 °C, no extras.
        let d = tempdir("mfold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mfold.fold"

[bio.mfold]
sequence        = "rna.seq"
output_basename = "fold"
"#,
        )
        .unwrap();
        let input = MfoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.sequence, PathBuf::from("rna.seq"));
        assert_eq!(input.output_basename, "fold");
        assert_eq!(input.temperature, 37.0);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_overrides() {
        // Cooler folding (mfold accepts T < 37) plus a couple of
        // KEY=VALUE extras. Validates the f64 path on a non-integer
        // temperature and the extra_args passthrough.
        let d = tempdir("mfold-full");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mfold.fold"

[bio.mfold]
sequence        = "tRNA.fa"
output_basename = "trna-fold"
temperature     = 25.5
extra_args      = ["MAX=10", "WIN=5"]
"#,
        )
        .unwrap();
        let input = MfoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.sequence, PathBuf::from("tRNA.fa"));
        assert_eq!(input.output_basename, "trna-fold");
        assert_eq!(input.temperature, 25.5);
        assert_eq!(
            input.extra_args,
            vec!["MAX=10".to_string(), "WIN=5".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_sequence() {
        // An empty sequence path is the most common typo; catch it at
        // parse time so the user gets a clear error rather than a
        // confusing mfold "file not found" later.
        let d = tempdir("mfold-empty");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mfold.fold"

[bio.mfold]
sequence        = ""
output_basename = "fold"
"#,
        )
        .unwrap();
        let err = MfoldInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("sequence"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
