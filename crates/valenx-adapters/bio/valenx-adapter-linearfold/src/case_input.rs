//! `[bio.linearfold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "linearfold.fold"
//!
//! [bio.linearfold]
//! sequence        = "rna.fa"
//! output_basename = "fold"
//! model           = "C"            # optional, defaults to "C"; "C" (CONTRAfold) | "V" (ViennaRNA)
//! beam_size       = 100            # optional, defaults to 100
//! extra_args      = []             # optional, defaults to []
//! ```
//!
//! LinearFold is the Baidu/OSU linear-time RNA secondary-structure
//! folder — the first algorithm to break the cubic-time barrier of
//! classic dynamic-programming folders, making genome-scale
//! single-sequence folding tractable. The CLI reads the sequence from
//! **stdin** and writes the predicted structure to **stdout**; the
//! adapter pipes the input file in and redirects stdout to a file
//! named after `output_basename`.
//!
//! `model` selects the scoring function:
//! - `"C"` — CONTRAfold parameters (the default; ML-trained scoring)
//! - `"V"` — ViennaRNA parameters (classical Turner thermodynamics)
//!
//! `beam_size` is the beam-search width: larger values approach
//! optimality, smaller values trade accuracy for speed. Default 100
//! matches the upstream recommendation.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical LinearFold scoring-model list. Module-public so the
/// adapter can surface the supported values to the UI without
/// redefining them.
pub const SUPPORTED_MODELS: &[&str] = &["C", "V"];

#[derive(Clone, Debug, PartialEq)]
pub struct LinearFoldInput {
    /// Path to the single-sequence input. Accepts FASTA (`.fa` /
    /// `.fasta`) or a plain text file containing the sequence on a
    /// single line. LinearFold reads from stdin; the adapter pipes
    /// this file in. Relative paths resolve against the case directory.
    pub sequence: PathBuf,
    /// Filename stem for the structure output. The redirected stdout
    /// lands at `<output_basename>.txt` in the workdir.
    pub output_basename: String,
    /// Scoring-model selector. `"C"` (CONTRAfold ML parameters,
    /// default) or `"V"` (ViennaRNA Turner thermodynamics). The
    /// upstream binary maps these to its `-V` / `-C` flags.
    pub model: String,
    /// Beam-search width. Larger values approach optimality at
    /// quadratic cost; smaller values trade accuracy for speed.
    /// Default 100.
    pub beam_size: u32,
    /// Additional CLI arguments appended to the linearfold invocation.
    /// Useful for `--sharpturn`, `--no-sharpturn`, or any future flags
    /// the upstream tool grows.
    pub extra_args: Vec<String>,
}

impl LinearFoldInput {
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
            .and_then(|v| v.get("linearfold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.linearfold] section",
                    case_toml.display()
                ))
            })?;

        let sequence = block
            .get("sequence")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.linearfold].sequence required (path to FASTA or plain-text \
                     sequence input)"
                ))
            })?;
        if sequence.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.linearfold].sequence must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.linearfold].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.linearfold].output_basename must not be empty"
            )));
        }

        let model = match block.get("model") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.linearfold].model must be a string"))
                })?;
                if !SUPPORTED_MODELS.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.linearfold].model `{s}` not recognised — \
                         expected one of {SUPPORTED_MODELS:?}"
                    )));
                }
                s.to_string()
            }
            None => "C".to_string(),
        };

        // beam_size accepts non-negative integers; 0 is a documented
        // upstream value meaning "exact (no beam pruning)".
        let beam_size = match block.get("beam_size") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.linearfold].beam_size must be an integer"
                    ))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.linearfold].beam_size must be >= 0, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.linearfold].beam_size `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 100,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.linearfold].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.linearfold].extra_args entries must be strings"
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
            model,
            beam_size,
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
        // Defaults: model "C" (CONTRAfold), beam_size 100, no extras.
        let d = tempdir("linearfold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "linearfold.fold"

[bio.linearfold]
sequence        = "rna.fa"
output_basename = "fold"
"#,
        )
        .unwrap();
        let input = LinearFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.sequence, PathBuf::from("rna.fa"));
        assert_eq!(input.output_basename, "fold");
        assert_eq!(input.model, "C");
        assert_eq!(input.beam_size, 100);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_vienna_model() {
        // Vienna parameters with a tighter beam for speed plus an
        // extra knob — exercises the model swap, beam override, and
        // extras passthrough together.
        let d = tempdir("linearfold-full");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "linearfold.fold"

[bio.linearfold]
sequence        = "rna.fa"
output_basename = "fold-vienna"
model           = "V"
beam_size       = 50
extra_args      = ["--no-sharpturn"]
"#,
        )
        .unwrap();
        let input = LinearFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model, "V");
        assert_eq!(input.beam_size, 50);
        assert_eq!(input.extra_args, vec!["--no-sharpturn".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_model() {
        // Models other than "C" / "V" are silently passed by argparse
        // but yield nonsense output; reject up front so the failure
        // is fast and obvious.
        let d = tempdir("linearfold-badmodel");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "linearfold.fold"

[bio.linearfold]
sequence        = "rna.fa"
output_basename = "fold"
model           = "T"
"#,
        )
        .unwrap();
        let err = LinearFoldInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
