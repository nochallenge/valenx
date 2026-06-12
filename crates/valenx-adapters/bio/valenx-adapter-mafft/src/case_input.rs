//! `[bio.mafft]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mafft.msa"
//!
//! [bio.mafft]
//! input      = "seqs.fa"
//! strategy   = "auto"             # optional, defaults to "auto"
//! threads    = 4                  # optional, defaults to 1; -1 = MAFFT auto-pick
//! extra_args = ["--reorder"]      # optional, defaults to []
//! ```
//!
//! `strategy` selects the MAFFT alignment algorithm:
//!
//! - `auto`   — let MAFFT pick based on input size (the
//!   recommended default; maps to the `--auto` flag)
//! - `linsi`  — L-INS-i, accuracy-first for sequences with one
//!   conserved domain (`--localpair --maxiterate 1000`)
//! - `ginsi`  — G-INS-i, accuracy-first for globally alignable
//!   sequences (`--globalpair --maxiterate 1000`)
//! - `einsi`  — E-INS-i, accuracy-first for sequences with multiple
//!   conserved domains separated by long unalignable regions
//! - `fftns1` — FFT-NS-1, fastest progressive alignment
//! - `fftns2` — FFT-NS-2, default progressive alignment
//! - `fftnsi` — FFT-NS-i, iterative refinement on the FFT-NS-2
//!   tree (the original MAFFT default)
//!
//! Threads accepts -1 (let MAFFT pick) but must not be 0 — that is
//! ambiguous in the MAFFT CLI and would silently degrade to single-
//! threaded execution.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical MAFFT strategy list. Module-public so the adapter can
/// surface the supported values to the UI without redefining them.
pub const SUPPORTED_STRATEGIES: &[&str] = &[
    "auto", "linsi", "ginsi", "einsi", "fftns1", "fftns2", "fftnsi",
];

#[derive(Clone, Debug, PartialEq)]
pub struct MafftInput {
    pub input: PathBuf,
    pub strategy: String,
    pub threads: i32,
    pub extra_args: Vec<String>,
}

impl MafftInput {
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
            .and_then(|v| v.get("mafft"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mafft] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.mafft].input required (path to multi-FASTA input)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mafft].input must not be empty"
            )));
        }

        let strategy = match block.get("strategy") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.mafft].strategy must be a string"))
                })?;
                if !SUPPORTED_STRATEGIES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mafft].strategy `{s}` not recognised — \
                         expected one of {SUPPORTED_STRATEGIES:?}"
                    )));
                }
                s.to_string()
            }
            None => "auto".to_string(),
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.mafft].threads must be an integer"))
                })?;
                if raw == 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mafft].threads must be >= 1 or -1 (auto), \
                         got 0"
                    )));
                }
                if raw < -1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mafft].threads must be >= 1 or -1 (auto), got {raw}"
                    )));
                }
                if raw > i32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mafft].threads `{raw}` exceeds i32::MAX"
                    )));
                }
                raw as i32
            }
            None => 1,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mafft].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.mafft].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input: PathBuf::from(input_str),
            strategy,
            threads,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn defaults_validate() {
        // Minimal config: only `input`. Everything else should fall
        // back to documented defaults — `auto` strategy, 1 thread,
        // no extras.
        let d = tempdir("mafft");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mafft.msa"

[bio.mafft]
input = "seqs.fa"
"#,
        )
        .unwrap();
        let input = MafftInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("seqs.fa"));
        assert_eq!(input.strategy, "auto");
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_strategy() {
        // A misspelt `"l-ins-i"` (which is the algorithm's familiar
        // name) should fail rather than be passed to MAFFT, which
        // doesn't recognise the dashed form on the CLI.
        let d = tempdir("mafft");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mafft.msa"

[bio.mafft]
input    = "seqs.fa"
strategy = "l-ins-i"
"#,
        )
        .unwrap();
        let err = MafftInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("linsi"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        // 0 is not a meaningful MAFFT thread count — it would silently
        // collapse to single-threaded; reject up front.
        let d = tempdir("mafft");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mafft.msa"

[bio.mafft]
input   = "seqs.fa"
threads = 0
"#,
        )
        .unwrap();
        let err = MafftInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_linsi_and_thread_auto() {
        // L-INS-i for accuracy with -1 threads (let MAFFT pick) and
        // a `--reorder` extra to sort the output by tree position.
        let d = tempdir("mafft");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mafft.msa"

[bio.mafft]
input      = "seqs.fa"
strategy   = "linsi"
threads    = -1
extra_args = ["--reorder"]
"#,
        )
        .unwrap();
        let input = MafftInput::from_case_dir(&d).unwrap();
        assert_eq!(input.strategy, "linsi");
        assert_eq!(input.threads, -1);
        assert_eq!(input.extra_args, vec!["--reorder".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }
}
