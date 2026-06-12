//! `[bio.muscle]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "muscle.msa"
//!
//! [bio.muscle]
//! input      = "seqs.fa"
//! mode       = "align"            # optional, defaults to "align"
//! threads    = 4                  # optional; if omitted MUSCLE auto-picks
//! extra_args = ["-perm", "abc"]   # optional, defaults to []
//! ```
//!
//! `mode` selects the MUSCLE 5 entry point:
//!
//! - `align`  — the default Probabilistic Consistency aligner;
//!   appropriate for inputs up to ~1k sequences.
//! - `super5` — the divide-and-conquer mode designed for very large
//!   inputs (10k+ sequences). Trades some accuracy for the ability
//!   to align datasets that `align` can't fit in memory.
//!
//! `threads` is `Option<u32>` rather than `u32` because MUSCLE 5
//! treats omission distinctly from any specific count: when no
//! `-threads` flag is passed it auto-selects based on the host
//! CPU. We preserve that distinction so a default-config case lets
//! MUSCLE pick.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical MUSCLE mode list. Module-public so the adapter can
/// surface the supported values to the UI without redefining them.
pub const SUPPORTED_MODES: &[&str] = &["align", "super5"];

#[derive(Clone, Debug, PartialEq)]
pub struct MuscleInput {
    pub input: PathBuf,
    pub mode: String,
    /// `None` -> let MUSCLE auto-pick the thread count.
    pub threads: Option<u32>,
    pub extra_args: Vec<String>,
}

impl MuscleInput {
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
            .and_then(|v| v.get("muscle"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.muscle] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.muscle].input required (path to multi-FASTA input)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.muscle].input must not be empty"
            )));
        }

        let mode = match block.get("mode") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.muscle].mode must be a string"))
                })?;
                if !SUPPORTED_MODES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.muscle].mode `{s}` not recognised — \
                         expected one of {SUPPORTED_MODES:?}"
                    )));
                }
                s.to_string()
            }
            None => "align".to_string(),
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.muscle].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.muscle].threads must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.muscle].threads `{raw}` exceeds u32::MAX"
                    )));
                }
                Some(raw as u32)
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.muscle].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.muscle].extra_args entries must be strings"
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
            mode,
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
        // Minimal config: only `input`. `mode` -> "align"; `threads`
        // -> None (let MUSCLE auto-pick); no extra args.
        let d = tempdir("muscle");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "muscle.msa"

[bio.muscle]
input = "seqs.fa"
"#,
        )
        .unwrap();
        let input = MuscleInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("seqs.fa"));
        assert_eq!(input.mode, "align");
        assert_eq!(input.threads, None);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_mode() {
        // Old MUSCLE 3.x had a `-quiet` flag and several modes; MUSCLE 5
        // collapsed everything to `align` / `super5`. Reject anything
        // outside that set fast.
        let d = tempdir("muscle");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "muscle.msa"

[bio.muscle]
input = "seqs.fa"
mode  = "fast"
"#,
        )
        .unwrap();
        let err = MuscleInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("super5"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_super5_and_threads() {
        // Large-input super5 path with explicit threading + a couple
        // of pass-through extras.
        let d = tempdir("muscle");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "muscle.msa"

[bio.muscle]
input      = "seqs.fa"
mode       = "super5"
threads    = 16
extra_args = ["-perm", "abc"]
"#,
        )
        .unwrap();
        let input = MuscleInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mode, "super5");
        assert_eq!(input.threads, Some(16));
        assert_eq!(
            input.extra_args,
            vec!["-perm".to_string(), "abc".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
