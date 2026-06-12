//! `[bio.raxml-ng]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "raxml-ng.tree"
//!
//! [bio.raxml-ng]
//! alignment  = "aln.fa"
//! model      = "GTR+G"
//! mode       = "search"          # "search" | "all" | "bootstrap"
//! bootstrap  = 100               # 0 OK for "search"; required >= 1 for "all" / "bootstrap"
//! threads    = 4                 # optional, defaults to 1
//! prefix     = "raxml_run"
//! extra_args = ["--seed", "42"]  # optional
//! ```
//!
//! `mode` selects the RAxML-NG operation:
//!
//! - `search`    — ML tree search only (`--search`); no bootstrap
//! - `bootstrap` — bootstrap replicates only (`--bootstrap`)
//! - `all`       — search + bootstrap + bootstrap-mapping (`--all`)
//!
//! `model` is required (no auto-detection in RAxML-NG, unlike
//! IQ-TREE's MFP). Common choices: `GTR+G` for nucleotides, `LG+G`
//! for amino-acids.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical RAxML-NG mode list. Module-public so the adapter can
/// surface the supported values to the UI.
pub const SUPPORTED_MODES: &[&str] = &["search", "all", "bootstrap"];

#[derive(Clone, Debug, PartialEq)]
pub struct RaxmlNgInput {
    pub alignment: PathBuf,
    pub model: String,
    pub mode: String,
    pub bootstrap: u32,
    pub threads: u32,
    pub prefix: String,
    pub extra_args: Vec<String>,
}

impl RaxmlNgInput {
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
            .and_then(|v| v.get("raxml-ng"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.raxml-ng] section",
                    case_toml.display()
                ))
            })?;

        let alignment_str = block
            .get("alignment")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.raxml-ng].alignment required (path to multi-FASTA / PHYLIP alignment)"
                ))
            })?;
        if alignment_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.raxml-ng].alignment must not be empty"
            )));
        }

        let model = block.get("model").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.raxml-ng].model required (e.g. \"GTR+G\" or \"LG+G\")"
            ))
        })?;
        if model.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.raxml-ng].model must not be empty"
            )));
        }
        let model = model.to_string();

        let mode = match block.get("mode") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.raxml-ng].mode must be a string"))
                })?;
                if !SUPPORTED_MODES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].mode `{s}` not recognised — \
                         expected one of {SUPPORTED_MODES:?}"
                    )));
                }
                s.to_string()
            }
            None => "search".to_string(),
        };

        let bootstrap = match block.get("bootstrap") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].bootstrap must be an integer"
                    ))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].bootstrap must be >= 0, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].bootstrap `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 0,
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].threads must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let prefix = match block.get("prefix") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.raxml-ng].prefix must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].prefix must not be empty"
                    )));
                }
                s.to_string()
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.raxml-ng].prefix required (basename for RAxML-NG outputs)"
                )));
            }
        };

        // Cross-field validation: bootstrap and all-mode runs both
        // need a positive replicate count or RAxML-NG will reject the
        // invocation. Catch up front with a more pointed error.
        if (mode == "bootstrap" || mode == "all") && bootstrap < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.raxml-ng] mode = \"{mode}\" requires bootstrap >= 1"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.raxml-ng].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.raxml-ng].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            alignment: PathBuf::from(alignment_str),
            model,
            mode,
            bootstrap,
            threads,
            prefix,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_search_minimal() {
        // Minimal: alignment + model + prefix; mode defaults to
        // search, bootstrap to 0 (unused for search), threads to 1.
        let d = tempdir("raxml-ng-raxml");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "raxml-ng.tree"

[bio.raxml-ng]
alignment = "aln.fa"
model     = "GTR+G"
prefix    = "run1"
"#,
        )
        .unwrap();
        let input = RaxmlNgInput::from_case_dir(&d).unwrap();
        assert_eq!(input.alignment, PathBuf::from("aln.fa"));
        assert_eq!(input.model, "GTR+G");
        assert_eq!(input.mode, "search");
        assert_eq!(input.bootstrap, 0);
        assert_eq!(input.threads, 1);
        assert_eq!(input.prefix, "run1");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_all_with_bootstrap() {
        // Full all-mode run with explicit bootstrap and threads.
        let d = tempdir("raxml-ng-raxml");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "raxml-ng.tree"

[bio.raxml-ng]
alignment  = "aln.fa"
model      = "LG+G"
mode       = "all"
bootstrap  = 200
threads    = 8
prefix     = "run1"
extra_args = ["--seed", "42"]
"#,
        )
        .unwrap();
        let input = RaxmlNgInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mode, "all");
        assert_eq!(input.bootstrap, 200);
        assert_eq!(input.threads, 8);
        assert_eq!(
            input.extra_args,
            vec!["--seed".to_string(), "42".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_mode() {
        // `evaluate` is an actual RAxML-NG mode (likelihood eval on
        // a fixed tree) but not one this adapter supports; reject.
        let d = tempdir("raxml-ng-raxml");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "raxml-ng.tree"

[bio.raxml-ng]
alignment = "aln.fa"
model     = "GTR+G"
mode      = "evaluate"
prefix    = "run1"
"#,
        )
        .unwrap();
        let err = RaxmlNgInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("search"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        let d = tempdir("raxml-ng-raxml");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "raxml-ng.tree"

[bio.raxml-ng]
alignment = "aln.fa"
model     = "GTR+G"
threads   = 0
prefix    = "run1"
"#,
        )
        .unwrap();
        let err = RaxmlNgInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_bootstrap_mode_without_replicates() {
        // mode = "bootstrap" without bootstrap >= 1 is a no-op that
        // RAxML-NG would reject anyway; surface the error here so
        // the user catches it before the run starts.
        let d = tempdir("raxml-ng-raxml");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "raxml-ng.tree"

[bio.raxml-ng]
alignment = "aln.fa"
model     = "GTR+G"
mode      = "bootstrap"
prefix    = "run1"
"#,
        )
        .unwrap();
        let err = RaxmlNgInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("bootstrap"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
