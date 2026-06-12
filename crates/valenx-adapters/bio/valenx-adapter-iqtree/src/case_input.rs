//! `[bio.iqtree]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "iqtree.tree"
//!
//! [bio.iqtree]
//! alignment  = "aln.fa"
//! model      = "MFP"          # optional, defaults to "MFP" (ModelFinder Plus)
//! bootstrap  = 1000           # optional, defaults to 1000 (UFBoot replicates); 0 = no bootstrap
//! threads    = "AUTO"         # optional, defaults to "AUTO"; numeric strings like "8" also accepted
//! prefix     = "iqtree_run"
//! extra_args = ["--quiet"]    # optional, defaults to []
//! ```
//!
//! `model` selects the substitution model. `MFP` is IQ-TREE's
//! ModelFinder Plus — it tries every model in the IQ-TREE catalogue
//! and picks the best by BIC. Specific names like `GTR+G`, `LG+I+G`,
//! or `WAG` are passed verbatim to IQ-TREE's `-m` flag.
//!
//! `bootstrap` enables ultrafast bootstrap (UFBoot) when > 0. Setting
//! to 0 skips bootstrapping (useful for fast topology-only runs).
//!
//! `threads` accepts the literal string `"AUTO"` (IQ-TREE auto-detects
//! the optimal core count) or a positive integer as a string. We keep
//! it as a string rather than a typed enum so `--prefix` propagation
//! is straightforward.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Returns true iff `s` is the literal `"AUTO"` (IQ-TREE's
/// auto-thread sentinel) or a positive integer in decimal form. We
/// match this by hand instead of pulling in a regex dep — IQ-TREE
/// accepts only this narrow set, so the predicate stays tiny.
pub fn is_valid_threads(s: &str) -> bool {
    if s == "AUTO" {
        return true;
    }
    !s.is_empty() && s.chars().all(|c| c.is_ascii_digit())
}

#[derive(Clone, Debug, PartialEq)]
pub struct IqTreeInput {
    pub alignment: PathBuf,
    pub model: String,
    pub bootstrap: u32,
    pub threads: String,
    pub prefix: String,
    pub extra_args: Vec<String>,
}

impl IqTreeInput {
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
            .and_then(|v| v.get("iqtree"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.iqtree] section",
                    case_toml.display()
                ))
            })?;

        let alignment_str = block
            .get("alignment")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.iqtree].alignment required (path to multi-FASTA / PHYLIP / NEXUS alignment)"
                ))
            })?;
        if alignment_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.iqtree].alignment must not be empty"
            )));
        }

        let model = match block.get("model") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.iqtree].model must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].model must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "MFP".to_string(),
        };

        let bootstrap = match block.get("bootstrap") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].bootstrap must be an integer"
                    ))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].bootstrap must be >= 0, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].bootstrap `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 1000,
        };

        let threads = match block.get("threads") {
            Some(v) => {
                // Accept either a string ("AUTO" or a numeric string) or
                // a TOML integer (we'll stringify it). IQ-TREE's CLI
                // takes both forms, so don't force the user to quote.
                if let Some(s) = v.as_str() {
                    if !is_valid_threads(s) {
                        return Err(AdapterError::Other(anyhow::anyhow!(
                            "[bio.iqtree].threads `{s}` invalid — \
                             expected \"AUTO\" or a positive integer"
                        )));
                    }
                    s.to_string()
                } else if let Some(i) = v.as_integer() {
                    if i < 1 {
                        return Err(AdapterError::Other(anyhow::anyhow!(
                            "[bio.iqtree].threads must be >= 1, got {i}"
                        )));
                    }
                    i.to_string()
                } else {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].threads must be a string \"AUTO\" or a positive integer"
                    )));
                }
            }
            None => "AUTO".to_string(),
        };

        let prefix = match block.get("prefix") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.iqtree].prefix must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].prefix must not be empty"
                    )));
                }
                s.to_string()
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.iqtree].prefix required (basename for IQ-TREE outputs)"
                )));
            }
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.iqtree].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.iqtree].extra_args entries must be strings"
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
    fn parses_minimal() {
        // Required fields only: alignment + prefix. Model defaults to
        // ModelFinder Plus, bootstrap to 1000, threads to AUTO.
        let d = tempdir("iqtree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "iqtree.tree"

[bio.iqtree]
alignment = "aln.fa"
prefix    = "run1"
"#,
        )
        .unwrap();
        let input = IqTreeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.alignment, PathBuf::from("aln.fa"));
        assert_eq!(input.prefix, "run1");
        assert_eq!(input.model, "MFP");
        assert_eq!(input.bootstrap, 1000);
        assert_eq!(input.threads, "AUTO");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_model_override() {
        // Explicit GTR+G model — common for nucleotide alignments
        // when the user already knows the right model.
        let d = tempdir("iqtree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "iqtree.tree"

[bio.iqtree]
alignment = "aln.fa"
model     = "GTR+G"
prefix    = "run1"
"#,
        )
        .unwrap();
        let input = IqTreeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model, "GTR+G");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_bootstrap() {
        // Lower bootstrap count + explicit numeric threads as string.
        let d = tempdir("iqtree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "iqtree.tree"

[bio.iqtree]
alignment  = "aln.fa"
bootstrap  = 100
threads    = "8"
prefix     = "run1"
extra_args = ["--quiet"]
"#,
        )
        .unwrap();
        let input = IqTreeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.bootstrap, 100);
        assert_eq!(input.threads, "8");
        assert_eq!(input.extra_args, vec!["--quiet".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_model() {
        // Empty `model = ""` would be silently dropped to default if we
        // didn't catch it; reject up front so the user sees the typo.
        let d = tempdir("iqtree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "iqtree.tree"

[bio.iqtree]
alignment = "aln.fa"
model     = ""
prefix    = "run1"
"#,
        )
        .unwrap();
        let err = IqTreeInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("model"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_invalid_threads() {
        // "auto" (lowercase) is NOT valid — IQ-TREE only recognises
        // the uppercase `AUTO` sentinel.
        let d = tempdir("iqtree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "iqtree.tree"

[bio.iqtree]
alignment = "aln.fa"
threads   = "auto"
prefix    = "run1"
"#,
        )
        .unwrap();
        let err = IqTreeInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
