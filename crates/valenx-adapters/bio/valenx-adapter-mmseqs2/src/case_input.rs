//! `[bio.mmseqs2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mmseqs2.search"
//!
//! [bio.mmseqs2]
//! action      = "easy-search"          # one of "easy-search" | "easy-cluster" | "easy-linsearch"
//! query       = "query.fa"
//! target      = "uniref90.fa"          # required for easy-search / easy-linsearch; ignored for easy-cluster
//! output      = "hits.m8"
//! sensitivity = 7.5                    # 1.0..=7.5, default 7.5
//! threads     = 8                      # optional, defaults to 1
//! extra_args  = ["-e", "1e-3"]         # optional, defaults to []
//! ```
//!
//! `action` selects which MMseqs2 high-level workflow the adapter
//! wraps:
//!
//! - `easy-search`    — exhaustive iterative profile search of `query`
//!   against `target`. Output is a BLAST-format-8 hit table.
//! - `easy-linsearch` — linear-time prefilter variant of `easy-search`,
//!   trades a little sensitivity for ~10x speed-up. Same I/O shape.
//! - `easy-cluster`   — single-pass clustering of `query` (no separate
//!   target). Output is the cluster representative table; the
//!   `_all_seqs.fasta`, `_cluster.tsv`, and `_rep_seq.fasta` files are
//!   sister artifacts MMseqs2 emits beside the basename.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical MMseqs2 action list. Module-public so the UI can surface
/// the supported values without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["easy-search", "easy-cluster", "easy-linsearch"];

/// Inclusive sensitivity bounds for `--s`. MMseqs2 internally clamps
/// the prefilter to [1.0, 7.5]; values outside the range either fall
/// back to the default or raise an error at runtime, so we reject them
/// up front here.
pub const SENSITIVITY_MIN: f64 = 1.0;
pub const SENSITIVITY_MAX: f64 = 7.5;

#[derive(Clone, Debug, PartialEq)]
pub struct Mmseqs2Input {
    pub action: String,
    pub query: PathBuf,
    pub target: Option<PathBuf>,
    pub output: PathBuf,
    pub sensitivity: f64,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl Mmseqs2Input {
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
            .and_then(|v| v.get("mmseqs2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mmseqs2] section",
                    case_toml.display()
                ))
            })?;

        let action = block
            .get("action")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.mmseqs2].action required (one of {SUPPORTED_ACTIONS:?})"
                ))
            })?;
        if !SUPPORTED_ACTIONS.contains(&action) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mmseqs2].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        let query_str = block
            .get("query")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.mmseqs2].query required")))?;
        if query_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mmseqs2].query must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.mmseqs2].output required")))?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mmseqs2].output must not be empty"
            )));
        }

        // `target` is required for easy-search / easy-linsearch; the
        // clustering action takes a single FASTA via `query` and has
        // no `target` slot at all.
        let target = block
            .get("target")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        match action {
            "easy-search" | "easy-linsearch" => {
                if target.is_none() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mmseqs2].target required for action `{action}`"
                    )));
                }
            }
            "easy-cluster" => {
                // No-op: target is meaningless for clustering. We
                // accept an inadvertent value silently rather than
                // erroring — MMseqs2 itself ignores it.
            }
            _ => unreachable!("action validated against SUPPORTED_ACTIONS above"),
        }

        let sensitivity = match block.get("sensitivity") {
            Some(v) => {
                let raw = v.as_float().or_else(|| v.as_integer().map(|i| i as f64));
                let raw = raw.ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mmseqs2].sensitivity must be a number"
                    ))
                })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mmseqs2].sensitivity must be finite, got {raw}"
                    )));
                }
                if !(SENSITIVITY_MIN..=SENSITIVITY_MAX).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mmseqs2].sensitivity must be in \
                         [{SENSITIVITY_MIN}, {SENSITIVITY_MAX}], got {raw}"
                    )));
                }
                raw
            }
            None => SENSITIVITY_MAX,
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.mmseqs2].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mmseqs2].threads must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mmseqs2].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.mmseqs2].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            action: action.to_string(),
            query: PathBuf::from(query_str),
            target,
            output: PathBuf::from(output_str),
            sensitivity,
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
    fn parses_easy_search_minimal() {
        // easy-search needs target; defaults to sensitivity=7.5,
        // threads=1, no extras.
        let d = tempdir("mmseqs2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mmseqs2.search"

[bio.mmseqs2]
action = "easy-search"
query  = "query.fa"
target = "uniref90.fa"
output = "hits.m8"
"#,
        )
        .unwrap();
        let input = Mmseqs2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "easy-search");
        assert_eq!(input.query, PathBuf::from("query.fa"));
        assert_eq!(input.target, Some(PathBuf::from("uniref90.fa")));
        assert_eq!(input.output, PathBuf::from("hits.m8"));
        assert!((input.sensitivity - 7.5).abs() < f64::EPSILON);
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_easy_cluster_no_target() {
        // Clustering is a one-FASTA workflow; `target` is meaningless
        // and must not be required.
        let d = tempdir("mmseqs2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mmseqs2.cluster"

[bio.mmseqs2]
action      = "easy-cluster"
query       = "all_proteins.fa"
output      = "clusters"
sensitivity = 5.7
threads     = 16
extra_args  = ["--min-seq-id", "0.5"]
"#,
        )
        .unwrap();
        let input = Mmseqs2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "easy-cluster");
        assert_eq!(input.query, PathBuf::from("all_proteins.fa"));
        assert_eq!(input.target, None);
        assert!((input.sensitivity - 5.7).abs() < 1e-9);
        assert_eq!(input.threads, 16);
        assert_eq!(
            input.extra_args,
            vec!["--min-seq-id".to_string(), "0.5".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `mmseqs createdb` is a real command but not one this adapter
        // wraps — must be rejected up front.
        let d = tempdir("mmseqs2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mmseqs2.createdb"

[bio.mmseqs2]
action = "createdb"
query  = "query.fa"
output = "queryDB"
"#,
        )
        .unwrap();
        let err = Mmseqs2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("easy-search"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_search_without_target() {
        // easy-search without a target DB is a useless invocation —
        // catch it before MMseqs2 itself bails.
        let d = tempdir("mmseqs2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mmseqs2.search"

[bio.mmseqs2]
action = "easy-search"
query  = "query.fa"
output = "hits.m8"
"#,
        )
        .unwrap();
        let err = Mmseqs2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("target required"), "msg: {msg}");
        assert!(msg.contains("easy-search"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_sensitivity_above_75() {
        // MMseqs2 caps sensitivity at 7.5; anything above is silently
        // clamped or rejected at runtime depending on version, so we
        // pre-empt with an InvalidCase here.
        let d = tempdir("mmseqs2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mmseqs2.search"

[bio.mmseqs2]
action      = "easy-search"
query       = "query.fa"
target      = "db.fa"
output      = "hits.m8"
sensitivity = 9.0
"#,
        )
        .unwrap();
        let err = Mmseqs2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("sensitivity"), "msg: {msg}");
        assert!(msg.contains("7.5"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
