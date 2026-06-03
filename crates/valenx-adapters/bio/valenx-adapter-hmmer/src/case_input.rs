//! `[bio.hmmer]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "hmmer.search"
//!
//! [bio.hmmer]
//! tool       = "hmmsearch"           # one of "hmmsearch", "hmmscan"
//! profile    = "Pfam-A.hmm"          # the .hmm file (or .hmm-database for hmmscan)
//! sequences  = "queries.fa"          # the FASTA (sequence DB for hmmsearch, query for hmmscan)
//! cpus       = 4                     # optional, defaults to 1
//! evalue     = 0.01                  # optional, defaults to 0.01
//! extra_args = ["--noali"]           # optional, defaults to []
//! ```
//!
//! `tool` selects which HMMER subcommand the adapter wraps:
//!
//! - `hmmsearch` — search a single profile against a sequence database.
//!   `<profile>` is the `.hmm` file, `<sequences>` is the FASTA db.
//! - `hmmscan` — search a query FASTA against a profile database.
//!   `<profile>` is the pressed `.hmm` profile-DB, `<sequences>` is the
//!   query FASTA.
//!
//! Both subcommands write a tabular summary (`--tblout`) and a verbose
//! human-readable report (`-o`); the adapter pins both to fixed names
//! in the workdir so `collect()` can reliably surface them.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical HMMER subcommand list. Module-public so the UI can
/// surface the supported values without redefining them here.
pub const SUPPORTED_TOOLS: &[&str] = &["hmmsearch", "hmmscan"];

#[derive(Clone, Debug, PartialEq)]
pub struct HmmerInput {
    pub tool: String,
    pub profile: PathBuf,
    pub sequences: PathBuf,
    pub cpus: u32,
    pub evalue: f64,
    pub extra_args: Vec<String>,
}

impl HmmerInput {
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
            .and_then(|v| v.get("hmmer"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.hmmer] section",
                    case_toml.display()
                ))
            })?;

        let tool = block.get("tool").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.hmmer].tool required (one of {SUPPORTED_TOOLS:?})"
            ))
        })?;
        if !SUPPORTED_TOOLS.contains(&tool) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hmmer].tool `{tool}` not recognised — \
                 expected one of {SUPPORTED_TOOLS:?}"
            )));
        }

        let profile_str = block
            .get("profile")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.hmmer].profile required (path to .hmm file or pressed profile DB)"
                ))
            })?;
        if profile_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hmmer].profile must not be empty"
            )));
        }

        let sequences_str = block
            .get("sequences")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.hmmer].sequences required (path to FASTA)"
                ))
            })?;
        if sequences_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hmmer].sequences must not be empty"
            )));
        }

        let cpus = match block.get("cpus") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.hmmer].cpus must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.hmmer].cpus must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        // E-values are floats; HMMER accepts integers like `1` too,
        // so admit both. Either way we validate finite + positive.
        let evalue = match block.get("evalue") {
            Some(v) => {
                let raw = if let Some(f) = v.as_float() {
                    f
                } else if let Some(i) = v.as_integer() {
                    i as f64
                } else {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.hmmer].evalue must be a number"
                    )));
                };
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.hmmer].evalue must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.hmmer].evalue must be positive, got {raw}"
                    )));
                }
                raw
            }
            None => 0.01,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.hmmer].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.hmmer].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            tool: tool.to_string(),
            profile: PathBuf::from(profile_str),
            sequences: PathBuf::from(sequences_str),
            cpus,
            evalue,
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
        // Only the three required keys: tool, profile, sequences.
        // Defaults: 1 cpu, evalue 0.01, no extras.
        let d = tempdir("hmmer");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hmmer.search"

[bio.hmmer]
tool      = "hmmsearch"
profile   = "Pfam-A.hmm"
sequences = "queries.fa"
"#,
        )
        .unwrap();
        let input = HmmerInput::from_case_dir(&d).unwrap();
        assert_eq!(input.tool, "hmmsearch");
        assert_eq!(input.profile, PathBuf::from("Pfam-A.hmm"));
        assert_eq!(input.sequences, PathBuf::from("queries.fa"));
        assert_eq!(input.cpus, 1);
        assert_eq!(input.evalue, 0.01);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_tool() {
        // `phmmer` exists in the HMMER suite but isn't one of the
        // wrapped subcommands — must be rejected up front.
        let d = tempdir("hmmer");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hmmer.search"

[bio.hmmer]
tool      = "phmmer"
profile   = "Pfam-A.hmm"
sequences = "queries.fa"
"#,
        )
        .unwrap();
        let err = HmmerInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("hmmsearch"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_negative_evalue() {
        // E-values are probabilities (capped at 1.0 inside HMMER) —
        // anything <= 0 is meaningless.
        let d = tempdir("hmmer");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hmmer.search"

[bio.hmmer]
tool      = "hmmsearch"
profile   = "Pfam-A.hmm"
sequences = "queries.fa"
evalue    = -0.5
"#,
        )
        .unwrap();
        let err = HmmerInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("positive"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("hmmer");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = HmmerInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.hmmer]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_hmmscan_with_overrides() {
        // hmmscan with explicit cpus + evalue + extras to exercise
        // the non-default branches in the parser.
        let d = tempdir("hmmer");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hmmer.search"

[bio.hmmer]
tool       = "hmmscan"
profile    = "Pfam-A.hmm"
sequences  = "query.fa"
cpus       = 8
evalue     = 1e-10
extra_args = ["--noali", "--cut_ga"]
"#,
        )
        .unwrap();
        let input = HmmerInput::from_case_dir(&d).unwrap();
        assert_eq!(input.tool, "hmmscan");
        assert_eq!(input.cpus, 8);
        assert_eq!(input.evalue, 1e-10);
        assert_eq!(
            input.extra_args,
            vec!["--noali".to_string(), "--cut_ga".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
