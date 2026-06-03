//! `[bio.diamond]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "diamond.blastp"
//!
//! [bio.diamond]
//! action      = "blastp"               # one of "blastp" | "blastx" | "makedb"
//! query       = "query.fa"             # query FASTA (blastp/blastx) OR input FASTA (makedb)
//! database    = "uniref90.dmnd"        # DIAMOND DB (blastp/blastx) OR output basename (makedb)
//! output      = "hits.m8"              # output hit table (blastp/blastx); set anything for makedb
//! sensitivity = "default"              # default | fast | sensitive | more-sensitive | very-sensitive | ultra-sensitive
//! threads     = 8                      # optional, defaults to 1
//! extra_args  = ["-e", "1e-5"]         # optional, defaults to []
//! ```
//!
//! `action` selects which DIAMOND mode the adapter wraps:
//!
//! - `blastp` — protein-vs-protein search.
//! - `blastx` — translated nucleotide-vs-protein search.
//! - `makedb` — build a `.dmnd` database from a FASTA. In this mode
//!   `query` is the input FASTA and `database` is the **output** DB
//!   basename (DIAMOND appends `.dmnd`); that's how DIAMOND's CLI is
//!   actually shaped — same field name, different role per action.
//!
//! Sensitivity values map directly to DIAMOND's `--<sensitivity>`
//! flags. `default` is special: DIAMOND has no `--default` flag —
//! omitting any sensitivity flag *is* the default — so the adapter
//! drops the flag entirely when this preset is selected.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical DIAMOND action list. Module-public so the UI can surface
/// the supported values without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["blastp", "blastx", "makedb"];

/// Canonical DIAMOND sensitivity preset list. Ordered roughly from
/// fastest to slowest / most sensitive. `default` represents "no
/// flag" — DIAMOND's built-in default — and is filtered out of the
/// command line in `lib.rs`.
pub const SUPPORTED_SENSITIVITIES: &[&str] = &[
    "default",
    "fast",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
];

#[derive(Clone, Debug, PartialEq)]
pub struct DiamondInput {
    pub action: String,
    pub query: PathBuf,
    pub database: PathBuf,
    pub output: PathBuf,
    pub sensitivity: String,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl DiamondInput {
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
            .and_then(|v| v.get("diamond"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.diamond] section",
                    case_toml.display()
                ))
            })?;

        let action = block
            .get("action")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.diamond].action required (one of {SUPPORTED_ACTIONS:?})"
                ))
            })?;
        if !SUPPORTED_ACTIONS.contains(&action) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.diamond].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        let query_str = block
            .get("query")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.diamond].query required")))?;
        if query_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.diamond].query must not be empty"
            )));
        }

        let database_str = block
            .get("database")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.diamond].database required"))
            })?;
        if database_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.diamond].database must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.diamond].output required")))?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.diamond].output must not be empty"
            )));
        }

        let sensitivity = match block.get("sensitivity") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.diamond].sensitivity must be a string"
                    ))
                })?;
                if !SUPPORTED_SENSITIVITIES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.diamond].sensitivity `{s}` not recognised — \
                         expected one of {SUPPORTED_SENSITIVITIES:?}"
                    )));
                }
                s.to_string()
            }
            None => "default".to_string(),
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.diamond].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.diamond].threads must be >= 1, got {raw}"
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
                        "[bio.diamond].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.diamond].extra_args entries must be strings"
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
            database: PathBuf::from(database_str),
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
    fn parses_blastp_minimal() {
        // Defaults: sensitivity = "default" (no flag), threads = 1, no
        // extras.
        let d = tempdir("diamond");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "diamond.blastp"

[bio.diamond]
action   = "blastp"
query    = "query.fa"
database = "uniref90.dmnd"
output   = "hits.m8"
"#,
        )
        .unwrap();
        let input = DiamondInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "blastp");
        assert_eq!(input.query, PathBuf::from("query.fa"));
        assert_eq!(input.database, PathBuf::from("uniref90.dmnd"));
        assert_eq!(input.output, PathBuf::from("hits.m8"));
        assert_eq!(input.sensitivity, "default");
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_makedb() {
        // makedb: query is the input FASTA, database is the *output*
        // DB basename. The schema field name is the same, the role
        // flips per action — that's DIAMOND's actual CLI.
        let d = tempdir("diamond");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "diamond.makedb"

[bio.diamond]
action   = "makedb"
query    = "uniref90.fasta"
database = "uniref90"
output   = "ignored.m8"
threads  = 16
"#,
        )
        .unwrap();
        let input = DiamondInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "makedb");
        assert_eq!(input.query, PathBuf::from("uniref90.fasta"));
        assert_eq!(input.database, PathBuf::from("uniref90"));
        assert_eq!(input.threads, 16);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_sensitivity_override() {
        let d = tempdir("diamond");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "diamond.blastx"

[bio.diamond]
action      = "blastx"
query       = "reads.fa"
database    = "nr.dmnd"
output      = "hits.m8"
sensitivity = "ultra-sensitive"
threads     = 32
extra_args  = ["-e", "1e-10"]
"#,
        )
        .unwrap();
        let input = DiamondInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "blastx");
        assert_eq!(input.sensitivity, "ultra-sensitive");
        assert_eq!(input.threads, 32);
        assert_eq!(
            input.extra_args,
            vec!["-e".to_string(), "1e-10".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `diamond cluster` is real DIAMOND but not one this adapter
        // wraps. Reject up front.
        let d = tempdir("diamond");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "diamond.cluster"

[bio.diamond]
action   = "cluster"
query    = "in.fa"
database = "db"
output   = "out"
"#,
        )
        .unwrap();
        let err = DiamondInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("blastp"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_sensitivity() {
        // `--super-sensitive` is a Bowtie2 preset; DIAMOND has its own
        // ladder. Surface the supported list in the error message.
        let d = tempdir("diamond");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "diamond.blastp"

[bio.diamond]
action      = "blastp"
query       = "query.fa"
database    = "db.dmnd"
output      = "hits.m8"
sensitivity = "super-sensitive"
"#,
        )
        .unwrap();
        let err = DiamondInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("ultra-sensitive"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
