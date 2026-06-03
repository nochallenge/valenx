//! `[bio.foldseek]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "foldseek.search"
//!
//! [bio.foldseek]
//! query           = "query.pdb"
//! database        = "db/pdb100"   # prefix; FoldSeek DB files live at db/pdb100_*
//! output_basename = "search_results"
//! threads         = 4              # optional, defaults to 1
//! extra_args      = []             # optional, defaults to []
//! ```
//!
//! `query` is the path to a structure file (PDB or mmCIF) — FoldSeek
//! accepts either. `database` is the **path stem** of a FoldSeek
//! structure database (built by `foldseek createdb` or downloaded with
//! `foldseek databases`). The on-disk files are sets named
//! `<prefix>`, `<prefix>.dbtype`, `<prefix>.index`, `<prefix>_ss`,
//! `<prefix>_ca`, etc.; we validate the parent directory exists but
//! don't enforce a specific extension to be present (same approach as
//! BLAST).
//!
//! `output_basename` is the stem of the result file FoldSeek writes —
//! the actual file lives at `<workdir>/<output_basename>.m8` (BLAST
//! tab-delimited format). Pinning the basename here keeps `prepare()`
//! and `collect()` agreed on what to look for.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct FoldseekInput {
    pub query: PathBuf,
    pub database: PathBuf,
    pub output_basename: String,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl FoldseekInput {
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
            .and_then(|v| v.get("foldseek"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.foldseek] section",
                    case_toml.display()
                ))
            })?;

        let query_str = block.get("query").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.foldseek].query required (path to structure query, PDB or mmCIF)"
            ))
        })?;
        if query_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.foldseek].query must not be empty"
            )));
        }

        let database_str = block
            .get("database")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.foldseek].database required (FoldSeek database prefix path)"
                ))
            })?;
        if database_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.foldseek].database must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.foldseek].output_basename required (results file stem)"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.foldseek].output_basename must not be empty"
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.foldseek].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.foldseek].threads must be >= 1, got {raw}"
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
                        "[bio.foldseek].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.foldseek].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            query: PathBuf::from(query_str),
            database: PathBuf::from(database_str),
            output_basename: output_basename.to_string(),
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
    fn parses_minimal_case() {
        // Only the three required keys: query, database, output_basename.
        // Defaults: threads 1, no extras.
        let d = tempdir("foldseek");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "foldseek.search"

[bio.foldseek]
query           = "query.pdb"
database        = "db/pdb100"
output_basename = "search_results"
"#,
        )
        .unwrap();
        let input = FoldseekInput::from_case_dir(&d).unwrap();
        assert_eq!(input.query, PathBuf::from("query.pdb"));
        assert_eq!(input.database, PathBuf::from("db/pdb100"));
        assert_eq!(input.output_basename, "search_results");
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_overrides() {
        // Threaded run with extra args (e.g. an alignment-mode tweak
        // and an E-value cutoff — both real FoldSeek flags).
        let d = tempdir("foldseek");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "foldseek.search"

[bio.foldseek]
query           = "query.pdb"
database        = "db/pdb100"
output_basename = "search_results"
threads         = 8
extra_args      = ["--alignment-type", "2", "-e", "0.001"]
"#,
        )
        .unwrap();
        let input = FoldseekInput::from_case_dir(&d).unwrap();
        assert_eq!(input.threads, 8);
        assert_eq!(
            input.extra_args,
            vec![
                "--alignment-type".to_string(),
                "2".to_string(),
                "-e".to_string(),
                "0.001".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("foldseek");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = FoldseekInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.foldseek]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
