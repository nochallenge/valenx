//! `[bio.blast]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "blast.search"
//!
//! [bio.blast]
//! program  = "blastn"
//! query    = "query.fa"
//! database = "db/nt"      # prefix; database files live at db/nt.nhr, db/nt.nin, db/nt.nsq
//! evalue   = 1e-5         # optional, defaults to 10.0 (BLAST default)
//! outfmt   = 6            # optional, defaults to 0 (pairwise)
//! threads  = 4            # optional, defaults to 1
//! extra_args = ["-task", "blastn-short"]   # optional, defaults to []
//! ```
//!
//! `program` selects which BLAST+ binary the adapter wraps. The five
//! supported programs span the full nucleotide/protein cross-product:
//!
//! - `blastn`  — nucleotide query  vs nucleotide database
//! - `blastp`  — protein    query  vs protein    database
//! - `blastx`  — nucleotide query  (translated) vs protein    database
//! - `tblastn` — protein    query  vs nucleotide database (translated)
//! - `tblastx` — nucleotide query  (translated) vs nucleotide database (translated)
//!
//! `database` is the **path stem** of a BLAST database (built by
//! `makeblastdb`). The actual files on disk are sets of three with
//! suffixes `.nhr/.nin/.nsq` (nucleotide DBs) or `.phr/.pin/.psq`
//! (protein DBs); we validate that the directory containing the prefix
//! exists but don't require any specific extension to be present, since
//! the suffix layout depends on which program is being run.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical BLAST+ program list. Module-public so the UI can surface
/// the supported values without redefining them here.
pub const SUPPORTED_PROGRAMS: &[&str] = &["blastn", "blastp", "blastx", "tblastn", "tblastx"];

#[derive(Clone, Debug, PartialEq)]
pub struct BlastInput {
    pub program: String,
    pub query: PathBuf,
    pub database: PathBuf,
    pub evalue: f64,
    pub outfmt: u8,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl BlastInput {
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
            .and_then(|v| v.get("blast"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.blast] section",
                    case_toml.display()
                ))
            })?;

        let program = block
            .get("program")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.blast].program required (one of {SUPPORTED_PROGRAMS:?})"
                ))
            })?;
        if !SUPPORTED_PROGRAMS.contains(&program) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.blast].program `{program}` not recognised — \
                 expected one of {SUPPORTED_PROGRAMS:?}"
            )));
        }

        let query_str = block.get("query").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.blast].query required (path to FASTA query file)"
            ))
        })?;
        if query_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.blast].query must not be empty"
            )));
        }

        let database_str = block
            .get("database")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.blast].database required (BLAST database prefix path)"
                ))
            })?;
        if database_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.blast].database must not be empty"
            )));
        }

        // E-values are floats; admit integers too (e.g. `evalue = 10`)
        // and validate finite + positive. BLAST's own default is 10.0.
        let evalue = match block.get("evalue") {
            Some(v) => {
                let raw = if let Some(f) = v.as_float() {
                    f
                } else if let Some(i) = v.as_integer() {
                    i as f64
                } else {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.blast].evalue must be a number"
                    )));
                };
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.blast].evalue must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.blast].evalue must be positive, got {raw}"
                    )));
                }
                raw
            }
            None => 10.0,
        };

        // outfmt is 0..=18 in BLAST+ 2.10+ (0=pairwise, 6=tab, 7=tab-with-comments,
        // 10=csv, 11=ASN.1 binary, 17=SAM, 18=VCF, etc.). We don't pin a hard upper
        // bound — let BLAST itself reject anything it doesn't know — but reject
        // values that won't fit in u8 (which would be nonsense regardless).
        let outfmt = match block.get("outfmt") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.blast].outfmt must be an integer"))
                })?;
                if !(0..=255).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.blast].outfmt must be 0..=255, got {raw}"
                    )));
                }
                raw as u8
            }
            None => 0,
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.blast].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.blast].threads must be >= 1, got {raw}"
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
                        "[bio.blast].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.blast].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            program: program.to_string(),
            query: PathBuf::from(query_str),
            database: PathBuf::from(database_str),
            evalue,
            outfmt,
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
        // Only the three required keys: program, query, database.
        // Defaults: evalue 10.0, outfmt 0, threads 1, no extras.
        let d = tempdir("blast");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "blast.search"

[bio.blast]
program  = "blastn"
query    = "query.fa"
database = "db/nt"
"#,
        )
        .unwrap();
        let input = BlastInput::from_case_dir(&d).unwrap();
        assert_eq!(input.program, "blastn");
        assert_eq!(input.query, PathBuf::from("query.fa"));
        assert_eq!(input.database, PathBuf::from("db/nt"));
        assert_eq!(input.evalue, 10.0);
        assert_eq!(input.outfmt, 0);
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_overrides() {
        // Protein-vs-protein with explicit threading, a strict E-value
        // cutoff, tab-separated output (outfmt 6), and a typical
        // task-tuning extra arg.
        let d = tempdir("blast");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "blast.search"

[bio.blast]
program    = "blastp"
query      = "query.fa"
database   = "db/nr"
evalue     = 1e-5
outfmt     = 6
threads    = 8
extra_args = ["-task", "blastp-short"]
"#,
        )
        .unwrap();
        let input = BlastInput::from_case_dir(&d).unwrap();
        assert_eq!(input.program, "blastp");
        assert_eq!(input.evalue, 1e-5);
        assert_eq!(input.outfmt, 6);
        assert_eq!(input.threads, 8);
        assert_eq!(
            input.extra_args,
            vec!["-task".to_string(), "blastp-short".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("blast");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = BlastInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.blast]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
