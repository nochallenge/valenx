//! `[bio.minimap2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "minimap2.align"
//!
//! [bio.minimap2]
//! reference  = "ref.fa"
//! reads      = ["reads.fq.gz"]   # one or more FASTQ / FASTA inputs
//! preset     = "map-ont"          # optional, defaults to "map-ont"
//! threads    = 4                  # optional, defaults to 1
//! extra_args = ["--secondary=no"] # optional, defaults to []
//! ```
//!
//! Unlike BWA, minimap2 builds its index on the fly — there is no
//! separate index step. The preset selects an internal scoring /
//! gap-extension profile tuned for the relevant read type or
//! genome-vs-genome workload:
//!
//! - `map-ont` — Oxford Nanopore long reads (default; the most
//!   common minimap2 invocation in modern pipelines)
//! - `map-pb`  — PacBio CLR long reads
//! - `map-hifi` — PacBio HiFi (very accurate long reads)
//! - `sr` — short-read alignment (Illumina-style)
//! - `asm5` / `asm10` / `asm20` — assembly-to-reference at varying
//!   divergence levels (5%, 10%, 20%)
//! - `splice` — RNA-seq spliced alignment
//! - `ava-pb` / `ava-ont` — all-vs-all overlap for assembly
//!
//! The preset list is canonicalised here so a typo in `case.toml`
//! fails fast in `validate()` rather than burning compute on a
//! minimap2 invocation that picks a wrong profile.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical preset list. Kept module-public so the adapter's
/// `validate()` (and any UI surface) can advertise the supported
/// values without redefining the slice.
pub const SUPPORTED_PRESETS: &[&str] = &[
    "map-ont", "map-pb", "map-hifi", "sr", "asm5", "asm10", "asm20", "splice", "ava-pb", "ava-ont",
];

#[derive(Clone, Debug, PartialEq)]
pub struct Minimap2Input {
    pub reference: PathBuf,
    pub reads: Vec<PathBuf>,
    pub preset: String,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl Minimap2Input {
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
            .and_then(|v| v.get("minimap2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.minimap2] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.minimap2].reference required"))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.minimap2].reference must not be empty"
            )));
        }

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.minimap2].reads required (array of FASTQ / FASTA paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.minimap2].reads entries must be strings"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.minimap2].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.minimap2].reads must contain at least one input file"
            )));
        }

        let preset = match block.get("preset") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.minimap2].preset must be a string"))
                })?;
                if !SUPPORTED_PRESETS.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.minimap2].preset `{s}` not recognised — \
                         expected one of {SUPPORTED_PRESETS:?}"
                    )));
                }
                s.to_string()
            }
            None => "map-ont".to_string(),
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.minimap2].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.minimap2].threads must be >= 1, got {raw}"
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
                        "[bio.minimap2].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.minimap2].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            reference: PathBuf::from(reference_str),
            reads,
            preset,
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
        // Minimal config: only the required fields. Everything else
        // should fall back to its documented default.
        let d = tempdir("minimap2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "minimap2.align"

[bio.minimap2]
reference = "ref.fa"
reads     = ["reads.fq.gz"]
"#,
        )
        .unwrap();
        let input = Minimap2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.reads, vec![PathBuf::from("reads.fq.gz")]);
        // Defaults: long-read ONT preset, 1 thread, no extra args.
        assert_eq!(input.preset, "map-ont");
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_preset() {
        // Typo in the preset name — a misspelt `"ont-map"` should
        // fail fast rather than be passed through to minimap2 (which
        // would accept it silently and pick a default profile).
        let d = tempdir("minimap2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "minimap2.align"

[bio.minimap2]
reference = "ref.fa"
reads     = ["reads.fq.gz"]
preset    = "ont-map"
"#,
        )
        .unwrap();
        let err = Minimap2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("map-ont"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_assembly_preset_and_extras() {
        // Genome-to-genome at 10% divergence with explicit threading
        // and a couple of `--secondary=no` style extras flowing through
        // verbatim.
        let d = tempdir("minimap2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "minimap2.align"

[bio.minimap2]
reference  = "ref.fa"
reads      = ["query_a.fa", "query_b.fa"]
preset     = "asm10"
threads    = 16
extra_args = ["--secondary=no", "-c"]
"#,
        )
        .unwrap();
        let input = Minimap2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.reads[0], PathBuf::from("query_a.fa"));
        assert_eq!(input.reads[1], PathBuf::from("query_b.fa"));
        assert_eq!(input.preset, "asm10");
        assert_eq!(input.threads, 16);
        assert_eq!(
            input.extra_args,
            vec!["--secondary=no".to_string(), "-c".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
