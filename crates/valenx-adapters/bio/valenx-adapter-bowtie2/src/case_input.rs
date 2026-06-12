//! `[bio.bowtie2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "bowtie2.align"
//!
//! [bio.bowtie2]
//! reference  = "ref.fa"
//! reads      = ["reads_R1.fq.gz", "reads_R2.fq.gz"]   # 1 (single-end) or 2 (paired-end)
//! threads    = 4                                       # optional, defaults to 1
//! skip_index = false                                   # optional, defaults to false
//! preset     = "sensitive"                             # optional, default "sensitive"
//! extra_args = ["--no-unal"]                           # optional, defaults to []
//! ```
//!
//! Bowtie2 takes either a single-end FASTQ or a paired-end pair of
//! FASTQs; it does **not** support 3+ read inputs as a single
//! invocation. `skip_index = true` reuses an existing Bowtie2 index
//! sitting next to the reference (`ref.fa.1.bt2`, `.2.bt2`, `.3.bt2`,
//! `.4.bt2`, `.rev.1.bt2`, `.rev.2.bt2`) — saves the bowtie2-build cost
//! on successive runs over the same reference.
//!
//! `preset` selects one of the four canonical end-to-end alignment
//! presets: `very-fast`, `fast`, `sensitive` (default), or
//! `very-sensitive`. These map directly to Bowtie2's `--<preset>`
//! flags.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical Bowtie2 alignment-preset list. Module-public so the UI
/// can surface the supported values without redefining them here.
pub const SUPPORTED_PRESETS: &[&str] = &["very-fast", "fast", "sensitive", "very-sensitive"];

#[derive(Clone, Debug, PartialEq)]
pub struct Bowtie2Input {
    pub reference: PathBuf,
    pub reads: Vec<PathBuf>,
    pub threads: u32,
    pub skip_index: bool,
    pub preset: String,
    pub extra_args: Vec<String>,
}

impl Bowtie2Input {
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
            .and_then(|v| v.get("bowtie2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.bowtie2] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.bowtie2].reference required"))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bowtie2].reference must not be empty"
            )));
        }

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.bowtie2].reads required (array of 1 or 2 FASTQ paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.bowtie2].reads entries must be strings"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.bowtie2].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() || reads.len() > 2 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bowtie2].reads must contain 1 (single-end) or 2 \
                 (paired-end) FASTQs, got {}",
                reads.len()
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.bowtie2].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bowtie2].threads must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let skip_index = block
            .get("skip_index")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let preset = match block.get("preset") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.bowtie2].preset must be a string"))
                })?;
                if !SUPPORTED_PRESETS.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bowtie2].preset `{s}` not recognised — \
                         expected one of {SUPPORTED_PRESETS:?}"
                    )));
                }
                s.to_string()
            }
            None => "sensitive".to_string(),
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.bowtie2].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.bowtie2].extra_args entries must be strings"
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
            threads,
            skip_index,
            preset,
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
        // Single-end alignment with all defaults: 1 thread, build the
        // index, "sensitive" preset, no extras.
        let d = tempdir("bowtie2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bowtie2.align"

[bio.bowtie2]
reference = "ref.fa"
reads     = ["reads.fq.gz"]
"#,
        )
        .unwrap();
        let input = Bowtie2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.reads, vec![PathBuf::from("reads.fq.gz")]);
        assert_eq!(input.threads, 1);
        assert!(!input.skip_index);
        assert_eq!(input.preset, "sensitive");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_paired_reads_and_preset() {
        // Paired-end, custom preset, multi-thread, pre-built index,
        // typical extra arg suppressing unaligned reads from the SAM.
        let d = tempdir("bowtie2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bowtie2.align"

[bio.bowtie2]
reference  = "ref.fa"
reads      = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
threads    = 8
skip_index = true
preset     = "very-sensitive"
extra_args = ["--no-unal"]
"#,
        )
        .unwrap();
        let input = Bowtie2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.reads[1], PathBuf::from("reads_R2.fq.gz"));
        assert_eq!(input.threads, 8);
        assert!(input.skip_index);
        assert_eq!(input.preset, "very-sensitive");
        assert_eq!(input.extra_args, vec!["--no-unal".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_three_reads() {
        // Bowtie2 is single-end (1 FASTQ) or paired-end (2 FASTQs).
        // Three is never a valid invocation.
        let d = tempdir("bowtie2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bowtie2.align"

[bio.bowtie2]
reference = "ref.fa"
reads     = ["a.fq", "b.fq", "c.fq"]
"#,
        )
        .unwrap();
        let err = Bowtie2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("1 (single-end) or 2 (paired-end)"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        let d = tempdir("bowtie2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bowtie2.align"

[bio.bowtie2]
reference = "ref.fa"
reads     = ["reads.fq"]
threads   = 0
"#,
        )
        .unwrap();
        let err = Bowtie2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads must be >= 1"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_preset() {
        // `--ultra-sensitive` is a DIAMOND preset, not Bowtie2 — must be
        // rejected up front so the user sees a helpful error rather
        // than waiting for bowtie2 to crash.
        let d = tempdir("bowtie2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bowtie2.align"

[bio.bowtie2]
reference = "ref.fa"
reads     = ["reads.fq"]
preset    = "ultra-sensitive"
"#,
        )
        .unwrap();
        let err = Bowtie2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("very-sensitive"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
