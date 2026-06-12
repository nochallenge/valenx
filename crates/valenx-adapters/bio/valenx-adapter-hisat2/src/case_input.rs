//! `[bio.hisat2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "hisat2.align"
//!
//! [bio.hisat2]
//! reference   = "ref.fa"
//! reads       = ["reads_R1.fq.gz", "reads_R2.fq.gz"]   # 1 (single-end) or 2 (paired-end)
//! threads     = 4                                       # optional, defaults to 1
//! skip_index  = false                                   # optional, defaults to false
//! strandness  = "unstranded"                            # optional, default "unstranded"
//! extra_args  = ["--no-unal"]                           # optional, defaults to []
//! ```
//!
//! HISAT2 accepts a single-end FASTQ or a paired-end pair of FASTQs.
//! `skip_index = true` reuses an existing HISAT2 graph index sitting
//! next to the reference (`<base>.1.ht2`, `.2.ht2`, … `.8.ht2`) so
//! successive runs over the same reference don't pay the
//! hisat2-build cost.
//!
//! `strandness` selects the RNA-seq strand-specificity protocol:
//! `unstranded` (default — no strand info), `F` / `R` (single-end:
//! sense / antisense to the transcript), or `FR` / `RF` (paired-end:
//! the standard dUTP / Illumina TruSeq strand orientations). These
//! map to HISAT2's `--rna-strandness` flag values.
//!
//! Strand-aware modes matter for downstream quantification — feeding
//! the wrong protocol to `featureCounts` flips antisense reads onto
//! the wrong gene. We surface the choice up front so the case file
//! captures it in provenance.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical HISAT2 strandness values. Module-public so the UI can
/// surface the supported values without redefining them here.
pub const SUPPORTED_STRANDNESS: &[&str] = &["unstranded", "F", "R", "FR", "RF"];

#[derive(Clone, Debug, PartialEq)]
pub struct Hisat2Input {
    pub reference: PathBuf,
    pub reads: Vec<PathBuf>,
    pub threads: u32,
    pub skip_index: bool,
    pub strandness: String,
    pub extra_args: Vec<String>,
}

impl Hisat2Input {
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
            .and_then(|v| v.get("hisat2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.hisat2] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.hisat2].reference required"))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hisat2].reference must not be empty"
            )));
        }

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.hisat2].reads required (array of 1 or 2 FASTQ paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.hisat2].reads entries must be strings"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.hisat2].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() || reads.len() > 2 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.hisat2].reads must contain 1 (single-end) or 2 \
                 (paired-end) FASTQs, got {}",
                reads.len()
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.hisat2].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.hisat2].threads must be >= 1, got {raw}"
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

        let strandness = match block.get("strandness") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.hisat2].strandness must be a string"))
                })?;
                if !SUPPORTED_STRANDNESS.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.hisat2].strandness `{s}` not recognised — \
                         expected one of {SUPPORTED_STRANDNESS:?}"
                    )));
                }
                s.to_string()
            }
            None => "unstranded".to_string(),
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.hisat2].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.hisat2].extra_args entries must be strings"
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
            strandness,
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
        // index, "unstranded" RNA-seq, no extras.
        let d = tempdir("hisat2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hisat2.align"

[bio.hisat2]
reference = "ref.fa"
reads     = ["reads.fq.gz"]
"#,
        )
        .unwrap();
        let input = Hisat2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.reads, vec![PathBuf::from("reads.fq.gz")]);
        assert_eq!(input.threads, 1);
        assert!(!input.skip_index);
        assert_eq!(input.strandness, "unstranded");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_paired_reads_and_strandness() {
        // Paired-end, FR strandness (Illumina TruSeq stranded), 8
        // threads, pre-built index, suppress unaligned reads from the
        // SAM.
        let d = tempdir("hisat2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hisat2.align"

[bio.hisat2]
reference  = "ref.fa"
reads      = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
threads    = 8
skip_index = true
strandness = "FR"
extra_args = ["--no-unal"]
"#,
        )
        .unwrap();
        let input = Hisat2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.reads[1], PathBuf::from("reads_R2.fq.gz"));
        assert_eq!(input.threads, 8);
        assert!(input.skip_index);
        assert_eq!(input.strandness, "FR");
        assert_eq!(input.extra_args, vec!["--no-unal".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_three_reads() {
        // HISAT2 is single-end (1 FASTQ) or paired-end (2 FASTQs).
        // Three is never a valid invocation.
        let d = tempdir("hisat2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hisat2.align"

[bio.hisat2]
reference = "ref.fa"
reads     = ["a.fq", "b.fq", "c.fq"]
"#,
        )
        .unwrap();
        let err = Hisat2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("1 (single-end) or 2 (paired-end)"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        let d = tempdir("hisat2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hisat2.align"

[bio.hisat2]
reference = "ref.fa"
reads     = ["reads.fq"]
threads   = 0
"#,
        )
        .unwrap();
        let err = Hisat2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads must be >= 1"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_strandness() {
        // `--fr-stranded` would be a salmon flag, not a valid HISAT2
        // strandness mode — must be rejected up front so the user
        // sees a clean error rather than a runtime crash.
        let d = tempdir("hisat2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "hisat2.align"

[bio.hisat2]
reference  = "ref.fa"
reads      = ["reads.fq"]
strandness = "fr-stranded"
"#,
        )
        .unwrap();
        let err = Hisat2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("unstranded"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
