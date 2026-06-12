//! `[bio.star]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "star.align"
//!
//! [bio.star]
//! genome_dir  = "star_index"                            # required: STAR index dir, or output dir for index-gen
//! reference   = "ref.fa"                                # required when skip_index = false
//! reads       = ["reads_R1.fq.gz", "reads_R2.fq.gz"]    # 1 (single-end) or 2 (paired-end)
//! threads     = 4                                        # optional, defaults to 1
//! skip_index  = false                                    # optional, defaults to false
//! output_type = "BAM_SortedByCoordinate"                 # optional, default BAM_SortedByCoordinate
//! sjdb_gtf    = "annotations.gtf"                        # optional, splice-junction annotation for index-gen
//! extra_args  = ["--outSAMattributes", "Standard"]       # optional, defaults to []
//! ```
//!
//! STAR is unusual in that the genome index isn't a small set of
//! files next to the FASTA — it's a multi-GB **directory** the user
//! either pre-builds or generates inline via `--runMode
//! genomeGenerate`. We model this with two paths:
//!
//! - `genome_dir` is always required: it's both the input to
//!   `--runMode alignReads --genomeDir` and the output of
//!   `--runMode genomeGenerate --genomeDir` when we build it
//!   ourselves.
//! - `reference` is only required when `skip_index = false`. When
//!   the user has a pre-built index they don't need a FASTA at all.
//!
//! `output_type` selects the alignment output format. The three
//! supported values map to STAR's `--outSAMtype` flag:
//!
//! - `BAM_Unsorted` → `--outSAMtype BAM Unsorted` (fastest, useful
//!   when piping into downstream tools that re-sort anyway).
//! - `BAM_SortedByCoordinate` → `--outSAMtype BAM
//!   SortedByCoordinate` (default, ready for `samtools index`).
//! - `SAM` → `--outSAMtype SAM` (text output, debugging).
//!
//! `sjdb_gtf` is the splice-junction annotation file used at
//! index-generation time (`--sjdbGTFfile`). Providing it gives STAR
//! prior knowledge of every annotated junction; the matching
//! `--sjdbOverhang 100` is set automatically (matches the standard
//! ~100-bp Illumina read).

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical STAR output-type values. Module-public so the UI can
/// surface the supported values without redefining them here.
pub const SUPPORTED_OUTPUT_TYPES: &[&str] = &["BAM_Unsorted", "BAM_SortedByCoordinate", "SAM"];

#[derive(Clone, Debug, PartialEq)]
pub struct StarInput {
    pub genome_dir: PathBuf,
    pub reference: Option<PathBuf>,
    pub reads: Vec<PathBuf>,
    pub threads: u32,
    pub skip_index: bool,
    pub output_type: String,
    pub sjdb_gtf: Option<PathBuf>,
    pub extra_args: Vec<String>,
}

impl StarInput {
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
            .and_then(|v| v.get("star"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.star] section",
                    case_toml.display()
                ))
            })?;

        let genome_dir_str = block
            .get("genome_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.star].genome_dir required"))
            })?;
        if genome_dir_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.star].genome_dir must not be empty"
            )));
        }

        let reference: Option<PathBuf> = match block.get("reference") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.star].reference must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.star].reference must not be empty"
                    )));
                }
                Some(PathBuf::from(s))
            }
            None => None,
        };

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.star].reads required (array of 1 or 2 FASTQ paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.star].reads entries must be strings"))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.star].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() || reads.len() > 2 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.star].reads must contain 1 (single-end) or 2 \
                 (paired-end) FASTQs, got {}",
                reads.len()
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.star].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.star].threads must be >= 1, got {raw}"
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

        let output_type = match block.get("output_type") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.star].output_type must be a string"))
                })?;
                if !SUPPORTED_OUTPUT_TYPES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.star].output_type `{s}` not recognised — \
                         expected one of {SUPPORTED_OUTPUT_TYPES:?}"
                    )));
                }
                s.to_string()
            }
            None => "BAM_SortedByCoordinate".to_string(),
        };

        let sjdb_gtf: Option<PathBuf> = match block.get("sjdb_gtf") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.star].sjdb_gtf must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.star].sjdb_gtf must not be empty"
                    )));
                }
                Some(PathBuf::from(s))
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.star].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.star].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        // When we're going to build the index ourselves, we need a
        // FASTA to build it *from*. The user can leave `reference`
        // off only if they're skipping index generation (`skip_index
        // = true`).
        if !skip_index && reference.is_none() {
            return Err(AdapterError::InvalidCase {
                case_path: case_toml,
                reason: "[bio.star].reference is required when skip_index = false \
                         (need a FASTA to build the genome index from)"
                    .into(),
            });
        }

        Ok(Self {
            genome_dir: PathBuf::from(genome_dir_str),
            reference,
            reads,
            threads,
            skip_index,
            output_type,
            sjdb_gtf,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_with_skip_index() {
        // Single-end alignment against a pre-built index: no
        // reference needed, default threads, default
        // BAM_SortedByCoordinate output.
        let d = tempdir("star");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "star.align"

[bio.star]
genome_dir = "star_index"
reads      = ["reads.fq.gz"]
skip_index = true
"#,
        )
        .unwrap();
        let input = StarInput::from_case_dir(&d).unwrap();
        assert_eq!(input.genome_dir, PathBuf::from("star_index"));
        assert!(input.reference.is_none());
        assert_eq!(input.reads, vec![PathBuf::from("reads.fq.gz")]);
        assert_eq!(input.threads, 1);
        assert!(input.skip_index);
        assert_eq!(input.output_type, "BAM_SortedByCoordinate");
        assert!(input.sjdb_gtf.is_none());
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_paired_reads_and_sjdb() {
        // Paired-end, building the index inline with a GTF, custom
        // BAM_Unsorted output, threading on.
        let d = tempdir("star");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "star.align"

[bio.star]
genome_dir  = "star_index"
reference   = "ref.fa"
reads       = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
threads     = 8
output_type = "BAM_Unsorted"
sjdb_gtf    = "annotations.gtf"
extra_args  = ["--outSAMattributes", "Standard"]
"#,
        )
        .unwrap();
        let input = StarInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.reads[1], PathBuf::from("reads_R2.fq.gz"));
        assert_eq!(input.threads, 8);
        assert!(!input.skip_index);
        assert_eq!(input.reference, Some(PathBuf::from("ref.fa")));
        assert_eq!(input.output_type, "BAM_Unsorted");
        assert_eq!(input.sjdb_gtf, Some(PathBuf::from("annotations.gtf")));
        assert_eq!(
            input.extra_args,
            vec!["--outSAMattributes".to_string(), "Standard".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_three_reads() {
        // STAR is single-end (1 FASTQ) or paired-end (2 FASTQs).
        // Three is never a valid invocation.
        let d = tempdir("star");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "star.align"

[bio.star]
genome_dir = "star_index"
reads      = ["a.fq", "b.fq", "c.fq"]
skip_index = true
"#,
        )
        .unwrap();
        let err = StarInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("1 (single-end) or 2 (paired-end)"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_output_type() {
        // `CRAM` is a samtools format, not a STAR output mode — must
        // be rejected up front so the user sees a clean error rather
        // than waiting for STAR to bail at runtime.
        let d = tempdir("star");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "star.align"

[bio.star]
genome_dir  = "star_index"
reads       = ["reads.fq"]
skip_index  = true
output_type = "CRAM"
"#,
        )
        .unwrap();
        let err = StarInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("BAM_SortedByCoordinate"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_index_generation_without_reference() {
        // Asking us to build the index (`skip_index = false`) without
        // a FASTA is incoherent — STAR's `genomeGenerate` mode needs
        // `--genomeFastaFiles` or it has nothing to index.
        let d = tempdir("star");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "star.align"

[bio.star]
genome_dir = "star_index"
reads      = ["reads.fq"]
"#,
        )
        .unwrap();
        let err = StarInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("reference is required when skip_index = false"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
