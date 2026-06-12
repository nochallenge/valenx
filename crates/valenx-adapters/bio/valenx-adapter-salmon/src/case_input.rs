//! `[bio.salmon]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "salmon.quant"
//!
//! [bio.salmon]
//! transcriptome = "transcripts.fa"
//! index_dir     = "salmon_index"                       # input/output dir for the salmon index
//! reads         = ["reads_R1.fq.gz", "reads_R2.fq.gz"] # 1 (single-end) or 2 (paired-end)
//! output_dir    = "salmon_quant"                       # quant output dir
//! threads       = 4                                     # optional, defaults to 1
//! skip_index    = false                                 # optional, defaults to false
//! libtype       = "A"                                   # optional, defaults to "A" (auto-detect)
//! extra_args    = ["--validateMappings"]                # optional, defaults to []
//! ```
//!
//! Salmon accepts a single-end FASTQ or a paired-end pair of FASTQs.
//! `skip_index = true` reuses an existing salmon index sitting at
//! `index_dir` so successive runs over the same transcriptome don't
//! pay the `salmon index` cost.
//!
//! `libtype` selects the library type Salmon uses to decide how reads
//! relate to transcripts. The default `"A"` asks Salmon to
//! auto-detect; common explicit values include `"IU"` (paired,
//! unstranded), `"ISF"` / `"ISR"` (paired, stranded — TruSeq /
//! dUTP), `"U"` / `"SF"` / `"SR"` (single-end variants). We don't
//! whitelist the supported values — Salmon prints a clear error if
//! the libtype is bogus, and the matrix of valid combinations is
//! large enough that re-deriving it here would just go stale.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SalmonInput {
    pub transcriptome: PathBuf,
    pub index_dir: PathBuf,
    pub reads: Vec<PathBuf>,
    pub output_dir: PathBuf,
    pub threads: u32,
    pub skip_index: bool,
    pub libtype: String,
    pub extra_args: Vec<String>,
}

impl SalmonInput {
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
            .and_then(|v| v.get("salmon"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.salmon] section",
                    case_toml.display()
                ))
            })?;

        let transcriptome_str = block
            .get("transcriptome")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.salmon].transcriptome required"))
            })?;
        if transcriptome_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.salmon].transcriptome must not be empty"
            )));
        }

        let index_dir_str = block
            .get("index_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.salmon].index_dir required"))
            })?;
        if index_dir_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.salmon].index_dir must not be empty"
            )));
        }

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.salmon].reads required (array of 1 or 2 FASTQ paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.salmon].reads entries must be strings"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.salmon].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() || reads.len() > 2 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.salmon].reads must contain 1 (single-end) or 2 \
                 (paired-end) FASTQs, got {}",
                reads.len()
            )));
        }

        let output_dir_str = block
            .get("output_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.salmon].output_dir required"))
            })?;
        if output_dir_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.salmon].output_dir must not be empty"
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.salmon].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.salmon].threads must be >= 1, got {raw}"
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

        let libtype = match block.get("libtype") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.salmon].libtype must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.salmon].libtype must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "A".to_string(),
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.salmon].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.salmon].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            transcriptome: PathBuf::from(transcriptome_str),
            index_dir: PathBuf::from(index_dir_str),
            reads,
            output_dir: PathBuf::from(output_dir_str),
            threads,
            skip_index,
            libtype,
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
        // Single-end quantification with all defaults: 1 thread, build
        // the index, "A" auto-detect libtype, no extras.
        let d = tempdir("salmon");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "transcripts.fa"
index_dir     = "salmon_index"
reads         = ["reads.fq.gz"]
output_dir    = "salmon_quant"
"#,
        )
        .unwrap();
        let input = SalmonInput::from_case_dir(&d).unwrap();
        assert_eq!(input.transcriptome, PathBuf::from("transcripts.fa"));
        assert_eq!(input.index_dir, PathBuf::from("salmon_index"));
        assert_eq!(input.reads, vec![PathBuf::from("reads.fq.gz")]);
        assert_eq!(input.output_dir, PathBuf::from("salmon_quant"));
        assert_eq!(input.threads, 1);
        assert!(!input.skip_index);
        assert_eq!(input.libtype, "A");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_paired_reads() {
        // Paired-end, 8 threads, skip index (pre-built), all extras.
        let d = tempdir("salmon");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "transcripts.fa"
index_dir     = "salmon_index"
reads         = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
output_dir    = "salmon_quant"
threads       = 8
skip_index    = true
extra_args    = ["--validateMappings"]
"#,
        )
        .unwrap();
        let input = SalmonInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.reads[1], PathBuf::from("reads_R2.fq.gz"));
        assert_eq!(input.threads, 8);
        assert!(input.skip_index);
        assert_eq!(input.libtype, "A");
        assert_eq!(input.extra_args, vec!["--validateMappings".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_libtype_isf() {
        // Stranded paired-end Illumina TruSeq: ISF (Inward, Stranded,
        // Forward). Most common explicit libtype in real RNA-seq
        // pipelines.
        let d = tempdir("salmon");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "transcripts.fa"
index_dir     = "salmon_index"
reads         = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
output_dir    = "salmon_quant"
libtype       = "ISF"
"#,
        )
        .unwrap();
        let input = SalmonInput::from_case_dir(&d).unwrap();
        assert_eq!(input.libtype, "ISF");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_three_reads() {
        // Salmon is single-end (1 FASTQ) or paired-end (2 FASTQs).
        // Three is never a valid invocation.
        let d = tempdir("salmon");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "transcripts.fa"
index_dir     = "salmon_index"
reads         = ["a.fq", "b.fq", "c.fq"]
output_dir    = "salmon_quant"
"#,
        )
        .unwrap();
        let err = SalmonInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("1 (single-end) or 2 (paired-end)"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        let d = tempdir("salmon");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "transcripts.fa"
index_dir     = "salmon_index"
reads         = ["reads.fq"]
output_dir    = "salmon_quant"
threads       = 0
"#,
        )
        .unwrap();
        let err = SalmonInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads must be >= 1"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
