//! `[bio.bwa]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "bwa.mem"
//!
//! [bio.bwa]
//! reference  = "ref.fa"
//! reads      = ["reads_R1.fq.gz", "reads_R2.fq.gz"]   # 1 (single-end) or 2 (paired-end)
//! threads    = 4                                       # optional, defaults to 1
//! skip_index = false                                   # optional, defaults to false
//! extra_args = ["-M"]                                  # optional, defaults to []
//! ```
//!
//! BWA-MEM accepts either a single-end FASTQ or a paired-end pair of
//! FASTQs; it does **not** support 3+ inputs as a single invocation.
//! `skip_index = true` lets the user re-use a pre-built BWT index
//! sitting next to the reference (`ref.fa.bwt`, `.pac`, `.ann`,
//! `.amb`, `.sa`) so successive runs over the same reference don't
//! pay the index-build cost.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct BwaInput {
    pub reference: PathBuf,
    pub reads: Vec<PathBuf>,
    pub threads: u32,
    pub skip_index: bool,
    pub extra_args: Vec<String>,
}

impl BwaInput {
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
            .and_then(|v| v.get("bwa"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.bwa] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.bwa].reference required")))?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bwa].reference must not be empty"
            )));
        }

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.bwa].reads required (array of 1 or 2 FASTQ paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.bwa].reads entries must be strings"))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.bwa].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() || reads.len() > 2 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bwa].reads must contain 1 (single-end) or 2 \
                 (paired-end) FASTQs, got {}",
                reads.len()
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.bwa].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.bwa].threads must be >= 1, got {raw}"
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

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.bwa].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.bwa].extra_args entries must be strings"
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
        let d = tempdir("bwa");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bwa.mem"

[bio.bwa]
reference = "ref.fa"
reads     = ["reads_R1.fq.gz"]
"#,
        )
        .unwrap();
        let input = BwaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.reads, vec![PathBuf::from("reads_R1.fq.gz")]);
        // Defaults: 1 thread, build the index, no extra args.
        assert_eq!(input.threads, 1);
        assert!(!input.skip_index);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_paired_reads_and_extra_args() {
        // Paired-end alignment with explicit threading, a pre-built
        // index reused from a prior run, and a typical Picard-friendly
        // BWA-MEM flag (`-M` marks shorter split hits as secondary so
        // downstream tools don't choke).
        let d = tempdir("bwa");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bwa.mem"

[bio.bwa]
reference  = "ref.fa"
reads      = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
threads    = 8
skip_index = true
extra_args = ["-M", "-K", "100000000"]
"#,
        )
        .unwrap();
        let input = BwaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.reads[0], PathBuf::from("reads_R1.fq.gz"));
        assert_eq!(input.reads[1], PathBuf::from("reads_R2.fq.gz"));
        assert_eq!(input.threads, 8);
        assert!(input.skip_index);
        assert_eq!(
            input.extra_args,
            vec!["-M".to_string(), "-K".to_string(), "100000000".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_three_reads() {
        // BWA-MEM is single-end (1 FASTQ) or paired-end (2 FASTQs).
        // Three is never a valid invocation.
        let d = tempdir("bwa");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bwa.mem"

[bio.bwa]
reference = "ref.fa"
reads     = ["a.fq", "b.fq", "c.fq"]
"#,
        )
        .unwrap();
        let err = BwaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("1 (single-end) or 2 (paired-end)"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("bwa");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = BwaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.bwa]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
