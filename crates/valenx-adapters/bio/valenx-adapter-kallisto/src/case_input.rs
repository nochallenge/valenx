//! `[bio.kallisto]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "kallisto.quant"
//!
//! [bio.kallisto]
//! transcriptome   = "transcripts.fa"
//! index           = "transcripts.idx"                    # single .idx file
//! reads           = ["reads_R1.fq.gz", "reads_R2.fq.gz"] # 1 (single-end) or 2 (paired-end)
//! output_dir      = "kallisto_quant"                     # quant output dir
//! threads         = 4                                     # optional, defaults to 1
//! skip_index      = false                                 # optional, defaults to false
//! fragment_length = 200.0                                 # required for single-end (mean fragment length)
//! fragment_sd     = 20.0                                  # required for single-end (fragment-length stdev)
//! extra_args      = ["--bias"]                            # optional, defaults to []
//! ```
//!
//! Kallisto accepts a single-end FASTQ or a paired-end pair of
//! FASTQs. Unlike Salmon, kallisto's index is a single file
//! (`<name>.idx`) rather than a directory.
//!
//! For single-end reads, kallisto can't infer the fragment-length
//! distribution from the data — the user must supply
//! `fragment_length` (the mean) and `fragment_sd` (the standard
//! deviation) explicitly. Both must be positive, finite floats.
//! For paired-end reads kallisto auto-detects from the data and
//! these fields are ignored if present.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct KallistoInput {
    pub transcriptome: PathBuf,
    pub index: PathBuf,
    pub reads: Vec<PathBuf>,
    pub output_dir: PathBuf,
    pub threads: u32,
    pub skip_index: bool,
    pub fragment_length: Option<f64>,
    pub fragment_sd: Option<f64>,
    pub extra_args: Vec<String>,
}

impl KallistoInput {
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
            .and_then(|v| v.get("kallisto"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.kallisto] section",
                    case_toml.display()
                ))
            })?;

        let transcriptome_str = block
            .get("transcriptome")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.kallisto].transcriptome required"))
            })?;
        if transcriptome_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.kallisto].transcriptome must not be empty"
            )));
        }

        let index_str = block
            .get("index")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.kallisto].index required")))?;
        if index_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.kallisto].index must not be empty"
            )));
        }

        let reads_arr = block
            .get("reads")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].reads required (array of 1 or 2 FASTQ paths)"
                ))
            })?;
        let mut reads: Vec<PathBuf> = Vec::with_capacity(reads_arr.len());
        for entry in reads_arr {
            let s = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].reads entries must be strings"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].reads entries must not be empty"
                )));
            }
            reads.push(PathBuf::from(s));
        }
        if reads.is_empty() || reads.len() > 2 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.kallisto].reads must contain 1 (single-end) or 2 \
                 (paired-end) FASTQs, got {}",
                reads.len()
            )));
        }

        let output_dir_str = block
            .get("output_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.kallisto].output_dir required"))
            })?;
        if output_dir_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.kallisto].output_dir must not be empty"
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.kallisto].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.kallisto].threads must be >= 1, got {raw}"
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

        // `fragment_length` / `fragment_sd` accept either an integer
        // or a float in TOML — kallisto's flag is a float so we
        // canonicalise to f64 here.
        let fragment_length = match block.get("fragment_length") {
            Some(v) => {
                let f = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.kallisto].fragment_length must be a number"
                        ))
                    })?;
                Some(f)
            }
            None => None,
        };
        let fragment_sd = match block.get("fragment_sd") {
            Some(v) => {
                let f = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.kallisto].fragment_sd must be a number"
                        ))
                    })?;
                Some(f)
            }
            None => None,
        };

        // Single-end runs need explicit fragment-length stats —
        // kallisto can't auto-detect them from a single FASTQ. Both
        // values must be positive, finite floats.
        if reads.len() == 1 {
            let l = fragment_length.ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].fragment_length is required for single-end \
                     runs (kallisto can't auto-detect from a single FASTQ)"
                ))
            })?;
            if !(l.is_finite() && l > 0.0) {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].fragment_length must be a positive finite \
                     number, got {l}"
                )));
            }
            let s = fragment_sd.ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].fragment_sd is required for single-end \
                     runs (kallisto can't auto-detect from a single FASTQ)"
                ))
            })?;
            if !(s.is_finite() && s > 0.0) {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.kallisto].fragment_sd must be a positive finite \
                     number, got {s}"
                )));
            }
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.kallisto].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.kallisto].extra_args entries must be strings"
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
            index: PathBuf::from(index_str),
            reads,
            output_dir: PathBuf::from(output_dir_str),
            threads,
            skip_index,
            fragment_length,
            fragment_sd,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_paired_minimal() {
        // Paired-end with all defaults: 1 thread, build the index, no
        // fragment stats (auto-detect from paired reads), no extras.
        let d = tempdir("kallisto");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome = "transcripts.fa"
index         = "transcripts.idx"
reads         = ["reads_R1.fq.gz", "reads_R2.fq.gz"]
output_dir    = "kallisto_quant"
"#,
        )
        .unwrap();
        let input = KallistoInput::from_case_dir(&d).unwrap();
        assert_eq!(input.transcriptome, PathBuf::from("transcripts.fa"));
        assert_eq!(input.index, PathBuf::from("transcripts.idx"));
        assert_eq!(input.reads.len(), 2);
        assert_eq!(input.output_dir, PathBuf::from("kallisto_quant"));
        assert_eq!(input.threads, 1);
        assert!(!input.skip_index);
        assert!(input.fragment_length.is_none());
        assert!(input.fragment_sd.is_none());
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_single_end_with_fragment_stats() {
        // Single-end with required fragment stats and a typical
        // bias-correction extra arg.
        let d = tempdir("kallisto");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome   = "transcripts.fa"
index           = "transcripts.idx"
reads           = ["reads.fq.gz"]
output_dir      = "kallisto_quant"
threads         = 4
skip_index      = true
fragment_length = 200.0
fragment_sd     = 20.0
extra_args      = ["--bias"]
"#,
        )
        .unwrap();
        let input = KallistoInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reads.len(), 1);
        assert_eq!(input.threads, 4);
        assert!(input.skip_index);
        assert_eq!(input.fragment_length, Some(200.0));
        assert_eq!(input.fragment_sd, Some(20.0));
        assert_eq!(input.extra_args, vec!["--bias".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_single_end_without_fragment_length() {
        // Single-end without `fragment_length` is incoherent —
        // kallisto can't pseudoalign a single FASTQ without knowing
        // the mean fragment length.
        let d = tempdir("kallisto");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome = "transcripts.fa"
index         = "transcripts.idx"
reads         = ["reads.fq.gz"]
output_dir    = "kallisto_quant"
fragment_sd   = 20.0
"#,
        )
        .unwrap();
        let err = KallistoInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("fragment_length is required for single-end"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_single_end_without_fragment_sd() {
        // Same story — kallisto needs both the mean *and* the stdev
        // for single-end pseudoalignment.
        let d = tempdir("kallisto");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome   = "transcripts.fa"
index           = "transcripts.idx"
reads           = ["reads.fq.gz"]
output_dir      = "kallisto_quant"
fragment_length = 200.0
"#,
        )
        .unwrap();
        let err = KallistoInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("fragment_sd is required for single-end"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_three_reads() {
        // Kallisto is single-end (1 FASTQ) or paired-end (2 FASTQs).
        // Three is never a valid invocation.
        let d = tempdir("kallisto");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome = "transcripts.fa"
index         = "transcripts.idx"
reads         = ["a.fq", "b.fq", "c.fq"]
output_dir    = "kallisto_quant"
"#,
        )
        .unwrap();
        let err = KallistoInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("1 (single-end) or 2 (paired-end)"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
