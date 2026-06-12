//! `[bio.clustalo]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "clustalo.align"
//!
//! [bio.clustalo]
//! input           = "input.fa"
//! output_basename = "alignment"
//! outfmt          = "clustal"           # optional, defaults to "clustal"
//! threads         = 4                   # optional, defaults to 1
//! extra_args      = ["--iter=2"]        # optional, defaults to []
//! ```
//!
//! `outfmt` selects Clustal Omega's output-format flag (passed through
//! as `--outfmt=<value>`). The full set Clustal Omega understands is
//! `fasta`, `clustal`, `msf`, `phylip`, `selex`, `stockholm`, and
//! `vienna`; we don't gate on the list since Clustal Omega itself emits
//! a clean error if an unknown value is passed, and downstream the
//! filename extension is picked to match (defaulting to `.aln` so the
//! file always lands somewhere `collect()` can find it).
//!
//! `threads` must be >= 1; 0 would silently degrade Clustal Omega to
//! single-threaded with no warning, so reject up front.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ClustaloInput {
    pub input: PathBuf,
    pub output_basename: String,
    pub outfmt: String,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl ClustaloInput {
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
            .and_then(|v| v.get("clustalo"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.clustalo] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.clustalo].input required (path to multi-FASTA input)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.clustalo].input must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.clustalo].output_basename required \
                     (basename for the alignment output, e.g. \"alignment\")"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.clustalo].output_basename must not be empty"
            )));
        }

        let outfmt = match block.get("outfmt") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.clustalo].outfmt must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.clustalo].outfmt must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "clustal".to_string(),
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.clustalo].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.clustalo].threads must be >= 1, got {raw}"
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
                        "[bio.clustalo].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.clustalo].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input: PathBuf::from(input_str),
            output_basename: output_basename.to_string(),
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
        // Only the two required keys: input, output_basename. Defaults:
        // outfmt = "clustal", threads = 1, no extras.
        let d = tempdir("clustalo");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "clustalo.align"

[bio.clustalo]
input           = "input.fa"
output_basename = "alignment"
"#,
        )
        .unwrap();
        let input = ClustaloInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.fa"));
        assert_eq!(input.output_basename, "alignment");
        assert_eq!(input.outfmt, "clustal");
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_overrides() {
        // FASTA output, 8 threads, an iter override extra. Confirms
        // every optional key threads through cleanly.
        let d = tempdir("clustalo");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "clustalo.align"

[bio.clustalo]
input           = "input.fa"
output_basename = "aligned"
outfmt          = "fasta"
threads         = 8
extra_args      = ["--iter=2", "--full"]
"#,
        )
        .unwrap();
        let input = ClustaloInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.fa"));
        assert_eq!(input.output_basename, "aligned");
        assert_eq!(input.outfmt, "fasta");
        assert_eq!(input.threads, 8);
        assert_eq!(
            input.extra_args,
            vec!["--iter=2".to_string(), "--full".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("clustalo");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = ClustaloInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.clustalo]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
