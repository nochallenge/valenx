//! `[bio.wgsim]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "wgsim.simulate"
//!
//! [bio.wgsim]
//! reference     = "ref.fa"
//! output1       = "reads_R1.fq"
//! output2       = "reads_R2.fq"
//! num_pairs     = 1000000
//! length1       = 70                      # optional, defaults to 70
//! length2       = 70                      # optional, defaults to 70
//! fragment_size = 500                     # optional, defaults to 500
//! error_rate    = 0.02                    # optional, defaults to 0.02
//! extra_args    = ["-r", "0.001"]         # optional, defaults to []
//! ```
//!
//! wgsim is the canonical "small + classic" short-read simulator
//! that ships with samtools. Always paired-end, always
//! position-uniform — handy for fast smoke-testing of pipelines
//! when ART's empirical error model is not required.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct WgsimInput {
    pub reference: PathBuf,
    pub output1: PathBuf,
    pub output2: PathBuf,
    pub num_pairs: u32,
    pub length1: u32,
    pub length2: u32,
    pub fragment_size: u32,
    pub error_rate: f64,
    pub extra_args: Vec<String>,
}

impl WgsimInput {
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
            .and_then(|v| v.get("wgsim"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.wgsim] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.wgsim].reference required (path to reference FASTA)"
                ))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.wgsim].reference must not be empty"
            )));
        }

        let output1_str = block
            .get("output1")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.wgsim].output1 required (path for the first read FASTQ)"
                ))
            })?;
        if output1_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.wgsim].output1 must not be empty"
            )));
        }

        let output2_str = block
            .get("output2")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.wgsim].output2 required (path for the second read FASTQ)"
                ))
            })?;
        if output2_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.wgsim].output2 must not be empty"
            )));
        }

        let num_pairs = match block.get("num_pairs") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.wgsim].num_pairs must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].num_pairs must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].num_pairs `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.wgsim].num_pairs required (number of paired reads to simulate)"
                )));
            }
        };

        let length1 = match block.get("length1") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.wgsim].length1 must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].length1 must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].length1 `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 70,
        };

        let length2 = match block.get("length2") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.wgsim].length2 must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].length2 must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].length2 `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 70,
        };

        let fragment_size = match block.get("fragment_size") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].fragment_size must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].fragment_size must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].fragment_size `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 500,
        };

        let error_rate = match block.get("error_rate") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.wgsim].error_rate must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].error_rate must be finite, got {raw}"
                    )));
                }
                if !(0.0..=1.0).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].error_rate must be in 0.0..=1.0, got {raw}"
                    )));
                }
                raw
            }
            None => 0.02,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.wgsim].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.wgsim].extra_args entries must be strings"
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
            output1: PathBuf::from(output1_str),
            output2: PathBuf::from(output2_str),
            num_pairs,
            length1,
            length2,
            fragment_size,
            error_rate,
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
        // The smallest reasonable config: a reference, two output
        // FASTQs, and a pair count. Everything else falls back to
        // wgsim's documented defaults — 70 bp reads, 500 bp inserts,
        // 2% error.
        let d = tempdir("wgsim");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "wgsim.simulate"

[bio.wgsim]
reference = "ref.fa"
output1   = "r1.fq"
output2   = "r2.fq"
num_pairs = 1000
"#,
        )
        .unwrap();
        let input = WgsimInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.output1, PathBuf::from("r1.fq"));
        assert_eq!(input.output2, PathBuf::from("r2.fq"));
        assert_eq!(input.num_pairs, 1000);
        assert_eq!(input.length1, 70);
        assert_eq!(input.length2, 70);
        assert_eq!(input.fragment_size, 500);
        assert_eq!(input.error_rate, 0.02);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Modern Illumina paired-end pattern: 150 bp reads, 350 bp
        // fragments, slightly elevated 3% error rate, extra `-r`
        // mutation rate flag.
        let d = tempdir("wgsim");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "wgsim.simulate"

[bio.wgsim]
reference     = "ref.fa"
output1       = "r1.fq"
output2       = "r2.fq"
num_pairs     = 1000000
length1       = 150
length2       = 150
fragment_size = 350
error_rate    = 0.03
extra_args    = ["-r", "0.001"]
"#,
        )
        .unwrap();
        let input = WgsimInput::from_case_dir(&d).unwrap();
        assert_eq!(input.num_pairs, 1_000_000);
        assert_eq!(input.length1, 150);
        assert_eq!(input.length2, 150);
        assert_eq!(input.fragment_size, 350);
        assert_eq!(input.error_rate, 0.03);
        assert_eq!(
            input.extra_args,
            vec!["-r".to_string(), "0.001".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_num_pairs() {
        // 0 pairs is nonsensical — would produce empty FASTQs and
        // confuse downstream pipelines. Reject up front.
        let d = tempdir("wgsim");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "wgsim.simulate"

[bio.wgsim]
reference = "ref.fa"
output1   = "r1.fq"
output2   = "r2.fq"
num_pairs = 0
"#,
        )
        .unwrap();
        let err = WgsimInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("num_pairs"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_error_rate_above_1() {
        // wgsim's `-e` is a per-base probability; it must lie in
        // [0, 1]. 1.5 would be caught by wgsim itself, but flagging
        // it pre-spawn keeps the error message friendlier.
        let d = tempdir("wgsim");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "wgsim.simulate"

[bio.wgsim]
reference  = "ref.fa"
output1    = "r1.fq"
output2    = "r2.fq"
num_pairs  = 1000
error_rate = 1.5
"#,
        )
        .unwrap();
        let err = WgsimInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("error_rate"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
