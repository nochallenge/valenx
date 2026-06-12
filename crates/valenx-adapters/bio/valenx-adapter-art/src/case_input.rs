//! `[bio.art]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "art.simulate"
//!
//! [bio.art]
//! reference          = "ref.fa"
//! output_prefix      = "sim"
//! sequencing_system  = "HS25"            # one of HS25, HSXt, MSv3, NS50, MinS
//! read_length        = 150
//! fold_coverage      = 30.0
//! paired_end         = false             # optional, defaults to false
//! fragment_mean      = 200.0             # optional (used iff paired_end)
//! fragment_sd        = 10.0              # optional (used iff paired_end)
//! extra_args         = ["-na"]           # optional, defaults to []
//! ```
//!
//! `sequencing_system` selects the ART error profile:
//!
//! - `HS25`  — HiSeq 2500 (125 / 150 bp)
//! - `HSXt`  — HiSeq X TruSeq (150 bp)
//! - `MSv3`  — MiSeq v3 (250 bp)
//! - `NS50`  — NextSeq 500 (75 bp)
//! - `MinS`  — MiniSeq TruSeq (50 bp)
//!
//! Paired-end runs additionally require `fragment_mean` and
//! `fragment_sd` (the insert-size distribution); single-end runs
//! ignore both. ART itself accepts the same `-m` / `-s` flags only
//! under `-p`.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical ART sequencing-system list. Module-public so the adapter
/// can surface the supported values to the UI without redefining
/// them.
pub const SUPPORTED_SEQUENCING_SYSTEMS: &[&str] = &["HS25", "HSXt", "MSv3", "NS50", "MinS"];

#[derive(Clone, Debug, PartialEq)]
pub struct ArtInput {
    pub reference: PathBuf,
    pub output_prefix: String,
    pub sequencing_system: String,
    pub read_length: u32,
    pub fold_coverage: f64,
    pub paired_end: bool,
    pub fragment_mean: f64,
    pub fragment_sd: f64,
    pub extra_args: Vec<String>,
}

impl ArtInput {
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
            .and_then(|v| v.get("art"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.art] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].reference required (path to reference FASTA)"
                ))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.art].reference must not be empty"
            )));
        }

        let output_prefix = block
            .get("output_prefix")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].output_prefix required (filename stem ART writes to)"
                ))
            })?;
        if output_prefix.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.art].output_prefix must not be empty"
            )));
        }

        let sequencing_system = block
            .get("sequencing_system")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].sequencing_system required (one of {SUPPORTED_SEQUENCING_SYSTEMS:?})"
                ))
            })?;
        if !SUPPORTED_SEQUENCING_SYSTEMS.contains(&sequencing_system) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.art].sequencing_system `{sequencing_system}` not recognised — \
                 expected one of {SUPPORTED_SEQUENCING_SYSTEMS:?}"
            )));
        }

        let read_length = match block.get("read_length") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.art].read_length must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.art].read_length must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.art].read_length `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].read_length required (read length in bp)"
                )));
            }
        };

        let fold_coverage = match block.get("fold_coverage") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.art].fold_coverage must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.art].fold_coverage must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.art].fold_coverage must be > 0.0, got {raw}"
                    )));
                }
                raw
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].fold_coverage required (e.g. 30.0 for 30x)"
                )));
            }
        };

        let paired_end = block
            .get("paired_end")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let fragment_mean = match block.get("fragment_mean") {
            Some(v) => v
                .as_float()
                .or_else(|| v.as_integer().map(|i| i as f64))
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.art].fragment_mean must be a number"))
                })?,
            None => 200.0,
        };
        let fragment_sd = match block.get("fragment_sd") {
            Some(v) => v
                .as_float()
                .or_else(|| v.as_integer().map(|i| i as f64))
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.art].fragment_sd must be a number"))
                })?,
            None => 10.0,
        };

        if paired_end {
            if !fragment_mean.is_finite() || fragment_mean <= 0.0 {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].fragment_mean must be finite and > 0.0 \
                     when paired_end = true, got {fragment_mean}"
                )));
            }
            if !fragment_sd.is_finite() || fragment_sd <= 0.0 {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.art].fragment_sd must be finite and > 0.0 \
                     when paired_end = true, got {fragment_sd}"
                )));
            }
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.art].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.art].extra_args entries must be strings"
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
            output_prefix: output_prefix.to_string(),
            sequencing_system: sequencing_system.to_string(),
            read_length,
            fold_coverage,
            paired_end,
            fragment_mean,
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
    fn parses_minimal_single_end() {
        // The smallest reasonable single-end config: HiSeq 2500
        // 150 bp reads at 30x coverage. paired_end defaults to false
        // so fragment_* are unused (their default values still load
        // for round-trip stability).
        let d = tempdir("art");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "art.simulate"

[bio.art]
reference         = "ref.fa"
output_prefix     = "sim"
sequencing_system = "HS25"
read_length       = 150
fold_coverage     = 30.0
"#,
        )
        .unwrap();
        let input = ArtInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.output_prefix, "sim");
        assert_eq!(input.sequencing_system, "HS25");
        assert_eq!(input.read_length, 150);
        assert_eq!(input.fold_coverage, 30.0);
        assert!(!input.paired_end);
        assert_eq!(input.fragment_mean, 200.0);
        assert_eq!(input.fragment_sd, 10.0);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_paired_end_with_fragment_stats() {
        // Paired-end MiSeq v3 250 bp at 50x with a non-default
        // insert distribution (typical TruSeq Nano: ~350 bp ± 40).
        let d = tempdir("art");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "art.simulate"

[bio.art]
reference         = "ref.fa"
output_prefix     = "sim"
sequencing_system = "MSv3"
read_length       = 250
fold_coverage     = 50.0
paired_end        = true
fragment_mean     = 350.0
fragment_sd       = 40.0
extra_args        = ["-na"]
"#,
        )
        .unwrap();
        let input = ArtInput::from_case_dir(&d).unwrap();
        assert!(input.paired_end);
        assert_eq!(input.sequencing_system, "MSv3");
        assert_eq!(input.read_length, 250);
        assert_eq!(input.fold_coverage, 50.0);
        assert_eq!(input.fragment_mean, 350.0);
        assert_eq!(input.fragment_sd, 40.0);
        assert_eq!(input.extra_args, vec!["-na".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_sequencing_system() {
        // ART's `-ss` flag accepts only the canonical short codes.
        // A familiar long name like "HiSeq2500" is not what ART
        // wants on the CLI — fail early so the user catches it.
        let d = tempdir("art");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "art.simulate"

[bio.art]
reference         = "ref.fa"
output_prefix     = "sim"
sequencing_system = "HiSeq2500"
read_length       = 150
fold_coverage     = 30.0
"#,
        )
        .unwrap();
        let err = ArtInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("HS25"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_read_length() {
        // 0 bp reads are nonsensical. Fail before we ever spawn ART.
        let d = tempdir("art");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "art.simulate"

[bio.art]
reference         = "ref.fa"
output_prefix     = "sim"
sequencing_system = "HS25"
read_length       = 0
fold_coverage     = 30.0
"#,
        )
        .unwrap();
        let err = ArtInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("read_length"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_paired_end_with_zero_fragment_mean() {
        // Paired-end mode requires a positive insert-size mean —
        // ART's `-m` won't accept 0.
        let d = tempdir("art");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "art.simulate"

[bio.art]
reference         = "ref.fa"
output_prefix     = "sim"
sequencing_system = "HS25"
read_length       = 150
fold_coverage     = 30.0
paired_end        = true
fragment_mean     = 0.0
fragment_sd       = 10.0
"#,
        )
        .unwrap();
        let err = ArtInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("fragment_mean"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
