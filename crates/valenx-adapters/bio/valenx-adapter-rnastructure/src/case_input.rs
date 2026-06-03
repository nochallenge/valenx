//! `[bio.rnastructure]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "rnastructure.fold"
//!
//! [bio.rnastructure]
//! input          = "rna.seq"
//! output         = "fold.ct"
//! max_structures = 20             # optional, defaults to 20
//! max_percent    = 10             # optional, defaults to 10 (percent)
//! temperature    = 310.15         # optional, defaults to 310.15 K (37 C)
//! extra_args     = ["--maxlength", "30"]
//! ```
//!
//! `RNAstructure`'s `Fold` reads a `.seq` (or `.fa`) sequence file
//! and writes a connectivity table (`.ct`) holding the predicted MFE
//! plus suboptimal structures. `max_structures` caps the number of
//! suboptimals returned; `max_percent` caps how far above MFE they
//! may sit (in % of the MFE). `temperature` is in Kelvin — the
//! Mathews-lab convention; the upstream Fold's default of 310.15 K
//! (37 °C) is preserved as our default.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RnaStructureInput {
    pub input: PathBuf,
    pub output: PathBuf,
    pub max_structures: u32,
    pub max_percent: u32,
    pub temperature: f64,
    pub extra_args: Vec<String>,
}

impl RnaStructureInput {
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
            .and_then(|v| v.get("rnastructure"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.rnastructure] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.rnastructure].input required (path to sequence file)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rnastructure].input must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.rnastructure].output required (CT-format output path)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rnastructure].output must not be empty"
            )));
        }

        let max_structures = match block.get("max_structures") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].max_structures must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].max_structures must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].max_structures `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 20,
        };

        let max_percent = match block.get("max_percent") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].max_percent must be an integer"
                    ))
                })?;
                if !(0..=100).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].max_percent must be in 0..=100, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 10,
        };

        let temperature = match block.get("temperature") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.rnastructure].temperature must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].temperature must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].temperature must be > 0.0 (Kelvin), got {raw}"
                    )));
                }
                raw
            }
            None => 310.15,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.rnastructure].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.rnastructure].extra_args entries must be strings"
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
            output: PathBuf::from(output_str),
            max_structures,
            max_percent,
            temperature,
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
        // Minimum-config case: input + output. Defaults must match
        // documented values (20 structures, 10%, 310.15 K).
        let d = tempdir("rnastructure");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rnastructure.fold"

[bio.rnastructure]
input  = "rna.seq"
output = "fold.ct"
"#,
        )
        .unwrap();
        let input = RnaStructureInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("rna.seq"));
        assert_eq!(input.output, PathBuf::from("fold.ct"));
        assert_eq!(input.max_structures, 20);
        assert_eq!(input.max_percent, 10);
        assert!((input.temperature - 310.15).abs() < 1e-9);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // All four knobs explicitly set, plus extras. Exercises the
        // integer / float parsing paths together.
        let d = tempdir("rnastructure");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rnastructure.fold"

[bio.rnastructure]
input          = "rna.seq"
output         = "fold.ct"
max_structures = 5
max_percent    = 25
temperature    = 298.15
extra_args     = ["--maxlength", "30"]
"#,
        )
        .unwrap();
        let input = RnaStructureInput::from_case_dir(&d).unwrap();
        assert_eq!(input.max_structures, 5);
        assert_eq!(input.max_percent, 25);
        assert!((input.temperature - 298.15).abs() < 1e-9);
        assert_eq!(
            input.extra_args,
            vec!["--maxlength".to_string(), "30".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_max_structures() {
        // 0 structures is a no-op invocation — Fold would refuse or
        // emit an empty CT. Reject up front so the user gets a clear
        // error from the adapter rather than an opaque Fold failure.
        let d = tempdir("rnastructure");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rnastructure.fold"

[bio.rnastructure]
input          = "rna.seq"
output         = "fold.ct"
max_structures = 0
"#,
        )
        .unwrap();
        let err = RnaStructureInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("max_structures"), "msg: {msg}");
        assert!(msg.contains(">= 1"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_max_percent_above_100() {
        // Percentages above 100 are nonsensical for "% above MFE".
        let d = tempdir("rnastructure");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rnastructure.fold"

[bio.rnastructure]
input       = "rna.seq"
output      = "fold.ct"
max_percent = 150
"#,
        )
        .unwrap();
        let err = RnaStructureInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("max_percent"), "msg: {msg}");
        assert!(msg.contains("0..=100"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
