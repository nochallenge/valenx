//! `[bio.badread]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "badread.simulate"
//!
//! [bio.badread]
//! reference     = "ref.fa"
//! output        = "reads.fq"
//! quantity      = "100M"                    # bytes of sequence to simulate
//! error_model   = "nanopore2023"            # one of nanopore2018, nanopore2020, nanopore2023, pacbio2016
//! identity_mean = 87.5                      # optional, defaults to 87.5
//! length_mean   = 15000.0                   # optional, defaults to 15000.0
//! length_sd     = 13000.0                   # optional, defaults to 13000.0
//! extra_args    = ["--seed", "42"]          # optional, defaults to []
//! ```
//!
//! `quantity` is Badread's `--quantity` argument: a count of bases
//! with an optional SI suffix (e.g. `100M` for 100 megabases, `5G`
//! for 5 gigabases). The format mirrors the underlying Badread CLI.
//!
//! `error_model` selects the per-platform error profile baked into
//! the Badread distribution:
//!
//! - `nanopore2018` — historical R9.4 chemistry, ~85% identity
//! - `nanopore2020` — late-2020 R9.4.1 chemistry, ~92% identity
//! - `nanopore2023` — R10.4.1 chemistry, ~99% identity
//! - `pacbio2016`  — historical PacBio CLR, ~87% identity

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical Badread error-model list. Module-public so the adapter
/// can surface the supported values to the UI without redefining
/// them.
pub const SUPPORTED_ERROR_MODELS: &[&str] =
    &["nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016"];

/// Returns true if `s` is a Badread `--quantity` literal: one or
/// more decimal digits followed by an optional `K`, `M`, `G`, or
/// `T` SI suffix. Empty string fails; whitespace fails; lower-case
/// suffix fails (Badread itself is case-sensitive on the CLI).
pub fn is_valid_quantity(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }
    let bytes = s.as_bytes();
    let last = bytes.len() - 1;
    let (digits_end, has_suffix) = match bytes[last] {
        b'K' | b'M' | b'G' | b'T' => (last, true),
        _ => (bytes.len(), false),
    };
    if digits_end == 0 {
        // Suffix-only (e.g. "M") is not a valid quantity.
        return false;
    }
    if has_suffix && digits_end == 0 {
        return false;
    }
    bytes[..digits_end].iter().all(|b| b.is_ascii_digit())
}

#[derive(Clone, Debug, PartialEq)]
pub struct BadreadInput {
    pub reference: PathBuf,
    pub output: PathBuf,
    pub quantity: String,
    pub error_model: String,
    pub identity_mean: f64,
    pub length_mean: f64,
    pub length_sd: f64,
    pub extra_args: Vec<String>,
}

impl BadreadInput {
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
            .and_then(|v| v.get("badread"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.badread] section",
                    case_toml.display()
                ))
            })?;

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.badread].reference required (path to reference FASTA)"
                ))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.badread].reference must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.badread].output required (path for the simulated FASTQ)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.badread].output must not be empty"
            )));
        }

        let quantity = block
            .get("quantity")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.badread].quantity required (e.g. \"100M\", \"5G\")"
                ))
            })?;
        if !is_valid_quantity(quantity) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.badread].quantity `{quantity}` invalid — expected \
                 a positive integer optionally followed by K, M, G, or T \
                 (e.g. \"100M\", \"5G\")"
            )));
        }

        let error_model = block
            .get("error_model")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.badread].error_model required (one of {SUPPORTED_ERROR_MODELS:?})"
                ))
            })?;
        if !SUPPORTED_ERROR_MODELS.contains(&error_model) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.badread].error_model `{error_model}` not recognised — \
                 expected one of {SUPPORTED_ERROR_MODELS:?}"
            )));
        }

        let identity_mean = match block.get("identity_mean") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.badread].identity_mean must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].identity_mean must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].identity_mean must be > 0.0, got {raw}"
                    )));
                }
                raw
            }
            None => 87.5,
        };

        let length_mean = match block.get("length_mean") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.badread].length_mean must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].length_mean must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].length_mean must be > 0.0, got {raw}"
                    )));
                }
                raw
            }
            None => 15000.0,
        };

        let length_sd = match block.get("length_sd") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.badread].length_sd must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].length_sd must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].length_sd must be > 0.0, got {raw}"
                    )));
                }
                raw
            }
            None => 13000.0,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.badread].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.badread].extra_args entries must be strings"
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
            output: PathBuf::from(output_str),
            quantity: quantity.to_string(),
            error_model: error_model.to_string(),
            identity_mean,
            length_mean,
            length_sd,
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
        // The smallest reasonable config: a reference, an output
        // FASTQ, the quantity, and the error model. Length and
        // identity defaults track Badread's modern Nanopore
        // assumptions (~99% identity, ~15 kb mean length).
        let d = tempdir("badread");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "badread.simulate"

[bio.badread]
reference   = "ref.fa"
output      = "reads.fq"
quantity    = "100M"
error_model = "nanopore2023"
"#,
        )
        .unwrap();
        let input = BadreadInput::from_case_dir(&d).unwrap();
        assert_eq!(input.reference, PathBuf::from("ref.fa"));
        assert_eq!(input.output, PathBuf::from("reads.fq"));
        assert_eq!(input.quantity, "100M");
        assert_eq!(input.error_model, "nanopore2023");
        assert_eq!(input.identity_mean, 87.5);
        assert_eq!(input.length_mean, 15000.0);
        assert_eq!(input.length_sd, 13000.0);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Older R9.4 chemistry at high coverage with a tighter
        // length distribution and a fixed seed for reproducibility.
        let d = tempdir("badread");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "badread.simulate"

[bio.badread]
reference     = "ref.fa"
output        = "reads.fq"
quantity      = "5G"
error_model   = "nanopore2018"
identity_mean = 85.0
length_mean   = 8000.0
length_sd     = 3000.0
extra_args    = ["--seed", "42"]
"#,
        )
        .unwrap();
        let input = BadreadInput::from_case_dir(&d).unwrap();
        assert_eq!(input.quantity, "5G");
        assert_eq!(input.error_model, "nanopore2018");
        assert_eq!(input.identity_mean, 85.0);
        assert_eq!(input.length_mean, 8000.0);
        assert_eq!(input.length_sd, 3000.0);
        assert_eq!(
            input.extra_args,
            vec!["--seed".to_string(), "42".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_error_model() {
        // A typo'd "nanopore2024" (or a wishful-thinking future
        // model) should fail rather than be silently passed to
        // Badread, which only ships the four canonical profiles.
        let d = tempdir("badread");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "badread.simulate"

[bio.badread]
reference   = "ref.fa"
output      = "reads.fq"
quantity    = "100M"
error_model = "nanopore2024"
"#,
        )
        .unwrap();
        let err = BadreadInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("nanopore2023"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_invalid_quantity() {
        // Lower-case "m" doesn't parse — Badread's CLI is case-
        // sensitive on the SI suffix; reject pre-spawn.
        let d = tempdir("badread");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "badread.simulate"

[bio.badread]
reference   = "ref.fa"
output      = "reads.fq"
quantity    = "100m"
error_model = "nanopore2023"
"#,
        )
        .unwrap();
        let err = BadreadInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("quantity"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_identity_mean() {
        // Zero identity is meaningless — every base would be wrong.
        let d = tempdir("badread");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "badread.simulate"

[bio.badread]
reference     = "ref.fa"
output        = "reads.fq"
quantity      = "100M"
error_model   = "nanopore2023"
identity_mean = 0.0
"#,
        )
        .unwrap();
        let err = BadreadInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("identity_mean"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
