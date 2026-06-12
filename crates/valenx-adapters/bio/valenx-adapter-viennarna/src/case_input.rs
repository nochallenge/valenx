//! `[bio.viennarna]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "viennarna.fold"
//!
//! [bio.viennarna]
//! input              = "rna.fa"
//! output             = "fold.out"
//! temperature        = 37.0           # optional, defaults to 37.0 (Celsius)
//! partition_function = false          # optional, defaults to false
//! allow_gu           = true           # optional, defaults to true
//! extra_args         = ["--noLP"]     # optional, defaults to []
//! ```
//!
//! ViennaRNA's `RNAfold` reads a FASTA-style input and writes the
//! minimum-free-energy structure (plus optional partition-function /
//! base-pair-probability output) to **stdout**. The adapter's
//! `prepare()` composes the invocation; `run()` redirects stdout to
//! `<workdir>/<output>`.
//!
//! `temperature` is the folding temperature in Celsius; `RNAfold`'s
//! `-T <C>` flag takes Celsius rather than Kelvin. `partition_function`
//! enables the `-p` flag (writes the dot-plot PostScript and partition-
//! function ensemble values). `allow_gu = false` disables non-canonical
//! G-U wobble pairing via `--noGU`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ViennaRnaInput {
    pub input: PathBuf,
    pub output: PathBuf,
    pub temperature: f64,
    pub partition_function: bool,
    pub allow_gu: bool,
    pub extra_args: Vec<String>,
}

impl ViennaRnaInput {
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
            .and_then(|v| v.get("viennarna"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.viennarna] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.viennarna].input required (path to FASTA input)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.viennarna].input must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.viennarna].output required (path for RNAfold's stdout, \
                     resolved against the workdir)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.viennarna].output must not be empty"
            )));
        }

        let temperature = match block.get("temperature") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.viennarna].temperature must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.viennarna].temperature must be finite, got {raw}"
                    )));
                }
                raw
            }
            None => 37.0,
        };

        let partition_function = block
            .get("partition_function")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let allow_gu = block
            .get("allow_gu")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.viennarna].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.viennarna].extra_args entries must be strings"
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
            temperature,
            partition_function,
            allow_gu,
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
        // Minimum-config case: just input + output. Defaults must
        // match the documented values: 37 °C, no partition function,
        // G-U pairing allowed.
        let d = tempdir("viennarna");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "viennarna.fold"

[bio.viennarna]
input  = "rna.fa"
output = "fold.out"
"#,
        )
        .unwrap();
        let input = ViennaRnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("rna.fa"));
        assert_eq!(input.output, PathBuf::from("fold.out"));
        assert_eq!(input.temperature, 37.0);
        assert!(!input.partition_function);
        assert!(input.allow_gu);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_partition_function() {
        // Partition-function mode (-p): emits the ensemble free energy
        // and the dot-plot PostScript next to the MFE structure.
        let d = tempdir("viennarna");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "viennarna.fold"

[bio.viennarna]
input              = "rna.fa"
output             = "fold.out"
partition_function = true
temperature        = 25.0
"#,
        )
        .unwrap();
        let input = ViennaRnaInput::from_case_dir(&d).unwrap();
        assert!(input.partition_function);
        assert_eq!(input.temperature, 25.0);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_no_gu() {
        // Disable non-canonical G-U wobble pairing — useful when
        // benchmarking against canonical-only thermodynamic models.
        let d = tempdir("viennarna");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "viennarna.fold"

[bio.viennarna]
input    = "rna.fa"
output   = "fold.out"
allow_gu = false
"#,
        )
        .unwrap();
        let input = ViennaRnaInput::from_case_dir(&d).unwrap();
        assert!(!input.allow_gu);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input() {
        // Empty string isn't a path — would silently end up next to
        // case.toml, not what the user wanted. Reject up front.
        let d = tempdir("viennarna");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "viennarna.fold"

[bio.viennarna]
input  = ""
output = "fold.out"
"#,
        )
        .unwrap();
        let err = ViennaRnaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input must not be empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
