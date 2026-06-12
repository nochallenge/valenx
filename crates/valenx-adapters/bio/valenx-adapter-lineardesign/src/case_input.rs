//! `[bio.lineardesign]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "lineardesign.design"
//!
//! [bio.lineardesign]
//! protein         = "spike_protein.fa"
//! output_basename = "designed_mrna"
//! lambda_param    = 1.0       # optional, defaults to 1.0
//! codon_usage     = "human"   # optional, defaults to "human"
//! extra_args      = []        # optional, defaults to []
//! ```
//!
//! LinearDesign is Baidu Research's joint codon + secondary-structure
//! mRNA design tool — given a target protein, it co-optimises codon
//! adaptation index (CAI) against minimum free energy (MFE) of the
//! resulting mRNA's secondary structure, with a single tunable
//! `--lambda` knob trading the two off. The design workhorse behind
//! the modern mRNA-vaccine era.
//!
//! `lambda_param` is named `lambda_param` here and in the TOML because
//! `lambda` is a Rust reserved keyword; the CLI flag emitted by
//! `prepare()` is the upstream `--lambda`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct LinearDesignInput {
    /// Path to the protein FASTA input. LinearDesign reads this file
    /// in place via `--aa <path>`; relative paths resolve against the
    /// case directory but the file is NOT staged into the workdir.
    pub protein: PathBuf,
    /// Filename stem for outputs. LinearDesign writes
    /// `<basename>*.fasta` (designed mRNA) and `<basename>*.txt`
    /// (CAI/MFE report) into the workdir.
    pub output_basename: String,
    /// CAI vs MFE tradeoff. 0.0 minimises MFE only (most-stable
    /// structure regardless of codon usage); large values approach
    /// pure CAI optimisation. Default 1.0 — balanced.
    pub lambda_param: f64,
    /// Codon-usage table key. LinearDesign ships with `"human"`,
    /// `"yeast"`, `"e_coli"`, etc.; the user can also supply their
    /// own table by path. Default `"human"`.
    pub codon_usage: String,
    /// Additional CLI arguments appended to the lineardesign
    /// invocation. Useful for `--verbose`, `--beamsize <N>`, or any
    /// future flags the upstream tool grows.
    pub extra_args: Vec<String>,
}

impl LinearDesignInput {
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
            .and_then(|v| v.get("lineardesign"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.lineardesign] section",
                    case_toml.display()
                ))
            })?;

        let protein = block
            .get("protein")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.lineardesign].protein required (path to FASTA protein file)"
                ))
            })?;
        if protein.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.lineardesign].protein must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.lineardesign].output_basename required"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.lineardesign].output_basename must not be empty"
            )));
        }

        // `lambda_param` is the CAI/MFE tradeoff knob — admit integers
        // (e.g. `lambda_param = 1`) and validate finite + non-negative.
        // 0.0 is a legitimate value (pure MFE optimisation).
        let lambda_param = match block.get("lambda_param") {
            Some(v) => {
                let raw = if let Some(f) = v.as_float() {
                    f
                } else if let Some(i) = v.as_integer() {
                    i as f64
                } else {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.lineardesign].lambda_param must be a number"
                    )));
                };
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.lineardesign].lambda_param must be finite, got {raw}"
                    )));
                }
                if raw < 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.lineardesign].lambda_param must be >= 0.0, got {raw}"
                    )));
                }
                raw
            }
            None => 1.0,
        };

        let codon_usage = block
            .get("codon_usage")
            .and_then(|v| v.as_str())
            .unwrap_or("human")
            .to_string();
        if codon_usage.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.lineardesign].codon_usage must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.lineardesign].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.lineardesign].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            protein: PathBuf::from(protein),
            output_basename: output_basename.to_string(),
            lambda_param,
            codon_usage,
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
        // Only the two required keys: protein and output_basename.
        // Defaults: lambda_param 1.0, codon_usage "human", no extras.
        let d = tempdir("lineardesign-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "lineardesign.design"

[bio.lineardesign]
protein         = "spike_protein.fa"
output_basename = "designed_mrna"
"#,
        )
        .unwrap();
        let input = LinearDesignInput::from_case_dir(&d).unwrap();
        assert_eq!(input.protein, PathBuf::from("spike_protein.fa"));
        assert_eq!(input.output_basename, "designed_mrna");
        assert_eq!(input.lambda_param, 1.0);
        assert_eq!(input.codon_usage, "human");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_overrides() {
        // Custom lambda (heavy CAI weighting), explicit codon usage,
        // and a verbose extra to show the multi-flag composition.
        let d = tempdir("lineardesign-full");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "lineardesign.design"

[bio.lineardesign]
protein         = "spike.fa"
output_basename = "spike-design"
lambda_param    = 3.0
codon_usage     = "yeast"
extra_args      = ["--verbose"]
"#,
        )
        .unwrap();
        let input = LinearDesignInput::from_case_dir(&d).unwrap();
        assert_eq!(input.protein, PathBuf::from("spike.fa"));
        assert_eq!(input.output_basename, "spike-design");
        assert_eq!(input.lambda_param, 3.0);
        assert_eq!(input.codon_usage, "yeast");
        assert_eq!(input.extra_args, vec!["--verbose".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_negative_lambda() {
        // The CAI/MFE tradeoff is non-negative by construction —
        // negative weights would invert the optimisation direction
        // and yield nonsense designs. Reject up front.
        let d = tempdir("lineardesign-neg");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "lineardesign.design"

[bio.lineardesign]
protein         = "spike.fa"
output_basename = "designed"
lambda_param    = -1.0
"#,
        )
        .unwrap();
        let err = LinearDesignInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("lambda_param"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
