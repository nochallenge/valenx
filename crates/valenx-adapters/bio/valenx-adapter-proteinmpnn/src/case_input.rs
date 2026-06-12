//! `[bio.proteinmpnn]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "proteinmpnn.design"
//!
//! [bio.proteinmpnn]
//! script             = "design_proteinmpnn.py"
//! python             = "python3"          # optional, default python3
//! input_pdb          = "backbone.pdb"
//! model_variant      = "vanilla"          # vanilla | soluble | ca-only
//! temperature        = 0.1                # optional, default 0.1, must be > 0 and finite
//! num_seq_per_target = 8                  # optional, default 8, must be >= 1
//! output_basename    = "design"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ProteinMpnnInput {
    pub script: PathBuf,
    pub python: String,
    pub input_pdb: PathBuf,
    /// One of [`MODEL_VARIANTS`] — picks which ProteinMPNN weights
    /// the user script should load.
    pub model_variant: String,
    /// Sampling temperature. Defaults to 0.1; lower values make the
    /// sequences less diverse.
    pub temperature: f64,
    /// How many sequences to sample per target chain. Defaults to 8.
    pub num_seq_per_target: u32,
    /// Stem the user script should write outputs under (typically
    /// `{output_basename}.fa`).
    pub output_basename: String,
}

const DEFAULT_TEMPERATURE: f64 = 0.1;
const DEFAULT_NUM_SEQ_PER_TARGET: u32 = 8;

/// Recognised `model_variant` values. `ca-only` is ProteinMPNN's
/// alpha-carbon-only weights for backbone-light inputs; `soluble`
/// biases for soluble-protein design; `vanilla` is the default.
pub const MODEL_VARIANTS: &[&str] = &["vanilla", "soluble", "ca-only"];

impl ProteinMpnnInput {
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
            .and_then(|v| v.get("proteinmpnn"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.proteinmpnn] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.proteinmpnn].script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.proteinmpnn].input_pdb required"))
            })?;
        let model_variant = block
            .get("model_variant")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.proteinmpnn].model_variant required"))
            })?;
        if !MODEL_VARIANTS.contains(&model_variant) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.proteinmpnn].model_variant `{model_variant}` not recognised; \
                 must be one of {MODEL_VARIANTS:?}"
            )));
        }
        let temperature = block
            .get("temperature")
            .and_then(|v| v.as_float())
            .or_else(|| {
                // TOML parses `0.5` as float and `1` as integer; allow
                // integer literals too so `temperature = 1` works.
                block
                    .get("temperature")
                    .and_then(|v| v.as_integer())
                    .map(|n| n as f64)
            })
            .unwrap_or(DEFAULT_TEMPERATURE);
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.proteinmpnn].temperature must be > 0 and finite, got {temperature}"
            )));
        }
        let num_seq_per_target = block
            .get("num_seq_per_target")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_SEQ_PER_TARGET);
        if num_seq_per_target < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.proteinmpnn].num_seq_per_target must be >= 1, got {num_seq_per_target}"
            )));
        }
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.proteinmpnn].output_basename required"
                ))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.proteinmpnn].output_basename must be non-empty"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_pdb: PathBuf::from(input_pdb),
            model_variant: model_variant.to_string(),
            temperature,
            num_seq_per_target,
            output_basename: output_basename.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("proteinmpnn-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "proteinmpnn.design"

[bio.proteinmpnn]
script          = "design.py"
input_pdb       = "backbone.pdb"
model_variant   = "vanilla"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = ProteinMpnnInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.input_pdb, PathBuf::from("backbone.pdb"));
        assert_eq!(input.model_variant, "vanilla");
        assert_eq!(input.output_basename, "design");
        // Defaults.
        assert_eq!(input.python, "python3");
        assert!((input.temperature - 0.1).abs() < f64::EPSILON);
        assert_eq!(input.num_seq_per_target, 8);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_soluble_variant() {
        let d = tempdir("proteinmpnn-soluble");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "proteinmpnn.design"

[bio.proteinmpnn]
script             = "design.py"
python             = "/opt/conda/bin/python"
input_pdb          = "backbone.pdb"
model_variant      = "soluble"
temperature        = 0.5
num_seq_per_target = 32
output_basename    = "soluble_run"
"#,
        )
        .unwrap();
        let input = ProteinMpnnInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model_variant, "soluble");
        assert!((input.temperature - 0.5).abs() < f64::EPSILON);
        assert_eq!(input.num_seq_per_target, 32);
        assert_eq!(input.python, "/opt/conda/bin/python");
        assert_eq!(input.output_basename, "soluble_run");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_variant() {
        let d = tempdir("proteinmpnn-badvariant");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "proteinmpnn.design"

[bio.proteinmpnn]
script          = "design.py"
input_pdb       = "backbone.pdb"
model_variant   = "made_up_variant"
output_basename = "design"
"#,
        )
        .unwrap();
        let err = ProteinMpnnInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("model_variant"));
        assert!(format!("{err}").contains("made_up_variant"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_temperature() {
        let d = tempdir("proteinmpnn-zerotemp");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "proteinmpnn.design"

[bio.proteinmpnn]
script          = "design.py"
input_pdb       = "backbone.pdb"
model_variant   = "vanilla"
temperature     = 0.0
output_basename = "design"
"#,
        )
        .unwrap();
        let err = ProteinMpnnInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("temperature"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
