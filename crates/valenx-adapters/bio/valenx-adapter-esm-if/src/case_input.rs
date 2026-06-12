//! `[bio.esm-if]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "esm-if.design"
//!
//! [bio.esm-if]
//! script          = "design_esmif.py"
//! python          = "python3"          # optional, default python3
//! input_pdb       = "backbone.pdb"
//! model           = "esm_if1_gvp4_t16_142M_UR50"
//! temperature     = 1.0                # optional, default 1.0, must be > 0 and finite
//! num_samples     = 8                  # optional, default 8, must be >= 1
//! output_basename = "design"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct EsmIfInput {
    pub script: PathBuf,
    pub python: String,
    pub input_pdb: PathBuf,
    /// Model identifier passed through to the user script (no
    /// whitelist — ESM-IF model names evolve fast, so we just
    /// require non-empty).
    pub model: String,
    /// Sampling temperature. Defaults to 1.0; lower values produce
    /// higher-likelihood, less-diverse sequences.
    pub temperature: f64,
    /// How many sequences to sample per input. Defaults to 8.
    pub num_samples: u32,
    /// Stem the user script should write outputs under (typically
    /// `{output_basename}.fa`).
    pub output_basename: String,
}

const DEFAULT_TEMPERATURE: f64 = 1.0;
const DEFAULT_NUM_SAMPLES: u32 = 8;

impl EsmIfInput {
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
            .and_then(|v| v.get("esm-if"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.esm-if] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esm-if].script required")))?;
        if script.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm-if].script must be non-empty"
            )));
        }
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esm-if].input_pdb required"))
            })?;
        if input_pdb.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm-if].input_pdb must be non-empty"
            )));
        }
        let model = block
            .get("model")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esm-if].model required")))?;
        if model.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm-if].model must be non-empty"
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
                "[bio.esm-if].temperature must be > 0 and finite, got {temperature}"
            )));
        }
        let num_samples = block
            .get("num_samples")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_SAMPLES);
        if num_samples < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm-if].num_samples must be >= 1, got {num_samples}"
            )));
        }
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esm-if].output_basename required"))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm-if].output_basename must be non-empty"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_pdb: PathBuf::from(input_pdb),
            model: model.to_string(),
            temperature,
            num_samples,
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
        let d = tempdir("esm-if-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm-if.design"

[bio.esm-if]
script          = "design.py"
input_pdb       = "backbone.pdb"
model           = "esm_if1_gvp4_t16_142M_UR50"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = EsmIfInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.input_pdb, PathBuf::from("backbone.pdb"));
        assert_eq!(input.model, "esm_if1_gvp4_t16_142M_UR50");
        assert_eq!(input.output_basename, "design");
        // Defaults.
        assert_eq!(input.python, "python3");
        assert!((input.temperature - 1.0).abs() < f64::EPSILON);
        assert_eq!(input.num_samples, 8);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_temperature_override() {
        let d = tempdir("esm-if-temp");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm-if.design"

[bio.esm-if]
script          = "design.py"
python          = "/opt/conda/bin/python"
input_pdb       = "backbone.pdb"
model           = "esm_if1_gvp4_t16_142M_UR50"
temperature     = 0.3
num_samples     = 32
output_basename = "esm_run"
"#,
        )
        .unwrap();
        let input = EsmIfInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/bin/python");
        assert!((input.temperature - 0.3).abs() < f64::EPSILON);
        assert_eq!(input.num_samples, 32);
        assert_eq!(input.output_basename, "esm_run");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_model() {
        let d = tempdir("esm-if-empty_model");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm-if.design"

[bio.esm-if]
script          = "design.py"
input_pdb       = "backbone.pdb"
model           = ""
output_basename = "design"
"#,
        )
        .unwrap();
        let err = EsmIfInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("model"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_samples() {
        let d = tempdir("esm-if-zero_samples");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm-if.design"

[bio.esm-if]
script          = "design.py"
input_pdb       = "backbone.pdb"
model           = "esm_if1_gvp4_t16_142M_UR50"
num_samples     = 0
output_basename = "design"
"#,
        )
        .unwrap();
        let err = EsmIfInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_samples"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
