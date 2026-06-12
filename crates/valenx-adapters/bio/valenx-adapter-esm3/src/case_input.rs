//! `[bio.esm3]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "esm3.generate"
//!
//! [bio.esm3]
//! script          = "run_esm3.py"
//! python          = "python3"          # optional, default python3
//! model_variant   = "open"             # one of "open" | "open-multimer" | "small"
//! mode            = "design"           # one of "design" | "inverse-fold" | "scaffold" | "predict"
//! num_samples     = 4                  # optional, default 4, must be >= 1
//! input_pdb       = "scaffold.pdb"     # required for inverse-fold / scaffold
//! input_fasta     = "query.fasta"      # required for predict
//! temperature     = 1.0                # optional, default 1.0, must be > 0 and finite
//! output_basename = "esm3_run"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct Esm3Input {
    pub script: PathBuf,
    pub python: String,
    pub model_variant: String,
    pub mode: String,
    pub num_samples: u32,
    pub input_pdb: Option<PathBuf>,
    pub input_fasta: Option<PathBuf>,
    pub temperature: f64,
    pub output_basename: String,
}

const DEFAULT_NUM_SAMPLES: u32 = 4;
const DEFAULT_TEMPERATURE: f64 = 1.0;

const ALLOWED_VARIANTS: &[&str] = &["open", "open-multimer", "small"];
const ALLOWED_MODES: &[&str] = &["design", "inverse-fold", "scaffold", "predict"];

impl Esm3Input {
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
            .and_then(|v| v.get("esm3"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.esm3] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esm3].script required")))?;
        if script.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].script must be non-empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let model_variant = block
            .get("model_variant")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esm3].model_variant required"))
            })?;
        if !ALLOWED_VARIANTS.contains(&model_variant) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].model_variant must be one of {ALLOWED_VARIANTS:?}, got `{model_variant}`"
            )));
        }

        let mode = block
            .get("mode")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esm3].mode required")))?;
        if !ALLOWED_MODES.contains(&mode) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].mode must be one of {ALLOWED_MODES:?}, got `{mode}`"
            )));
        }

        let num_samples = block
            .get("num_samples")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_SAMPLES);
        if num_samples < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].num_samples must be >= 1, got {num_samples}"
            )));
        }

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        let input_fasta = block
            .get("input_fasta")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);

        // Mode-conditional inputs.
        if (mode == "inverse-fold" || mode == "scaffold") && input_pdb.is_none() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].input_pdb required when mode = `{mode}`"
            )));
        }
        if mode == "predict" && input_fasta.is_none() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].input_fasta required when mode = `predict`"
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
                "[bio.esm3].temperature must be > 0 and finite, got {temperature}"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esm3].output_basename required"))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esm3].output_basename must be non-empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            model_variant: model_variant.to_string(),
            mode: mode.to_string(),
            num_samples,
            input_pdb,
            input_fasta,
            temperature,
            output_basename: output_basename.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_design_minimal() {
        let d = tempdir("esm3-design_min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm3.generate"

[bio.esm3]
script          = "run_esm3.py"
model_variant   = "open"
mode            = "design"
output_basename = "esm3_design"
"#,
        )
        .unwrap();
        let input = Esm3Input::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("run_esm3.py"));
        assert_eq!(input.model_variant, "open");
        assert_eq!(input.mode, "design");
        assert_eq!(input.output_basename, "esm3_design");
        // Defaults.
        assert_eq!(input.python, "python3");
        assert_eq!(input.num_samples, 4);
        assert!((input.temperature - 1.0).abs() < f64::EPSILON);
        assert!(input.input_pdb.is_none());
        assert!(input.input_fasta.is_none());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_inverse_fold_with_pdb() {
        let d = tempdir("esm3-inv_fold");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm3.generate"

[bio.esm3]
script          = "run_esm3.py"
python          = "/opt/conda/envs/esm3/bin/python"
model_variant   = "open-multimer"
mode            = "inverse-fold"
num_samples     = 16
input_pdb       = "scaffold.pdb"
temperature     = 0.7
output_basename = "esm3_invfold"
"#,
        )
        .unwrap();
        let input = Esm3Input::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/esm3/bin/python");
        assert_eq!(input.model_variant, "open-multimer");
        assert_eq!(input.mode, "inverse-fold");
        assert_eq!(input.num_samples, 16);
        assert!((input.temperature - 0.7).abs() < f64::EPSILON);
        assert_eq!(input.input_pdb, Some(PathBuf::from("scaffold.pdb")));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_inverse_fold_without_pdb() {
        let d = tempdir("esm3-inv_no_pdb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm3.generate"

[bio.esm3]
script          = "run_esm3.py"
model_variant   = "open"
mode            = "inverse-fold"
output_basename = "esm3_invfold"
"#,
        )
        .unwrap();
        let err = Esm3Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("input_pdb"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_predict_without_fasta() {
        let d = tempdir("esm3-predict_no_fasta");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm3.generate"

[bio.esm3]
script          = "run_esm3.py"
model_variant   = "open"
mode            = "predict"
output_basename = "esm3_pred"
"#,
        )
        .unwrap();
        let err = Esm3Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("input_fasta"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_mode() {
        let d = tempdir("esm3-bad_mode");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esm3.generate"

[bio.esm3]
script          = "run_esm3.py"
model_variant   = "open"
mode            = "hallucinate"
output_basename = "esm3_run"
"#,
        )
        .unwrap();
        let err = Esm3Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("mode"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
