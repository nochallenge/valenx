//! `[bio.chroma]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "chroma.design"
//!
//! [bio.chroma]
//! script          = "design_chroma.py"
//! python          = "python3"          # optional, default python3
//! num_samples     = 4                  # optional, default 4, must be >= 1
//! length          = 100                # required, must be >= 1
//! temperature     = 1.0                # optional, default 1.0, must be > 0 and finite
//! output_basename = "design"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ChromaInput {
    pub script: PathBuf,
    pub python: String,
    /// How many design samples to draw. Defaults to 4.
    pub num_samples: u32,
    /// Length (in residues) of each sampled backbone.
    pub length: u32,
    /// Sampling temperature. Defaults to 1.0; lower values produce
    /// less-diverse designs.
    pub temperature: f64,
    /// Stem the user script should write outputs under
    /// (`{output_basename}_0.pdb`, `{output_basename}_0.fa`, …).
    pub output_basename: String,
}

const DEFAULT_NUM_SAMPLES: u32 = 4;
const DEFAULT_TEMPERATURE: f64 = 1.0;

impl ChromaInput {
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
            .and_then(|v| v.get("chroma"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.chroma] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.chroma].script required")))?;
        if script.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chroma].script must be non-empty"
            )));
        }
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let num_samples = block
            .get("num_samples")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_SAMPLES);
        if num_samples < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chroma].num_samples must be >= 1, got {num_samples}"
            )));
        }
        let length = block
            .get("length")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.chroma].length required")))?;
        if length < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chroma].length must be >= 1, got {length}"
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
                "[bio.chroma].temperature must be > 0 and finite, got {temperature}"
            )));
        }
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chroma].output_basename required"))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chroma].output_basename must be non-empty"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            num_samples,
            length,
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
    fn parses_minimal() {
        let d = tempdir("chroma-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chroma.design"

[bio.chroma]
script          = "design.py"
length          = 100
output_basename = "design"
"#,
        )
        .unwrap();
        let input = ChromaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.length, 100);
        assert_eq!(input.output_basename, "design");
        // Defaults.
        assert_eq!(input.python, "python3");
        assert_eq!(input.num_samples, 4);
        assert!((input.temperature - 1.0).abs() < f64::EPSILON);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        let d = tempdir("chroma-overrides");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chroma.design"

[bio.chroma]
script          = "design.py"
python          = "/opt/conda/bin/python"
num_samples     = 16
length          = 200
temperature     = 0.5
output_basename = "chroma_run"
"#,
        )
        .unwrap();
        let input = ChromaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/bin/python");
        assert_eq!(input.num_samples, 16);
        assert_eq!(input.length, 200);
        assert!((input.temperature - 0.5).abs() < f64::EPSILON);
        assert_eq!(input.output_basename, "chroma_run");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_length() {
        let d = tempdir("chroma-zerolength");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chroma.design"

[bio.chroma]
script          = "design.py"
length          = 0
output_basename = "design"
"#,
        )
        .unwrap();
        let err = ChromaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("length"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_temperature() {
        let d = tempdir("chroma-zerotemp");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chroma.design"

[bio.chroma]
script          = "design.py"
length          = 100
temperature     = 0.0
output_basename = "design"
"#,
        )
        .unwrap();
        let err = ChromaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("temperature"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
