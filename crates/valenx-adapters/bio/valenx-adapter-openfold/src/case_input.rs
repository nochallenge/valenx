//! `[bio.openfold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "openfold.predict"
//!
//! [bio.openfold]
//! script        = "predict_openfold.py"
//! python        = "python3"          # optional, default python3
//! query_fasta   = "query.fasta"
//! model_preset  = "model_1_ptm"      # validated against the AF2 / multimer set
//! use_templates = false              # optional, default false
//! num_recycles  = 3                  # optional, default 3, range 1..=12
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct OpenFoldInput {
    pub script: PathBuf,
    pub python: String,
    pub query_fasta: PathBuf,
    /// One of the OpenFold / AlphaFold 2 model presets — validated
    /// against [`MODEL_PRESETS`] on parse.
    pub model_preset: String,
    /// Whether the OpenFold pipeline should pull in template
    /// structures during prediction. Off by default — heavy and not
    /// always available.
    pub use_templates: bool,
    /// Number of OpenFold recycling rounds. OpenFold defaults to 3.
    pub num_recycles: u32,
}

const DEFAULT_NUM_RECYCLES: u32 = 3;
const NUM_RECYCLES_MIN: u32 = 1;
const NUM_RECYCLES_MAX: u32 = 12;

/// Set of valid `model_preset` values, mirroring AlphaFold 2 + the
/// multimer v3 family that OpenFold ships pretrained weights for.
pub const MODEL_PRESETS: &[&str] = &[
    "model_1",
    "model_2",
    "model_3",
    "model_4",
    "model_5",
    "model_1_ptm",
    "model_2_ptm",
    "model_3_ptm",
    "model_4_ptm",
    "model_5_ptm",
    "model_1_multimer_v3",
    "model_2_multimer_v3",
    "model_3_multimer_v3",
    "model_4_multimer_v3",
    "model_5_multimer_v3",
];

impl OpenFoldInput {
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
            .and_then(|v| v.get("openfold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.openfold] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.openfold].script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let query_fasta = block
            .get("query_fasta")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.openfold].query_fasta required"))
            })?;
        let model_preset = block
            .get("model_preset")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.openfold].model_preset required"))
            })?;
        if !MODEL_PRESETS.contains(&model_preset) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.openfold].model_preset `{model_preset}` not recognised; \
                 must be one of {MODEL_PRESETS:?}"
            )));
        }
        let use_templates = block
            .get("use_templates")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);
        let num_recycles = block
            .get("num_recycles")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_RECYCLES);
        if !(NUM_RECYCLES_MIN..=NUM_RECYCLES_MAX).contains(&num_recycles) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.openfold].num_recycles must be in \
                 {NUM_RECYCLES_MIN}..={NUM_RECYCLES_MAX}, got {num_recycles}"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            query_fasta: PathBuf::from(query_fasta),
            model_preset: model_preset.to_string(),
            use_templates,
            num_recycles,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_defaults() {
        let d = tempdir("openfold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openfold.predict"

[bio.openfold]
script       = "predict.py"
query_fasta  = "query.fasta"
model_preset = "model_1_ptm"
"#,
        )
        .unwrap();
        let input = OpenFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("predict.py"));
        assert_eq!(input.query_fasta, PathBuf::from("query.fasta"));
        assert_eq!(input.model_preset, "model_1_ptm");
        // Defaults: python3 / no templates / 3 recycles.
        assert_eq!(input.python, "python3");
        assert!(!input.use_templates);
        assert_eq!(input.num_recycles, 3);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_template_use() {
        // Bumping `use_templates` and `num_recycles` round-trips.
        let d = tempdir("openfold-templates");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openfold.predict"

[bio.openfold]
script        = "predict.py"
query_fasta   = "query.fasta"
model_preset  = "model_5_multimer_v3"
use_templates = true
num_recycles  = 6
"#,
        )
        .unwrap();
        let input = OpenFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model_preset, "model_5_multimer_v3");
        assert!(input.use_templates);
        assert_eq!(input.num_recycles, 6);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_model_preset() {
        let d = tempdir("openfold-badpreset");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openfold.predict"

[bio.openfold]
script       = "predict.py"
query_fasta  = "q.fasta"
model_preset = "made_up_model"
"#,
        )
        .unwrap();
        let err = OpenFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("model_preset"));
        assert!(format!("{err}").contains("made_up_model"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("openfold-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = OpenFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.openfold]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
