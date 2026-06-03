//! `[bio.scvi]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "scvi.train"
//!
//! [bio.scvi]
//! script       = "train.py"
//! python       = "python3"          # optional, default python3
//! input_h5ad   = "raw.h5ad"
//! output_h5ad  = "with_latent.h5ad"
//! model        = "scvi"             # "scvi" | "scanvi" | "totalvi" | "linear-scvi"
//! n_latent     = 10                 # optional, default 10
//! n_hidden     = 128                # optional, default 128
//! n_layers     = 2                  # optional, default 2
//! max_epochs   = 400                # optional, default 400
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ScviInput {
    pub script: PathBuf,
    pub python: String,
    pub input_h5ad: PathBuf,
    pub output_h5ad: String,
    /// Which scvi-tools model variant to train. Pinned to the small
    /// set of stable single-cell models the upstream library exposes.
    pub model: String,
    /// Latent-space dimension for the variational autoencoder.
    /// scvi-tools defaults to 10; bumping to 30 helps complex tissues.
    pub n_latent: u32,
    /// Hidden-layer width for both encoder and decoder MLPs.
    pub n_hidden: u32,
    /// Number of hidden layers in encoder / decoder.
    pub n_layers: u32,
    /// Maximum training epochs. scvi-tools defaults to 400 with an
    /// early-stopping callback that typically halts well before.
    pub max_epochs: u32,
}

const DEFAULT_MODEL: &str = "scvi";
const DEFAULT_N_LATENT: u32 = 10;
const DEFAULT_N_HIDDEN: u32 = 128;
const DEFAULT_N_LAYERS: u32 = 2;
const DEFAULT_MAX_EPOCHS: u32 = 400;

/// The set of model strings the user can specify under `[bio.scvi]`.
/// Each maps to an scvi-tools class:
/// - `scvi` → `scvi.model.SCVI`
/// - `scanvi` → `scvi.model.SCANVI` (semi-supervised, label-aware)
/// - `totalvi` → `scvi.model.TOTALVI` (CITE-seq joint RNA + protein)
/// - `linear-scvi` → `scvi.model.LinearSCVI` (linear decoder for
///   interpretable latent factors)
const ALLOWED_MODELS: &[&str] = &["scvi", "scanvi", "totalvi", "linear-scvi"];

impl ScviInput {
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
            .and_then(|v| v.get("scvi"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.scvi] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .map(|s| s.trim().to_string())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.scvi].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].script must be non-empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_h5ad = block
            .get("input_h5ad")
            .and_then(|v| v.as_str())
            .map(|s| s.trim().to_string())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.scvi].input_h5ad required"))
            })?;
        if input_h5ad.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].input_h5ad must be non-empty"
            )));
        }

        let output_h5ad = block
            .get("output_h5ad")
            .and_then(|v| v.as_str())
            .map(|s| s.trim().to_string())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.scvi].output_h5ad required"))
            })?;
        if output_h5ad.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].output_h5ad must be non-empty"
            )));
        }

        let model = block
            .get("model")
            .and_then(|v| v.as_str())
            .unwrap_or(DEFAULT_MODEL)
            .to_string();
        if !ALLOWED_MODELS.contains(&model.as_str()) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].model must be one of {ALLOWED_MODELS:?}, got `{model}`"
            )));
        }

        let n_latent = block
            .get("n_latent")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_N_LATENT);
        if n_latent < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].n_latent must be >= 1, got {n_latent}"
            )));
        }

        let n_hidden = block
            .get("n_hidden")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_N_HIDDEN);
        if n_hidden < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].n_hidden must be >= 1, got {n_hidden}"
            )));
        }

        let n_layers = block
            .get("n_layers")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_N_LAYERS);
        if n_layers < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].n_layers must be >= 1, got {n_layers}"
            )));
        }

        let max_epochs = block
            .get("max_epochs")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_MAX_EPOCHS);
        if max_epochs < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scvi].max_epochs must be >= 1, got {max_epochs}"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_h5ad: PathBuf::from(input_h5ad),
            output_h5ad,
            model,
            n_latent,
            n_hidden,
            n_layers,
            max_epochs,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("scvi-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scvi.train"

[bio.scvi]
script      = "train.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "with_latent.h5ad"
"#,
        )
        .unwrap();
        let input = ScviInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("train.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_h5ad, PathBuf::from("raw.h5ad"));
        assert_eq!(input.output_h5ad, "with_latent.h5ad");
        // Defaults pinned to scvi-tools' canonical hyperparameters.
        assert_eq!(input.model, "scvi");
        assert_eq!(input.n_latent, 10);
        assert_eq!(input.n_hidden, 128);
        assert_eq!(input.n_layers, 2);
        assert_eq!(input.max_epochs, 400);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_scanvi_model() {
        // SCANVI is the semi-supervised label-aware variant — common
        // when the user has a partially-annotated reference dataset.
        // Round-trip the override cleanly.
        let d = tempdir("scvi-scanvi");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scvi.train"

[bio.scvi]
script      = "train.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "annotated.h5ad"
model       = "scanvi"
n_latent    = 30
n_hidden    = 256
n_layers    = 3
max_epochs  = 200
"#,
        )
        .unwrap();
        let input = ScviInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model, "scanvi");
        assert_eq!(input.n_latent, 30);
        assert_eq!(input.n_hidden, 256);
        assert_eq!(input.n_layers, 3);
        assert_eq!(input.max_epochs, 200);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_model() {
        let d = tempdir("scvi-badmodel");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scvi.train"

[bio.scvi]
script      = "train.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "out.h5ad"
model       = "destvi"
"#,
        )
        .unwrap();
        let err = ScviInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("model"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_latent() {
        let d = tempdir("scvi-zerolat");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scvi.train"

[bio.scvi]
script      = "train.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "out.h5ad"
n_latent    = 0
"#,
        )
        .unwrap();
        let err = ScviInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("n_latent"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
