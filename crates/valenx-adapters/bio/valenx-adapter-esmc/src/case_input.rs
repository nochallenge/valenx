//! `[bio.esmc]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "esmc.embed"
//!
//! [bio.esmc]
//! script          = "embed_esmc.py"
//! python          = "python3"          # optional, default python3
//! input_fasta     = "query.fasta"
//! model_variant   = "esmc-300m"        # one of "esmc-300m" | "esmc-600m"
//! pooling         = "per-residue"      # one of "per-residue" | "mean"
//! output_basename = "embeddings"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct EsmcInput {
    pub script: PathBuf,
    pub python: String,
    pub input_fasta: PathBuf,
    pub model_variant: String,
    pub pooling: String,
    pub output_basename: String,
}

const ALLOWED_VARIANTS: &[&str] = &["esmc-300m", "esmc-600m"];
const ALLOWED_POOLING: &[&str] = &["per-residue", "mean"];

impl EsmcInput {
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
            .and_then(|v| v.get("esmc"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.esmc] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esmc].script required")))?;
        if script.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esmc].script must be non-empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_fasta = block
            .get("input_fasta")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esmc].input_fasta required"))
            })?;
        if input_fasta.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esmc].input_fasta must be non-empty"
            )));
        }

        let model_variant = block
            .get("model_variant")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esmc].model_variant required"))
            })?;
        if !ALLOWED_VARIANTS.contains(&model_variant) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esmc].model_variant must be one of {ALLOWED_VARIANTS:?}, got `{model_variant}`"
            )));
        }

        let pooling = block
            .get("pooling")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esmc].pooling required")))?;
        if !ALLOWED_POOLING.contains(&pooling) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esmc].pooling must be one of {ALLOWED_POOLING:?}, got `{pooling}`"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esmc].output_basename required"))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esmc].output_basename must be non-empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_fasta: PathBuf::from(input_fasta),
            model_variant: model_variant.to_string(),
            pooling: pooling.to_string(),
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
        let d = tempdir("esmc-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmc.embed"

[bio.esmc]
script          = "embed.py"
input_fasta     = "query.fasta"
model_variant   = "esmc-300m"
pooling         = "per-residue"
output_basename = "embeddings"
"#,
        )
        .unwrap();
        let input = EsmcInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("embed.py"));
        assert_eq!(input.input_fasta, PathBuf::from("query.fasta"));
        assert_eq!(input.model_variant, "esmc-300m");
        assert_eq!(input.pooling, "per-residue");
        assert_eq!(input.output_basename, "embeddings");
        assert_eq!(input.python, "python3");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_600m_and_mean_pooling() {
        let d = tempdir("esmc-600m");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmc.embed"

[bio.esmc]
script          = "embed.py"
python          = "/opt/conda/envs/esmc/bin/python"
input_fasta     = "many.fasta"
model_variant   = "esmc-600m"
pooling         = "mean"
output_basename = "esmc_600_mean"
"#,
        )
        .unwrap();
        let input = EsmcInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/esmc/bin/python");
        assert_eq!(input.model_variant, "esmc-600m");
        assert_eq!(input.pooling, "mean");
        assert_eq!(input.output_basename, "esmc_600_mean");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_variant() {
        let d = tempdir("esmc-bad_variant");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmc.embed"

[bio.esmc]
script          = "embed.py"
input_fasta     = "query.fasta"
model_variant   = "esmc-3b"
pooling         = "mean"
output_basename = "embeddings"
"#,
        )
        .unwrap();
        let err = EsmcInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("model_variant"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_pooling() {
        let d = tempdir("esmc-bad_pool");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmc.embed"

[bio.esmc]
script          = "embed.py"
input_fasta     = "query.fasta"
model_variant   = "esmc-300m"
pooling         = "max"
output_basename = "embeddings"
"#,
        )
        .unwrap();
        let err = EsmcInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("pooling"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
