//! `[bio.colabfold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "colabfold.predict"
//!
//! [bio.colabfold]
//! input_fasta   = "query.fasta"
//! num_recycles  = 3              # optional, default 3, range 1..=12
//! num_models    = 5              # optional, default 5, range 1..=5
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ColabFoldInput {
    pub input_fasta: PathBuf,
    /// Number of AlphaFold2 recycling rounds. ColabFold defaults to
    /// 3; bumping past ~6 has diminishing returns on most targets.
    pub num_recycles: u32,
    /// Number of models to ensemble. ColabFold defaults to 5; pinning
    /// to 1 is a useful shortcut for fast turnaround during sequence
    /// triage.
    pub num_models: u32,
}

const DEFAULT_NUM_RECYCLES: u32 = 3;
const DEFAULT_NUM_MODELS: u32 = 5;
const NUM_RECYCLES_MIN: u32 = 1;
const NUM_RECYCLES_MAX: u32 = 12;
const NUM_MODELS_MIN: u32 = 1;
const NUM_MODELS_MAX: u32 = 5;

impl ColabFoldInput {
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
            .and_then(|v| v.get("colabfold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.colabfold] section",
                    case_toml.display()
                ))
            })?;
        let input_fasta = block
            .get("input_fasta")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.colabfold].input_fasta required"))
            })?;
        let num_recycles = block
            .get("num_recycles")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_RECYCLES);
        if !(NUM_RECYCLES_MIN..=NUM_RECYCLES_MAX).contains(&num_recycles) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.colabfold].num_recycles must be in \
                 {NUM_RECYCLES_MIN}..={NUM_RECYCLES_MAX}, got {num_recycles}"
            )));
        }
        let num_models = block
            .get("num_models")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_MODELS);
        if !(NUM_MODELS_MIN..=NUM_MODELS_MAX).contains(&num_models) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.colabfold].num_models must be in \
                 {NUM_MODELS_MIN}..={NUM_MODELS_MAX}, got {num_models}"
            )));
        }
        Ok(Self {
            input_fasta: PathBuf::from(input_fasta),
            num_recycles,
            num_models,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_defaults() {
        let d = tempdir("colabfold");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "colabfold.predict"

[bio.colabfold]
input_fasta = "query.fasta"
"#,
        )
        .unwrap();
        let input = ColabFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input_fasta, PathBuf::from("query.fasta"));
        // Defaults match ColabFold's CLI defaults.
        assert_eq!(input.num_recycles, 3);
        assert_eq!(input.num_models, 5);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("colabfold");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = ColabFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.colabfold]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_input_fasta() {
        let d = tempdir("colabfold");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "colabfold.predict"

[bio.colabfold]
num_recycles = 5
"#,
        )
        .unwrap();
        let err = ColabFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("input_fasta"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_num_recycles_and_num_models_overrides() {
        let d = tempdir("colabfold");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "colabfold.predict"

[bio.colabfold]
input_fasta  = "target.fa"
num_recycles = 6
num_models   = 1
"#,
        )
        .unwrap();
        let input = ColabFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.num_recycles, 6);
        assert_eq!(input.num_models, 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_out_of_range_num_recycles() {
        let d = tempdir("colabfold");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "colabfold.predict"

[bio.colabfold]
input_fasta  = "x.fasta"
num_recycles = 99
"#,
        )
        .unwrap();
        let err = ColabFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_recycles"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_out_of_range_num_models() {
        // Range is 1..=5; 6 is invalid.
        let d = tempdir("colabfold");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "colabfold.predict"

[bio.colabfold]
input_fasta = "x.fasta"
num_models  = 6
"#,
        )
        .unwrap();
        let err = ColabFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_models"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
