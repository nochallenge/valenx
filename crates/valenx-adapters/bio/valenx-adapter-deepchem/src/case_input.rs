//! `[bio.deepchem]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "deepchem.script"
//!
//! [bio.deepchem]
//! script      = "train.py"
//! python      = "python3"                # optional, defaults to python3
//! smiles      = ["CCO", "c1ccccc1"]      # optional, inline SMILES list
//! dataset_csv = "molecules.csv"          # optional
//! checkpoint  = "model.pt"               # optional
//! ```
//!
//! DeepChem is a PyTorch-backed cheminformatics framework — every site
//! authors its own training / inference script (the API surface is too
//! broad for a one-size-fits-all CLI). This adapter stages the
//! user-supplied Python script along with optional dataset and
//! checkpoint files, exposes the parsed knobs via `valenx_params.json`,
//! and runs `python <script>` headlessly.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct DeepChemInput {
    pub script: PathBuf,
    pub python: String,
    /// Optional inline SMILES list. Empty vec when the case omits the
    /// key — scripts that read SMILES from files or datasets ignore it.
    pub smiles: Vec<String>,
    /// Optional CSV dataset file (typical DeepChem input format).
    pub dataset_csv: Option<PathBuf>,
    /// Optional model checkpoint to resume / fine-tune from.
    pub checkpoint: Option<PathBuf>,
}

impl DeepChemInput {
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
            .and_then(|v| v.get("deepchem"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.deepchem] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.deepchem].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.deepchem].script must not be empty"
            )));
        }
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let smiles: Vec<String> = block
            .get("smiles")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_str().map(|s| s.to_string()))
                    .collect()
            })
            .unwrap_or_default();
        let dataset_csv = block
            .get("dataset_csv")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        let checkpoint = block
            .get("checkpoint")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        Ok(Self {
            script: PathBuf::from(script),
            python,
            smiles,
            dataset_csv,
            checkpoint,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("deepchem");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "deepchem.script"

[bio.deepchem]
script = "train.py"
"#,
        )
        .unwrap();
        let input = DeepChemInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("train.py"));
        assert_eq!(input.python, "python3");
        assert!(input.smiles.is_empty());
        assert!(input.dataset_csv.is_none());
        assert!(input.checkpoint.is_none());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_smiles_and_dataset() {
        // A typical DeepChem case wires together inline SMILES (small
        // bench / smoke set), a tabular dataset for the bulk training
        // signal, and an optional checkpoint to fine-tune from.
        let d = tempdir("deepchem");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "deepchem.script"

[bio.deepchem]
script      = "predict.py"
python      = "python3.11"
smiles      = ["CCO", "c1ccccc1", "O"]
dataset_csv = "molecules.csv"
checkpoint  = "best.pt"
"#,
        )
        .unwrap();
        let input = DeepChemInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("predict.py"));
        assert_eq!(input.python, "python3.11");
        assert_eq!(
            input.smiles,
            vec!["CCO".to_string(), "c1ccccc1".to_string(), "O".to_string()]
        );
        assert_eq!(input.dataset_csv, Some(PathBuf::from("molecules.csv")));
        assert_eq!(input.checkpoint, Some(PathBuf::from("best.pt")));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("deepchem");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = DeepChemInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.deepchem]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
