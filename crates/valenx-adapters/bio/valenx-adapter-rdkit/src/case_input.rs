//! `[bio.rdkit]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "rdkit.script"
//!
//! [bio.rdkit]
//! script  = "screen.py"
//! python  = "python3"          # optional, defaults to python3
//! smiles  = ["CCO", "c1ccccc1"]  # optional, inline SMILES list
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RdkitInput {
    pub script: PathBuf,
    pub python: String,
    /// Optional inline SMILES list. Empty vec when the case omits the
    /// key — scripts that read SMILES from files can ignore it.
    pub smiles: Vec<String>,
}

impl RdkitInput {
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
            .and_then(|v| v.get("rdkit"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.rdkit] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.rdkit].script required")))?;
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
        Ok(Self {
            script: PathBuf::from(script),
            python,
            smiles,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case() {
        let d = tempdir("rdkit");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rdkit.script"

[bio.rdkit]
script = "screen.py"
"#,
        )
        .unwrap();
        let input = RdkitInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("screen.py"));
        assert_eq!(input.python, "python3");
        assert!(input.smiles.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("rdkit");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = RdkitInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.rdkit]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn round_trips_smiles_list() {
        // RDKit users often want to pass a small set of inline
        // molecules straight from the case.toml — verify the array
        // round-trips cleanly. (Ethanol, benzene, water, caffeine.)
        let d = tempdir("rdkit");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rdkit.script"

[bio.rdkit]
script = "screen.py"
smiles = ["CCO", "c1ccccc1", "O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
"#,
        )
        .unwrap();
        let input = RdkitInput::from_case_dir(&d).unwrap();
        assert_eq!(
            input.smiles,
            vec![
                "CCO".to_string(),
                "c1ccccc1".to_string(),
                "O".to_string(),
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
