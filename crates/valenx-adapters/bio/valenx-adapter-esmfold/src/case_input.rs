//! `[bio.esmfold]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "esmfold.predict"
//!
//! [bio.esmfold]
//! script        = "predict_esmfold.py"
//! python        = "python3"        # optional, default python3
//! query_fasta   = "query.fasta"
//! output_pdb    = "prediction.pdb"
//! num_recycles  = 4                # optional, default 4, range 1..=12
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct EsmFoldInput {
    pub script: PathBuf,
    pub python: String,
    pub query_fasta: PathBuf,
    pub output_pdb: PathBuf,
    /// Number of ESMFold recycling rounds. ESMFold defaults to 4;
    /// bumping past ~8 has diminishing returns on most targets.
    pub num_recycles: u32,
}

const DEFAULT_NUM_RECYCLES: u32 = 4;
const NUM_RECYCLES_MIN: u32 = 1;
const NUM_RECYCLES_MAX: u32 = 12;

impl EsmFoldInput {
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
            .and_then(|v| v.get("esmfold"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.esmfold] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.esmfold].script required")))?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let query_fasta = block
            .get("query_fasta")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esmfold].query_fasta required"))
            })?;
        let output_pdb = block
            .get("output_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.esmfold].output_pdb required"))
            })?;
        let num_recycles = block
            .get("num_recycles")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_RECYCLES);
        if !(NUM_RECYCLES_MIN..=NUM_RECYCLES_MAX).contains(&num_recycles) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.esmfold].num_recycles must be in \
                 {NUM_RECYCLES_MIN}..={NUM_RECYCLES_MAX}, got {num_recycles}"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            query_fasta: PathBuf::from(query_fasta),
            output_pdb: PathBuf::from(output_pdb),
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
        let d = tempdir("esmfold-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmfold.predict"

[bio.esmfold]
script      = "predict.py"
query_fasta = "query.fasta"
output_pdb  = "prediction.pdb"
"#,
        )
        .unwrap();
        let input = EsmFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("predict.py"));
        assert_eq!(input.query_fasta, PathBuf::from("query.fasta"));
        assert_eq!(input.output_pdb, PathBuf::from("prediction.pdb"));
        // Default python interpreter, default recycle count.
        assert_eq!(input.python, "python3");
        assert_eq!(input.num_recycles, 4);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("esmfold-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = EsmFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.esmfold]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_python_override() {
        // ESMFold is typically pinned to a conda env distinct from
        // the system Python — the user supplies the env's interpreter
        // path. Round-trip cleanly through the case-input parser.
        let d = tempdir("esmfold-py");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmfold.predict"

[bio.esmfold]
script      = "predict.py"
python      = "/opt/conda/envs/esmfold/bin/python"
query_fasta = "query.fasta"
output_pdb  = "out.pdb"
"#,
        )
        .unwrap();
        let input = EsmFoldInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/esmfold/bin/python");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_out_of_range_num_recycles() {
        let d = tempdir("esmfold-oor");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "esmfold.predict"

[bio.esmfold]
script       = "predict.py"
query_fasta  = "q.fasta"
output_pdb   = "out.pdb"
num_recycles = 99
"#,
        )
        .unwrap();
        let err = EsmFoldInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_recycles"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
