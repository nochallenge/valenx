//! `[bio.alphafold3]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "alphafold3.predict"
//!
//! [bio.alphafold3]
//! run_script             = "/path/to/alphafold3/run_alphafold.py"
//! python                 = "python3"
//! input_json             = "job.json"
//! model_dir              = "/path/to/af3-weights"
//! db_dir                 = "/path/to/af3-databases"
//! num_diffusion_samples  = 5            # optional, default 5, range 1..=64
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct AlphaFold3Input {
    pub run_script: PathBuf,
    pub python: String,
    /// AF3's JSON job spec — describes the complex (protein chains,
    /// nucleic acids, ligands) AF3 will fold.
    pub input_json: PathBuf,
    /// Directory holding the AF3 model weights. The user is
    /// responsible for downloading these under the AF3 weights
    /// licence (CC-BY-NC-4.0).
    pub model_dir: PathBuf,
    /// Directory holding AF3's reference databases.
    pub db_dir: PathBuf,
    /// Number of diffusion samples per inference. AF3 defaults to 5.
    pub num_diffusion_samples: u32,
}

const DEFAULT_NUM_DIFFUSION_SAMPLES: u32 = 5;
const SAMPLES_MIN: u32 = 1;
const SAMPLES_MAX: u32 = 64;

impl AlphaFold3Input {
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
            .and_then(|v| v.get("alphafold3"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.alphafold3] section",
                    case_toml.display()
                ))
            })?;
        let run_script = block
            .get("run_script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold3].run_script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let input_json = block
            .get("input_json")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold3].input_json required"))
            })?;
        let model_dir = block
            .get("model_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold3].model_dir required"))
            })?;
        let db_dir = block
            .get("db_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold3].db_dir required"))
            })?;
        let num_diffusion_samples = block
            .get("num_diffusion_samples")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_DIFFUSION_SAMPLES);
        if !(SAMPLES_MIN..=SAMPLES_MAX).contains(&num_diffusion_samples) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.alphafold3].num_diffusion_samples must be in \
                 {SAMPLES_MIN}..={SAMPLES_MAX}, got {num_diffusion_samples}"
            )));
        }
        Ok(Self {
            run_script: PathBuf::from(run_script),
            python,
            input_json: PathBuf::from(input_json),
            model_dir: PathBuf::from(model_dir),
            db_dir: PathBuf::from(db_dir),
            num_diffusion_samples,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case() {
        let d = tempdir("alphafold3-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script = "/opt/alphafold3/run_alphafold.py"
input_json = "job.json"
model_dir  = "/opt/af3-weights"
db_dir     = "/opt/af3-db"
"#,
        )
        .unwrap();
        let input = AlphaFold3Input::from_case_dir(&d).unwrap();
        assert_eq!(
            input.run_script,
            PathBuf::from("/opt/alphafold3/run_alphafold.py")
        );
        assert_eq!(input.input_json, PathBuf::from("job.json"));
        assert_eq!(input.model_dir, PathBuf::from("/opt/af3-weights"));
        assert_eq!(input.db_dir, PathBuf::from("/opt/af3-db"));
        // Defaults.
        assert_eq!(input.python, "python3");
        assert_eq!(input.num_diffusion_samples, 5);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("alphafold3-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = AlphaFold3Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.alphafold3]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_input_json() {
        let d = tempdir("alphafold3-nojson");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script = "/opt/alphafold3/run_alphafold.py"
model_dir  = "/opt/af3-weights"
db_dir     = "/opt/af3-db"
"#,
        )
        .unwrap();
        let err = AlphaFold3Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("input_json"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_out_of_range_num_diffusion_samples() {
        let d = tempdir("alphafold3-oor");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script             = "/opt/alphafold3/run_alphafold.py"
input_json             = "job.json"
model_dir              = "/opt/af3-weights"
db_dir                 = "/opt/af3-db"
num_diffusion_samples  = 999
"#,
        )
        .unwrap();
        let err = AlphaFold3Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_diffusion_samples"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
