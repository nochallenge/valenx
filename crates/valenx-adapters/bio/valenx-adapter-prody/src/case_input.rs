//! `[bio.prody]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "prody.analyze"
//!
//! [bio.prody]
//! script          = "analyse.py"
//! python          = "python3"           # optional, defaults to python3
//! input_pdb       = "1ake.pdb"
//! output_basename = "modes"
//! num_modes       = 20                  # optional, defaults to 20
//! cutoff          = 15.0                # optional, defaults to 15.0 angstroms
//! ```
//!
//! ProDy is a Python library for protein dynamics: elastic-network
//! models (ENM / GNM / ANM), normal-mode analysis, ensemble PCA, and
//! the NMD trajectory format consumed by VMD / NMWiz. The user
//! authors an `analyse.py` driver that loads `input_pdb`, builds the
//! ENM, computes the leading `num_modes` modes within the residue
//! contact `cutoff`, and writes the results into the workdir.
//!
//! `cutoff` is the contact-distance cutoff in angstroms — pairs of
//! residues within this distance get a non-zero spring in the
//! elastic-network Hessian. 15 A is the canonical default for
//! anisotropic network models on alpha-carbon coarse-grained
//! representations. Must be strictly positive and finite.
//!
//! `num_modes` is the number of low-frequency normal modes to
//! compute and surface; ProDy returns them in ascending eigenvalue
//! order. Must be >= 1.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ProdyInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute).
    pub script: PathBuf,
    /// Python interpreter to invoke. Defaults to `python3`.
    pub python: String,
    /// Path to the input `.pdb` ProDy will operate on (relative to
    /// the case directory, or absolute).
    pub input_pdb: PathBuf,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.npz` (ENM modes), `<basename>*.nmd` (NMD
    /// trajectory), and `<basename>*.csv` (analysis tables) into
    /// the workdir.
    pub output_basename: String,
    /// Number of leading low-frequency normal modes to compute.
    /// Must be >= 1.
    pub num_modes: u32,
    /// Contact-distance cutoff in angstroms for the elastic-network
    /// Hessian. Must be strictly positive and finite.
    pub cutoff: f64,
}

impl ProdyInput {
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
            .and_then(|v| v.get("prody"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.prody] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.prody].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.prody].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.prody].input_pdb required (path to .pdb file)"
                ))
            })?;
        if input_pdb.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.prody].input_pdb must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.prody].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.prody].output_basename must not be empty"
            )));
        }

        let num_modes = match block.get("num_modes") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.prody].num_modes must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.prody].num_modes must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.prody].num_modes `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 20,
        };

        let cutoff = match block.get("cutoff") {
            Some(toml::Value::Float(f)) => *f,
            Some(toml::Value::Integer(i)) => *i as f64,
            Some(_) => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.prody].cutoff must be a number"
                )));
            }
            None => 15.0,
        };
        if !cutoff.is_finite() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.prody].cutoff must be finite, got {cutoff}"
            )));
        }
        if cutoff <= 0.0 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.prody].cutoff must be > 0.0, got {cutoff}"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_pdb: PathBuf::from(input_pdb),
            output_basename: output_basename.to_string(),
            num_modes,
            cutoff,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("prody-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "prody.analyze"

[bio.prody]
script          = "analyse.py"
input_pdb       = "1ake.pdb"
output_basename = "modes"
"#,
        )
        .unwrap();
        let input = ProdyInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("analyse.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_pdb, PathBuf::from("1ake.pdb"));
        assert_eq!(input.output_basename, "modes");
        // Defaults: 20 modes, 15 A cutoff.
        assert_eq!(input.num_modes, 20);
        assert!((input.cutoff - 15.0).abs() < 1e-9);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Bigger ENM (50 modes) at a tighter contact cutoff (10 A)
        // with a pinned conda interpreter.
        let d = tempdir("prody-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "prody.analyze"

[bio.prody]
script          = "anm.py"
python          = "/opt/conda/envs/prody/bin/python"
input_pdb       = "ribosome.pdb"
output_basename = "ribo_anm"
num_modes       = 50
cutoff          = 10.0
"#,
        )
        .unwrap();
        let input = ProdyInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/prody/bin/python");
        assert_eq!(input.num_modes, 50);
        assert!((input.cutoff - 10.0).abs() < 1e-9);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_modes() {
        // Zero normal modes is meaningless — ProDy would compute
        // the Hessian and surface nothing. Reject up front so the
        // failure is fast and obvious.
        let d = tempdir("prody-zero");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "prody.analyze"

[bio.prody]
script          = "analyse.py"
input_pdb       = "1ake.pdb"
output_basename = "modes"
num_modes       = 0
"#,
        )
        .unwrap();
        let err = ProdyInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("num_modes"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_cutoff() {
        // Zero contact cutoff would zero out every spring in the
        // elastic-network Hessian — ProDy would surface a singular
        // matrix and crash on diagonalisation. Reject at validation
        // time so the user catches the typo instantly.
        let d = tempdir("prody-zerocut");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "prody.analyze"

[bio.prody]
script          = "analyse.py"
input_pdb       = "1ake.pdb"
output_basename = "modes"
cutoff          = 0.0
"#,
        )
        .unwrap();
        let err = ProdyInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("cutoff"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
