//! `[bio.plumed]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "plumed.analyze"
//!
//! [bio.plumed]
//! plumed_dat      = "plumed.dat"
//! trajectory      = "traj.xtc"
//! output_basename = "colvar"
//! kt              = 2.494                  # optional, defaults to kJ/mol at 300 K
//! extra_args      = ["--mc", "10000"]      # optional, defaults to []
//! ```
//!
//! PLUMED is the enhanced-sampling / free-energy plug-in that wraps
//! every major MD engine (GROMACS, LAMMPS, AMBER, NAMD). The
//! `plumed driver` sub-command runs PLUMED standalone over a pre-
//! computed trajectory: read frames, evaluate the collective
//! variables defined in `plumed.dat`, write COLVAR / bias files.
//!
//! `kt` is `k_B T` in the units PLUMED's input file expects (PLUMED
//! defaults to kJ/mol). 2.494 kJ/mol is room temperature (300 K)
//! and a reasonable default for protein / nucleic-acid simulations.
//! Must be strictly positive and finite — a zero or NaN `kt` would
//! crash PLUMED's reweighting on the first frame.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PlumedInput {
    /// Path to the PLUMED input file (`plumed.dat`) describing the
    /// collective variables and bias to compute. Relative paths
    /// resolve against the case directory.
    pub plumed_dat: PathBuf,
    /// Path to the input trajectory (XTC / DCD / TRR). The `plumed
    /// driver --mf_xtc` flag we wire here expects XTC; users can
    /// swap to `--mf_dcd` / `--mf_trr` via `extra_args`.
    pub trajectory: PathBuf,
    /// Filename stem used for the COLVAR / bias files PLUMED writes.
    /// Anchors `collect()`'s artefact filter so we don't surface
    /// unrelated files in the workdir.
    pub output_basename: String,
    /// `k_B T` in PLUMED's energy units (kJ/mol by default). 2.494
    /// is room temperature (300 K). Must be strictly positive and
    /// finite.
    pub kt: f64,
    /// Additional CLI arguments appended to the `plumed driver`
    /// invocation. Useful for `--mc <steps>`, `--multi <N>`, or
    /// switching the trajectory format flag.
    pub extra_args: Vec<String>,
}

impl PlumedInput {
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
            .and_then(|v| v.get("plumed"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.plumed] section",
                    case_toml.display()
                ))
            })?;

        let plumed_dat = block
            .get("plumed_dat")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.plumed].plumed_dat required"))
            })?;
        if plumed_dat.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.plumed].plumed_dat must not be empty"
            )));
        }

        let trajectory = block
            .get("trajectory")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.plumed].trajectory required"))
            })?;
        if trajectory.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.plumed].trajectory must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.plumed].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.plumed].output_basename must not be empty"
            )));
        }

        // `kt`: accept either float or integer (the user might write
        // `kt = 2` as shorthand). Reject bool / string / array.
        let kt = match block.get("kt") {
            Some(toml::Value::Float(f)) => *f,
            Some(toml::Value::Integer(i)) => *i as f64,
            Some(_) => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.plumed].kt must be a number"
                )));
            }
            None => 2.494,
        };
        if !kt.is_finite() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.plumed].kt must be finite, got {kt}"
            )));
        }
        if kt <= 0.0 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.plumed].kt must be > 0.0, got {kt}"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.plumed].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.plumed].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            plumed_dat: PathBuf::from(plumed_dat),
            trajectory: PathBuf::from(trajectory),
            output_basename: output_basename.to_string(),
            kt,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("plumed-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "plumed.analyze"

[bio.plumed]
plumed_dat      = "plumed.dat"
trajectory      = "traj.xtc"
output_basename = "colvar"
"#,
        )
        .unwrap();
        let input = PlumedInput::from_case_dir(&d).unwrap();
        assert_eq!(input.plumed_dat, PathBuf::from("plumed.dat"));
        assert_eq!(input.trajectory, PathBuf::from("traj.xtc"));
        assert_eq!(input.output_basename, "colvar");
        // Default kt = 2.494 kJ/mol (room temperature), no extras.
        assert!((input.kt - 2.494).abs() < 1e-9);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Higher-temperature reweighting (350 K ~ 2.910 kJ/mol) plus
        // a Monte Carlo step count and a switch to DCD trajectory.
        let d = tempdir("plumed-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "plumed.analyze"

[bio.plumed]
plumed_dat      = "metad.dat"
trajectory      = "metad.xtc"
output_basename = "metad_colvar"
kt              = 2.910
extra_args      = ["--mc", "100000"]
"#,
        )
        .unwrap();
        let input = PlumedInput::from_case_dir(&d).unwrap();
        assert_eq!(input.plumed_dat, PathBuf::from("metad.dat"));
        assert_eq!(input.output_basename, "metad_colvar");
        assert!((input.kt - 2.910).abs() < 1e-9);
        assert_eq!(
            input.extra_args,
            vec!["--mc".to_string(), "100000".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_kt() {
        // kt = 0 collapses every Boltzmann factor to 1 and breaks
        // PLUMED's reweighting math. Reject so the user catches the
        // typo at validation time rather than after a long PLUMED
        // spin-up.
        let d = tempdir("plumed-zerokt");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "plumed.analyze"

[bio.plumed]
plumed_dat      = "plumed.dat"
trajectory      = "traj.xtc"
output_basename = "colvar"
kt              = 0.0
"#,
        )
        .unwrap();
        let err = PlumedInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("kt"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_trajectory() {
        // `plumed driver --mf_xtc` requires a trajectory file path —
        // empty string would surface a tool-side "file not found"
        // after PLUMED has already loaded its plugins. Reject up
        // front.
        let d = tempdir("plumed-notraj");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "plumed.analyze"

[bio.plumed]
plumed_dat      = "plumed.dat"
trajectory      = ""
output_basename = "colvar"
"#,
        )
        .unwrap();
        let err = PlumedInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("trajectory"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
