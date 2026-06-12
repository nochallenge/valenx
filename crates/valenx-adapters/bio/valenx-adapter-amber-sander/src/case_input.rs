//! `[bio.sander]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "sander.simulate"
//!
//! [bio.sander]
//! topology        = "system.prmtop"
//! coordinates     = "system.inpcrd"
//! config          = "production.in"
//! output_basename = "production"
//! extra_args      = []                      # optional, defaults to []
//! ```
//!
//! sander is the OSS portion of AMBER's molecular-dynamics engine
//! (the proprietary `pmemd.cuda` ships under the AMBER license; sander
//! itself is GPL-3.0). It ingests three files at runtime: an Amber
//! topology (`-p`, `.prmtop` / `.parm7`), starting coordinates (`-c`,
//! `.inpcrd` / `.rst7`), and an mdin control deck (`-i`, `.in` /
//! `.mdin`) describing integrator, thermostat, restraints, and output
//! cadence. Outputs (`mdout`, restart `.rst`, NetCDF trajectory `.nc`,
//! and `mdinfo` checkpoint) are named via a single `output_basename`
//! prefix, so `prepare()` can wire them through `-o` / `-r` / `-x` and
//! `collect()` can match them by stem.
//!
//! All four named fields are required; `extra_args` defaults to empty.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SanderInput {
    /// Path to the Amber topology (`.prmtop` / `.parm7`). sander
    /// reads it via `-p <topology>`. Relative paths resolve against
    /// the case directory.
    pub topology: PathBuf,
    /// Path to the starting coordinates (`.inpcrd` / `.rst7`). sander
    /// reads them via `-c <coordinates>`.
    pub coordinates: PathBuf,
    /// Path to the mdin control deck (`.in` / `.mdin`). sander reads
    /// it via `-i <config>`.
    pub config: PathBuf,
    /// Common stem for sander's output files. `prepare()` wires the
    /// `<basename>.out` (mdout), `<basename>.rst` (restart), and
    /// `<basename>.nc` (NetCDF trajectory) flags from this; `collect()`
    /// also picks up `<basename>.mdinfo` checkpoints.
    pub output_basename: String,
    /// Additional CLI arguments appended to the sander invocation.
    /// Useful for `-ref <restraints>`, `-inf <mdinfo path override>`,
    /// `-AllowSmallBox`, MPI flags, etc.
    pub extra_args: Vec<String>,
}

impl SanderInput {
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
            .and_then(|v| v.get("sander"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.sander] section",
                    case_toml.display()
                ))
            })?;

        let topology = block
            .get("topology")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.sander].topology required"))
            })?;
        if topology.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.sander].topology must not be empty"
            )));
        }

        let coordinates = block
            .get("coordinates")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.sander].coordinates required"))
            })?;
        if coordinates.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.sander].coordinates must not be empty"
            )));
        }

        let config = block
            .get("config")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.sander].config required")))?;
        if config.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.sander].config must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.sander].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.sander].output_basename must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.sander].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.sander].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            topology: PathBuf::from(topology),
            coordinates: PathBuf::from(coordinates),
            config: PathBuf::from(config),
            output_basename: output_basename.to_string(),
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
        // Required-fields-only form. extras default to empty.
        let d = tempdir("amber-sander-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "sander.simulate"

[bio.sander]
topology        = "system.prmtop"
coordinates     = "system.inpcrd"
config          = "production.in"
output_basename = "production"
"#,
        )
        .unwrap();
        let input = SanderInput::from_case_dir(&d).unwrap();
        assert_eq!(input.topology, PathBuf::from("system.prmtop"));
        assert_eq!(input.coordinates, PathBuf::from("system.inpcrd"));
        assert_eq!(input.config, PathBuf::from("production.in"));
        assert_eq!(input.output_basename, "production");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_extras() {
        // Realistic re-run with restraints and a custom mdinfo path
        // forwarded verbatim to sander.
        let d = tempdir("amber-sander-extras");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "sander.simulate"

[bio.sander]
topology        = "complex.parm7"
coordinates     = "equil.rst7"
config          = "heat.mdin"
output_basename = "heat"
extra_args      = ["-ref", "equil.rst7", "-AllowSmallBox"]
"#,
        )
        .unwrap();
        let input = SanderInput::from_case_dir(&d).unwrap();
        assert_eq!(input.topology, PathBuf::from("complex.parm7"));
        assert_eq!(input.coordinates, PathBuf::from("equil.rst7"));
        assert_eq!(input.config, PathBuf::from("heat.mdin"));
        assert_eq!(input.output_basename, "heat");
        assert_eq!(
            input.extra_args,
            vec![
                "-ref".to_string(),
                "equil.rst7".to_string(),
                "-AllowSmallBox".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_output_basename() {
        // sander needs `output_basename` to wire `-o` / `-r` / `-x`;
        // if it's absent we'd silently overwrite some default file
        // name (e.g. `mdout`) and `collect()` would have nothing to
        // match. Reject up front so the failure is fast and obvious.
        let d = tempdir("amber-sander-nobasename");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "sander.simulate"

[bio.sander]
topology    = "system.prmtop"
coordinates = "system.inpcrd"
config      = "production.in"
"#,
        )
        .unwrap();
        let err = SanderInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("output_basename"),
            "error should reference missing output_basename; got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
