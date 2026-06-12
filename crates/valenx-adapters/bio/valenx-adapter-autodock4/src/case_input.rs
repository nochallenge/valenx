//! `[bio.autodock4]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "autodock4.dock"
//!
//! [bio.autodock4]
//! gpf             = "receptor.gpf"     # Grid Parameter File for autogrid4
//! dpf             = "ligand.dpf"       # Docking Parameter File for autodock4
//! skip_grid       = false              # optional, default false
//! grid_log        = "autogrid4.glg"    # optional, default "autogrid4.glg"
//! dock_log        = "autodock4.dlg"    # optional, default "autodock4.dlg"
//! extra_grid_args = []                 # optional, default []
//! extra_dock_args = []                 # optional, default []
//! ```
//!
//! AutoDock 4 splits the docking job in two stages:
//!
//! 1. `autogrid4 -p <gpf> -l <grid_log>` — pre-computes affinity
//!    grids around the receptor. Output: a set of `*.map` / `*.fld`
//!    files alongside the receptor.
//! 2. `autodock4 -p <dpf> -l <dock_log>` — runs the genetic-algorithm
//!    docking against those grids. Output: a single `.dlg` log with
//!    embedded poses + a clustered `.dlg` summary.
//!
//! `skip_grid = true` lets the user re-use existing grids from a
//! prior run (the same maps work for any number of ligands docked
//! against the same receptor); we still validate the parameter files
//! exist when the adapter actually runs.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct AutoDock4Input {
    pub gpf: PathBuf,
    pub dpf: PathBuf,
    pub skip_grid: bool,
    pub grid_log: String,
    pub dock_log: String,
    pub extra_grid_args: Vec<String>,
    pub extra_dock_args: Vec<String>,
}

impl AutoDock4Input {
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
            .and_then(|v| v.get("autodock4"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.autodock4] section",
                    case_toml.display()
                ))
            })?;

        let gpf = require_nonempty_string(block, "gpf").map(PathBuf::from)?;
        let dpf = require_nonempty_string(block, "dpf").map(PathBuf::from)?;

        let skip_grid = block
            .get("skip_grid")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let grid_log = match block.get("grid_log") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.autodock4].grid_log must be a string"
                    ))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.autodock4].grid_log must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "autogrid4.glg".to_string(),
        };

        let dock_log = match block.get("dock_log") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.autodock4].dock_log must be a string"
                    ))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.autodock4].dock_log must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "autodock4.dlg".to_string(),
        };

        let extra_grid_args = parse_optional_string_array(block, "extra_grid_args")?;
        let extra_dock_args = parse_optional_string_array(block, "extra_dock_args")?;

        Ok(Self {
            gpf,
            dpf,
            skip_grid,
            grid_log,
            dock_log,
            extra_grid_args,
            extra_dock_args,
        })
    }
}

fn require_nonempty_string(block: &toml::Value, key: &str) -> Result<String, AdapterError> {
    let s = block.get(key).and_then(|v| v.as_str()).ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!("[bio.autodock4].{key} required (string)"))
    })?;
    if s.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.autodock4].{key} must not be empty"
        )));
    }
    Ok(s.to_string())
}

fn parse_optional_string_array(
    block: &toml::Value,
    key: &str,
) -> Result<Vec<String>, AdapterError> {
    match block.get(key) {
        Some(arr) => {
            let arr = arr.as_array().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.autodock4].{key} must be an array of strings"
                ))
            })?;
            let mut out = Vec::with_capacity(arr.len());
            for entry in arr {
                let s = entry.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.autodock4].{key} entries must be strings"
                    ))
                })?;
                out.push(s.to_string());
            }
            Ok(out)
        }
        None => Ok(Vec::new()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("autodock4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "autodock4.dock"

[bio.autodock4]
gpf = "receptor.gpf"
dpf = "ligand.dpf"
"#,
        )
        .unwrap();
        let input = AutoDock4Input::from_case_dir(&d).unwrap();
        assert_eq!(input.gpf, PathBuf::from("receptor.gpf"));
        assert_eq!(input.dpf, PathBuf::from("ligand.dpf"));
        assert!(!input.skip_grid);
        assert!(input.extra_grid_args.is_empty());
        assert!(input.extra_dock_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_skip_grid() {
        // Re-using maps from a prior autogrid4 run: a common pattern
        // when docking many ligands against the same receptor pocket.
        let d = tempdir("autodock4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "autodock4.dock"

[bio.autodock4]
gpf             = "receptor.gpf"
dpf             = "ligand.dpf"
skip_grid       = true
grid_log        = "grid.log"
dock_log        = "dock.log"
extra_dock_args = ["-x"]
"#,
        )
        .unwrap();
        let input = AutoDock4Input::from_case_dir(&d).unwrap();
        assert!(input.skip_grid);
        assert_eq!(input.grid_log, "grid.log");
        assert_eq!(input.dock_log, "dock.log");
        assert_eq!(input.extra_dock_args, vec!["-x".to_string()]);
        assert!(input.extra_grid_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("autodock4");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = AutoDock4Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.autodock4]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn defaults_for_log_filenames() {
        // Default log filenames match the AutoDock User Guide (autogrid4
        // → `.glg`, autodock4 → `.dlg`); verifying defaults so users
        // upgrading don't have to re-set them in every case file.
        let d = tempdir("autodock4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "autodock4.dock"

[bio.autodock4]
gpf = "receptor.gpf"
dpf = "ligand.dpf"
"#,
        )
        .unwrap();
        let input = AutoDock4Input::from_case_dir(&d).unwrap();
        assert_eq!(input.grid_log, "autogrid4.glg");
        assert_eq!(input.dock_log, "autodock4.dlg");
        let _ = std::fs::remove_dir_all(&d);
    }
}
