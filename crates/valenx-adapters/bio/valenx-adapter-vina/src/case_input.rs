//! `[bio.vina]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "vina.dock"
//!
//! [bio.vina]
//! receptor       = "receptor.pdbqt"
//! ligand         = "ligand.pdbqt"
//! output         = "out.pdbqt"
//! center         = [10.0, 12.0, 18.0]   # x, y, z (Å) — search-box centre
//! size           = [20.0, 20.0, 20.0]   # x, y, z (Å) — search-box edges
//! exhaustiveness = 8                    # optional, default 8 (1..=32)
//! num_modes      = 9                    # optional, default 9
//! energy_range   = 3.0                  # optional, default 3.0 (kcal/mol)
//! cpu            = 0                    # optional, default 0 (= auto-detect)
//! extra_args     = []                   # optional, default []
//! ```
//!
//! AutoDock Vina is the de facto open-source flexible-ligand docking
//! engine: it places a ligand (PDBQT) into a rigid receptor (PDBQT)
//! within a user-defined search box and reports the lowest-energy
//! poses ranked by predicted binding affinity. The grid centre +
//! grid edge lengths fully describe the search volume; the three
//! tuning knobs trade speed for thoroughness:
//!
//! * `exhaustiveness` — Monte-Carlo restart count. Vina's CLI clamps
//!   this to 1..=32 with a warning beyond that band, so we reject
//!   anything outside the supported range up front.
//! * `num_modes` — maximum number of poses to write. Must be ≥ 1
//!   (zero would emit an empty output file).
//! * `energy_range` — kcal/mol cutoff above the best pose for which
//!   modes are still reported. Vina rejects non-finite values, so we
//!   require strictly positive and finite.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct VinaInput {
    pub receptor: PathBuf,
    pub ligand: PathBuf,
    pub output: PathBuf,
    pub center: [f64; 3],
    pub size: [f64; 3],
    pub exhaustiveness: u32,
    pub num_modes: u32,
    pub energy_range: f64,
    pub cpu: u32,
    pub extra_args: Vec<String>,
    /// Backing engine: `"native"` (default, uses valenx-dock) or
    /// `"external"` (subprocess to the `vina` binary).
    pub engine: String,
}

impl VinaInput {
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
            .and_then(|v| v.get("vina"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.vina] section",
                    case_toml.display()
                ))
            })?;

        let receptor = require_nonempty_path(block, "receptor")?;
        let ligand = require_nonempty_path(block, "ligand")?;
        let output = require_nonempty_path(block, "output")?;

        let center = require_f64_triplet(block, "center")?;
        let size = require_f64_triplet(block, "size")?;
        // `require_f64_triplet` already rejected NaN / infinity, so a
        // simple `<= 0.0` check is sufficient here.
        for (axis, v) in ["x", "y", "z"].iter().zip(size.iter()) {
            if *v <= 0.0 {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.vina].size[{axis}] must be > 0.0, got {v}"
                )));
            }
        }

        let exhaustiveness = match block.get("exhaustiveness") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].exhaustiveness must be an integer"
                    ))
                })?;
                if !(1..=32).contains(&raw) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].exhaustiveness must be in 1..=32, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 8,
        };

        let num_modes = match block.get("num_modes") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.vina].num_modes must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].num_modes must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 9,
        };

        let energy_range = match block.get("energy_range") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.vina].energy_range must be a number"
                        ))
                    })?;
                if !raw.is_finite() || raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].energy_range must be finite and > 0.0, got {raw}"
                    )));
                }
                raw
            }
            None => 3.0,
        };

        let cpu = match block.get("cpu") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.vina].cpu must be an integer"))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].cpu must be >= 0 (0 = auto), got {raw}"
                    )));
                }
                raw as u32
            }
            None => 0,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.vina].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        let engine = match block.get("engine") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.vina].engine must be a string"))
                })?;
                if s != "native" && s != "external" {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.vina].engine must be \"native\" or \"external\", got `{s}`"
                    )));
                }
                s.to_string()
            }
            None => "native".to_string(),
        };

        Ok(Self {
            receptor,
            ligand,
            output,
            center,
            size,
            exhaustiveness,
            num_modes,
            energy_range,
            cpu,
            extra_args,
            engine,
        })
    }
}

fn require_nonempty_path(block: &toml::Value, key: &str) -> Result<PathBuf, AdapterError> {
    let s = block.get(key).and_then(|v| v.as_str()).ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!("[bio.vina].{key} required (string)"))
    })?;
    if s.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.vina].{key} must not be empty"
        )));
    }
    Ok(PathBuf::from(s))
}

fn require_f64_triplet(block: &toml::Value, key: &str) -> Result<[f64; 3], AdapterError> {
    let arr = block.get(key).and_then(|v| v.as_array()).ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!(
            "[bio.vina].{key} required (array of three numbers)"
        ))
    })?;
    if arr.len() != 3 {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.vina].{key} must have exactly 3 elements, got {}",
            arr.len()
        )));
    }
    let mut out = [0.0_f64; 3];
    for (i, v) in arr.iter().enumerate() {
        let raw = v
            .as_float()
            .or_else(|| v.as_integer().map(|i| i as f64))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.vina].{key}[{i}] must be a number"))
            })?;
        if !raw.is_finite() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.vina].{key}[{i}] must be finite, got {raw}"
            )));
        }
        out[i] = raw;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor = "rec.pdbqt"
ligand   = "lig.pdbqt"
output   = "out.pdbqt"
center   = [0.0, 0.0, 0.0]
size     = [20.0, 20.0, 20.0]
"#,
        )
        .unwrap();
        let input = VinaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.receptor, PathBuf::from("rec.pdbqt"));
        assert_eq!(input.ligand, PathBuf::from("lig.pdbqt"));
        assert_eq!(input.output, PathBuf::from("out.pdbqt"));
        assert_eq!(input.center, [0.0, 0.0, 0.0]);
        assert_eq!(input.size, [20.0, 20.0, 20.0]);
        // Defaults: exhaustiveness=8, num_modes=9, energy_range=3.0,
        // cpu=0 (auto), no extra args.
        assert_eq!(input.exhaustiveness, 8);
        assert_eq!(input.num_modes, 9);
        assert!((input.energy_range - 3.0).abs() < 1e-12);
        assert_eq!(input.cpu, 0);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor       = "rec.pdbqt"
ligand         = "lig.pdbqt"
output         = "poses.pdbqt"
center         = [12.5, -3.0, 8.25]
size           = [25.0, 25.0, 30.0]
exhaustiveness = 16
num_modes      = 20
energy_range   = 5.5
cpu            = 4
extra_args     = ["--seed", "12345"]
"#,
        )
        .unwrap();
        let input = VinaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.center, [12.5, -3.0, 8.25]);
        assert_eq!(input.size, [25.0, 25.0, 30.0]);
        assert_eq!(input.exhaustiveness, 16);
        assert_eq!(input.num_modes, 20);
        assert!((input.energy_range - 5.5).abs() < 1e-12);
        assert_eq!(input.cpu, 4);
        assert_eq!(
            input.extra_args,
            vec!["--seed".to_string(), "12345".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_size() {
        // A search box with a zero-length axis is geometrically degenerate.
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor = "rec.pdbqt"
ligand   = "lig.pdbqt"
output   = "out.pdbqt"
center   = [0.0, 0.0, 0.0]
size     = [20.0, 0.0, 20.0]
"#,
        )
        .unwrap();
        let err = VinaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("size[y]"), "msg: {msg}");
        assert!(msg.contains("> 0.0"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_exhaustiveness_above_32() {
        // Vina's CLI caps useful exhaustiveness around 32; beyond that
        // the runtime cost balloons without measurable accuracy gain.
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor       = "rec.pdbqt"
ligand         = "lig.pdbqt"
output         = "out.pdbqt"
center         = [0.0, 0.0, 0.0]
size           = [20.0, 20.0, 20.0]
exhaustiveness = 64
"#,
        )
        .unwrap();
        let err = VinaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("1..=32"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_modes() {
        // num_modes=0 would emit an empty output file — never useful.
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor  = "rec.pdbqt"
ligand    = "lig.pdbqt"
output    = "out.pdbqt"
center    = [0.0, 0.0, 0.0]
size      = [20.0, 20.0, 20.0]
num_modes = 0
"#,
        )
        .unwrap();
        let err = VinaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_modes must be >= 1"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn defaults_to_native_engine() {
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor = "rec.pdbqt"
ligand   = "lig.pdbqt"
output   = "out.pdbqt"
center   = [0.0, 0.0, 0.0]
size     = [20.0, 20.0, 20.0]
"#,
        )
        .unwrap();
        let input = VinaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.engine, "native");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn accepts_external_engine() {
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor = "rec.pdbqt"
ligand   = "lig.pdbqt"
output   = "out.pdbqt"
center   = [0.0, 0.0, 0.0]
size     = [20.0, 20.0, 20.0]
engine   = "external"
"#,
        )
        .unwrap();
        let input = VinaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.engine, "external");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_engine() {
        let d = tempdir("vina");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor = "rec.pdbqt"
ligand   = "lig.pdbqt"
output   = "out.pdbqt"
center   = [0.0, 0.0, 0.0]
size     = [20.0, 20.0, 20.0]
engine   = "tensorflow"
"#,
        )
        .unwrap();
        let err = VinaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("engine"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
