//! `[bio.rfdiffusion]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "rfdiffusion.design"
//!
//! [bio.rfdiffusion]
//! script           = "design_rfdiffusion.py"
//! python           = "python3"          # optional, default python3
//! input_pdb        = "scaffold.pdb"
//! mode             = "motif"            # motif | binder | unconditional | partial-diffusion
//! num_designs      = 8                  # optional, default 8, must be >= 1
//! diffusion_steps  = 50                 # optional, default 50, must be >= 1
//! output_basename  = "design"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RfDiffusionInput {
    pub script: PathBuf,
    pub python: String,
    pub input_pdb: PathBuf,
    /// One of [`MODES`] — design strategy passed through to the user
    /// script via `valenx_params.json`.
    pub mode: String,
    /// Number of designs to sample. Defaults to 8.
    pub num_designs: u32,
    /// Number of diffusion timesteps. Defaults to 50.
    pub diffusion_steps: u32,
    /// Stem the user script should write outputs under
    /// (`{output_basename}_0.pdb`, `{output_basename}_1.pdb`, …).
    pub output_basename: String,
}

const DEFAULT_NUM_DESIGNS: u32 = 8;
const DEFAULT_DIFFUSION_STEPS: u32 = 50;

/// Recognised `mode` values. `partial-diffusion` is RFdiffusion's
/// "denoise an existing backbone" workflow; the other three are the
/// standard de novo strategies.
pub const MODES: &[&str] = &["motif", "binder", "unconditional", "partial-diffusion"];

impl RfDiffusionInput {
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
            .and_then(|v| v.get("rfdiffusion"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.rfdiffusion] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rfdiffusion].script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rfdiffusion].input_pdb required"))
            })?;
        let mode = block.get("mode").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!("[bio.rfdiffusion].mode required"))
        })?;
        if !MODES.contains(&mode) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfdiffusion].mode `{mode}` not recognised; must be one of {MODES:?}"
            )));
        }
        let num_designs = block
            .get("num_designs")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_DESIGNS);
        if num_designs < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfdiffusion].num_designs must be >= 1, got {num_designs}"
            )));
        }
        let diffusion_steps = block
            .get("diffusion_steps")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_DIFFUSION_STEPS);
        if diffusion_steps < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfdiffusion].diffusion_steps must be >= 1, got {diffusion_steps}"
            )));
        }
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.rfdiffusion].output_basename required"
                ))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfdiffusion].output_basename must be non-empty"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_pdb: PathBuf::from(input_pdb),
            mode: mode.to_string(),
            num_designs,
            diffusion_steps,
            output_basename: output_basename.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("rfdiffusion-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfdiffusion.design"

[bio.rfdiffusion]
script          = "design.py"
input_pdb       = "scaffold.pdb"
mode            = "motif"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = RfDiffusionInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.input_pdb, PathBuf::from("scaffold.pdb"));
        assert_eq!(input.mode, "motif");
        assert_eq!(input.output_basename, "design");
        // Defaults.
        assert_eq!(input.python, "python3");
        assert_eq!(input.num_designs, 8);
        assert_eq!(input.diffusion_steps, 50);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_binder_mode() {
        let d = tempdir("rfdiffusion-binder");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfdiffusion.design"

[bio.rfdiffusion]
script          = "design.py"
python          = "/opt/conda/bin/python"
input_pdb       = "target.pdb"
mode            = "binder"
num_designs     = 32
diffusion_steps = 100
output_basename = "binder_run"
"#,
        )
        .unwrap();
        let input = RfDiffusionInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mode, "binder");
        assert_eq!(input.num_designs, 32);
        assert_eq!(input.diffusion_steps, 100);
        assert_eq!(input.python, "/opt/conda/bin/python");
        assert_eq!(input.output_basename, "binder_run");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_mode() {
        let d = tempdir("rfdiffusion-badmode");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfdiffusion.design"

[bio.rfdiffusion]
script          = "design.py"
input_pdb       = "scaffold.pdb"
mode            = "made_up_mode"
output_basename = "design"
"#,
        )
        .unwrap();
        let err = RfDiffusionInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("mode"));
        assert!(format!("{err}").contains("made_up_mode"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_designs() {
        let d = tempdir("rfdiffusion-zerod");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfdiffusion.design"

[bio.rfdiffusion]
script          = "design.py"
input_pdb       = "scaffold.pdb"
mode            = "motif"
num_designs     = 0
output_basename = "design"
"#,
        )
        .unwrap();
        let err = RfDiffusionInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("num_designs"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
