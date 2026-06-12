//! `[bio.rfantibody]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "rfantibody.design"
//!
//! [bio.rfantibody]
//! script           = "design_rfantibody.py"
//! python           = "python3"          # optional, default python3
//! framework_pdb    = "framework.pdb"
//! target_pdb       = "antigen.pdb"
//! design_loops     = ["H3"]             # subset of CANONICAL_CDRS
//! num_designs      = 8                  # optional, default 8, must be >= 1
//! diffusion_steps  = 50                 # optional, default 50, must be >= 1
//! output_basename  = "design"
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RfAntibodyInput {
    pub script: PathBuf,
    pub python: String,
    pub framework_pdb: PathBuf,
    pub target_pdb: PathBuf,
    /// Which CDR loops to redesign. Each entry must be in
    /// [`CANONICAL_CDRS`] (case-sensitive).
    pub design_loops: Vec<String>,
    /// Number of antibody designs to sample. Defaults to 8.
    pub num_designs: u32,
    /// Number of diffusion timesteps. Defaults to 50.
    pub diffusion_steps: u32,
    /// Stem the user script should write outputs under
    /// (`{output_basename}_0.pdb`, `{output_basename}_1.pdb`, …).
    pub output_basename: String,
}

const DEFAULT_NUM_DESIGNS: u32 = 8;
const DEFAULT_DIFFUSION_STEPS: u32 = 50;

/// Recognised CDR loop names. Case-sensitive — these are the canonical
/// IMGT-style identifiers (heavy chain H1/H2/H3, light chain L1/L2/L3).
pub const CANONICAL_CDRS: &[&str] = &["H1", "H2", "H3", "L1", "L2", "L3"];

impl RfAntibodyInput {
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
            .and_then(|v| v.get("rfantibody"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.rfantibody] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rfantibody].script required"))
            })?;
        if script.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].script must be non-empty"
            )));
        }
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let framework_pdb = block
            .get("framework_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rfantibody].framework_pdb required"))
            })?;
        if framework_pdb.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].framework_pdb must be non-empty"
            )));
        }
        let target_pdb = block
            .get("target_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rfantibody].target_pdb required"))
            })?;
        if target_pdb.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].target_pdb must be non-empty"
            )));
        }
        let design_loops_raw = block
            .get("design_loops")
            .and_then(|v| v.as_array())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.rfantibody].design_loops required (array of CDR names)"
                ))
            })?;
        let mut design_loops: Vec<String> = Vec::with_capacity(design_loops_raw.len());
        for entry in design_loops_raw {
            let name = entry.as_str().ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.rfantibody].design_loops entries must be strings"
                ))
            })?;
            if !CANONICAL_CDRS.contains(&name) {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.rfantibody].design_loops entry `{name}` not recognised; \
                     must be one of {CANONICAL_CDRS:?}"
                )));
            }
            design_loops.push(name.to_string());
        }
        if design_loops.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].design_loops must be non-empty"
            )));
        }
        let num_designs = block
            .get("num_designs")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_NUM_DESIGNS);
        if num_designs < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].num_designs must be >= 1, got {num_designs}"
            )));
        }
        let diffusion_steps = block
            .get("diffusion_steps")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_DIFFUSION_STEPS);
        if diffusion_steps < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].diffusion_steps must be >= 1, got {diffusion_steps}"
            )));
        }
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rfantibody].output_basename required"))
            })?;
        if output_basename.trim().is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rfantibody].output_basename must be non-empty"
            )));
        }
        Ok(Self {
            script: PathBuf::from(script),
            python,
            framework_pdb: PathBuf::from(framework_pdb),
            target_pdb: PathBuf::from(target_pdb),
            design_loops,
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
    fn parses_minimal_with_h3() {
        let d = tempdir("rfantibody-min_h3");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfantibody.design"

[bio.rfantibody]
script          = "design.py"
framework_pdb   = "framework.pdb"
target_pdb      = "antigen.pdb"
design_loops    = ["H3"]
output_basename = "design"
"#,
        )
        .unwrap();
        let input = RfAntibodyInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.framework_pdb, PathBuf::from("framework.pdb"));
        assert_eq!(input.target_pdb, PathBuf::from("antigen.pdb"));
        assert_eq!(input.design_loops, vec!["H3".to_string()]);
        assert_eq!(input.output_basename, "design");
        // Defaults.
        assert_eq!(input.python, "python3");
        assert_eq!(input.num_designs, 8);
        assert_eq!(input.diffusion_steps, 50);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_multiple_cdrs() {
        let d = tempdir("rfantibody-multi");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfantibody.design"

[bio.rfantibody]
script          = "design.py"
python          = "/opt/conda/bin/python"
framework_pdb   = "framework.pdb"
target_pdb      = "antigen.pdb"
design_loops    = ["H1", "H2", "H3", "L3"]
num_designs     = 32
diffusion_steps = 100
output_basename = "ab_run"
"#,
        )
        .unwrap();
        let input = RfAntibodyInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/bin/python");
        assert_eq!(
            input.design_loops,
            vec![
                "H1".to_string(),
                "H2".to_string(),
                "H3".to_string(),
                "L3".to_string()
            ]
        );
        assert_eq!(input.num_designs, 32);
        assert_eq!(input.diffusion_steps, 100);
        assert_eq!(input.output_basename, "ab_run");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_invalid_cdr() {
        let d = tempdir("rfantibody-bad_cdr");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfantibody.design"

[bio.rfantibody]
script          = "design.py"
framework_pdb   = "framework.pdb"
target_pdb      = "antigen.pdb"
design_loops    = ["H3", "X9"]
output_basename = "design"
"#,
        )
        .unwrap();
        let err = RfAntibodyInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("design_loops"));
        assert!(format!("{err}").contains("X9"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_design_loops() {
        let d = tempdir("rfantibody-empty_loops");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rfantibody.design"

[bio.rfantibody]
script          = "design.py"
framework_pdb   = "framework.pdb"
target_pdb      = "antigen.pdb"
design_loops    = []
output_basename = "design"
"#,
        )
        .unwrap();
        let err = RfAntibodyInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("design_loops"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
