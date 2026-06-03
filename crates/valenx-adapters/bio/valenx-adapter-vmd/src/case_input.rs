//! `[bio.vmd]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "vmd.script"
//!
//! [bio.vmd]
//! script     = "trajectory.tcl"
//! headless   = true                # optional, defaults to true
//! structure  = "topology.psf"      # optional — if present, passed before -e
//! extra_args = ["-args", "frame=0"]  # optional, defaults to []
//! ```
//!
//! VMD `.tcl` scripts drive trajectory loading, representation
//! styling, frame iteration, and image rendering through
//! `render`/`render TachyonInternal` calls. The adapter stages the
//! script (and optional structure file) into the workdir and
//! invokes VMD with `-dispdev text -e <script>`.
//!
//! `-dispdev text` selects the headless renderer (no GUI window /
//! OpenGL context). Defaults to true so headless CI runs the happy
//! path; flip `headless = false` when capturing an interactive
//! session on a workstation. `structure` is the optional topology
//! file that VMD will load before executing the script (typical for
//! MD setups: load the `.psf` so the script can `mol addfile` the
//! trajectory).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct VmdInput {
    /// Path to the `.tcl` driver script (relative to the case
    /// directory, or absolute).
    pub script: PathBuf,
    /// Pass `-dispdev text` to VMD — runs without an OpenGL context
    /// or display. Defaults to true.
    pub headless: bool,
    /// Optional topology / structure file VMD loads before executing
    /// the script. Useful for the canonical MD setup of "load the
    /// topology, then `mol addfile` the trajectory in the script."
    pub structure: Option<PathBuf>,
    /// Additional CLI arguments appended after the script path.
    pub extra_args: Vec<String>,
}

impl VmdInput {
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
            .and_then(|v| v.get("vmd"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.vmd] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.vmd].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.vmd].script must not be empty"
            )));
        }

        let headless = block
            .get("headless")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        let structure = block
            .get("structure")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.vmd].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.vmd].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            script: PathBuf::from(script),
            headless,
            structure,
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
        // Script-only form. Headless defaults to true; no structure
        // file; no extras.
        let d = tempdir("vmd");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vmd.script"

[bio.vmd]
script = "trajectory.tcl"
"#,
        )
        .unwrap();
        let input = VmdInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("trajectory.tcl"));
        assert!(input.headless);
        assert_eq!(input.structure, None);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_structure() {
        // Canonical MD setup: load the topology PSF before the
        // script so `mol addfile` works in the script body.
        let d = tempdir("vmd");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vmd.script"

[bio.vmd]
script    = "render_traj.tcl"
structure = "topology.psf"
"#,
        )
        .unwrap();
        let input = VmdInput::from_case_dir(&d).unwrap();
        assert_eq!(input.structure, Some(PathBuf::from("topology.psf")));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn defaults_to_headless() {
        // No headless key — must default to true so CI takes the
        // `-dispdev text` path.
        let d = tempdir("vmd");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vmd.script"

[bio.vmd]
script = "any.tcl"
"#,
        )
        .unwrap();
        let input = VmdInput::from_case_dir(&d).unwrap();
        assert!(input.headless, "headless must default to true");
        let _ = std::fs::remove_dir_all(&d);
    }
}
