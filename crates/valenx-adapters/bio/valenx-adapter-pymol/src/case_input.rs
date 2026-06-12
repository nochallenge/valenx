//! `[bio.pymol]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "pymol.script"
//!
//! [bio.pymol]
//! script     = "render.pml"
//! nogui      = true                     # optional, defaults to true
//! quiet      = true                     # optional, defaults to true
//! extra_args = ["-W", "1024", "-H", "768"]  # optional, defaults to []
//! ```
//!
//! PyMOL `.pml` scripts are line-oriented sequences of PyMOL commands
//! (`load 1abc.pdb`, `bg_color white`, `ray 1024, 768`, `png
//! out.png`). The adapter stages the script into the workdir and
//! invokes PyMOL with `-c` (headless / no-GUI) and `-q` (quiet
//! banner). Both flags default to true so headless CI takes the happy
//! path; flip either to false for the rare interactive / verbose
//! workstation run.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PymolInput {
    /// Path to the `.pml` command script (relative to the case
    /// directory, or absolute).
    pub script: PathBuf,
    /// Pass `-c` to PyMOL — runs without launching the GUI window.
    /// Defaults to true.
    pub nogui: bool,
    /// Pass `-q` to PyMOL — suppresses the startup banner. Defaults
    /// to true so log output stays focused on the script's own
    /// commands.
    pub quiet: bool,
    /// Additional CLI arguments appended after the script path. Used
    /// for `-W <width> -H <height>` and similar that the script
    /// itself can't set.
    pub extra_args: Vec<String>,
}

impl PymolInput {
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
            .and_then(|v| v.get("pymol"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.pymol] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.pymol].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pymol].script must not be empty"
            )));
        }

        let nogui = block.get("nogui").and_then(|v| v.as_bool()).unwrap_or(true);
        let quiet = block.get("quiet").and_then(|v| v.as_bool()).unwrap_or(true);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.pymol].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.pymol].extra_args entries must be strings"
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
            nogui,
            quiet,
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
        // The script-only minimal form. Both flags fall back to their
        // headless / quiet defaults; no extras.
        let d = tempdir("pymol");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pymol.script"

[bio.pymol]
script = "render.pml"
"#,
        )
        .unwrap();
        let input = PymolInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("render.pml"));
        assert!(input.nogui);
        assert!(input.quiet);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn defaults_to_nogui_quiet() {
        // No nogui / quiet keys at all — both must default true so
        // the headless / quiet path is the implicit one, matching
        // ChimeraX's `nogui = true` default.
        let d = tempdir("pymol");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pymol.script"

[bio.pymol]
script = "snapshot.pml"
"#,
        )
        .unwrap();
        let input = PymolInput::from_case_dir(&d).unwrap();
        assert!(input.nogui, "nogui must default to true");
        assert!(input.quiet, "quiet must default to true");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_explicit_overrides() {
        // Interactive workstation path: GUI on, banner on, plus
        // `-W` / `-H` extras for the window size. Pin every override
        // through.
        let d = tempdir("pymol");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pymol.script"

[bio.pymol]
script     = "interactive.pml"
nogui      = false
quiet      = false
extra_args = ["-W", "1280", "-H", "720"]
"#,
        )
        .unwrap();
        let input = PymolInput::from_case_dir(&d).unwrap();
        assert!(!input.nogui);
        assert!(!input.quiet);
        assert_eq!(
            input.extra_args,
            vec![
                "-W".to_string(),
                "1280".to_string(),
                "-H".to_string(),
                "720".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
