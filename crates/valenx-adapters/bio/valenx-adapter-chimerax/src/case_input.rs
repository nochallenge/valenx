//! `[bio.chimerax]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "chimerax.script"
//!
//! [bio.chimerax]
//! script = "render.cxc"
//! nogui  = true              # optional, defaults to true
//! ```
//!
//! ChimeraX's `.cxc` command scripts are line-oriented sequences of
//! ChimeraX commands (`open 1abc; cartoon; save snapshot.png`). The
//! adapter stages the script into the workdir and invokes ChimeraX
//! with `--script`. `nogui` defaults to true so headless CI runs the
//! happy path; flip to false to record an interactive session (rare).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ChimeraXInput {
    pub script: PathBuf,
    pub nogui: bool,
}

impl ChimeraXInput {
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
            .and_then(|v| v.get("chimerax"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.chimerax] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chimerax].script required"))
            })?;
        let nogui = block.get("nogui").and_then(|v| v.as_bool()).unwrap_or(true);
        Ok(Self {
            script: PathBuf::from(script),
            nogui,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_nogui_default() {
        let d = tempdir("chimerax");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chimerax.script"

[bio.chimerax]
script = "render.cxc"
"#,
        )
        .unwrap();
        let input = ChimeraXInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("render.cxc"));
        // Default: headless. Interactive sessions are the rare path.
        assert!(input.nogui);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("chimerax");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = ChimeraXInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.chimerax]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_nogui_override() {
        // Interactive recording mode: rare in CI but a legitimate
        // path for capturing demo sessions on a workstation.
        let d = tempdir("chimerax");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chimerax.script"

[bio.chimerax]
script = "interactive.cxc"
nogui  = false
"#,
        )
        .unwrap();
        let input = ChimeraXInput::from_case_dir(&d).unwrap();
        assert!(!input.nogui);
        let _ = std::fs::remove_dir_all(&d);
    }
}
