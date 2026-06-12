//! `[bio.physicell]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "physicell.simulate"
//!
//! [bio.physicell]
//! binary     = "/path/to/project"             # required, the per-project compiled binary
//! config     = "config/PhysiCell_settings.xml" # required, run-time XML configuration
//! extra_args = []                              # optional, defaults to []
//! ```
//!
//! Unlike a typical CLI tool, PhysiCell models compile to a per-project
//! C++ binary (the build is driven by the project's own `Makefile`).
//! The user therefore tells the adapter where the compiled binary
//! lives and which XML configuration to feed it. Run-time arguments
//! beyond the config path are forwarded via `extra_args`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PhysiCellInput {
    pub binary: PathBuf,
    pub config: PathBuf,
    pub extra_args: Vec<String>,
}

impl PhysiCellInput {
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
            .and_then(|v| v.get("physicell"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.physicell] section",
                    case_toml.display()
                ))
            })?;

        let binary_str = block
            .get("binary")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.physicell].binary required (path to compiled PhysiCell project binary)"
                ))
            })?;
        if binary_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.physicell].binary must not be empty"
            )));
        }

        let config_str = block
            .get("config")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.physicell].config required (path to PhysiCell_settings.xml)"
                ))
            })?;
        if config_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.physicell].config must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.physicell].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.physicell].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            binary: PathBuf::from(binary_str),
            config: PathBuf::from(config_str),
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
        let d = tempdir("physicell");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "physicell.simulate"

[bio.physicell]
binary = "./project"
config = "config/PhysiCell_settings.xml"
"#,
        )
        .unwrap();
        let input = PhysiCellInput::from_case_dir(&d).unwrap();
        assert_eq!(input.binary, PathBuf::from("./project"));
        assert_eq!(input.config, PathBuf::from("config/PhysiCell_settings.xml"));
        // Defaults: no extras.
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_binary() {
        // PhysiCell models compile per-project; an empty binary path
        // means the user hasn't pointed us at a built binary at all.
        let d = tempdir("physicell");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "physicell.simulate"

[bio.physicell]
binary = ""
config = "config/PhysiCell_settings.xml"
"#,
        )
        .unwrap();
        let err = PhysiCellInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("binary"), "msg: {msg}");
        assert!(msg.contains("empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_config() {
        let d = tempdir("physicell");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "physicell.simulate"

[bio.physicell]
binary = "./project"
config = ""
"#,
        )
        .unwrap();
        let err = PhysiCellInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("config"), "msg: {msg}");
        assert!(msg.contains("empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
