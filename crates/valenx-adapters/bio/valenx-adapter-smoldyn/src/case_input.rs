//! `[bio.smoldyn]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "smoldyn.simulate"
//!
//! [bio.smoldyn]
//! config     = "system.txt"
//! extra_args = []                       # optional, defaults to []
//! ```
//!
//! Smoldyn is the Andrews lab's spatial stochastic reaction-diffusion
//! simulator — particles diffuse in 1D/2D/3D continuous space, react
//! according to user-defined chemistry, and bounce off geometric
//! surfaces. The whole simulation (geometry, species, reactions,
//! diffusion coefficients, output rules) is described in a single
//! plain-text configuration file (conventionally `*.txt`) that
//! Smoldyn reads as the only positional argument.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SmoldynInput {
    /// Path to the Smoldyn configuration file. Smoldyn reads it as
    /// the sole positional argument: `smoldyn <config>`. Relative
    /// paths resolve against the case directory.
    pub config: PathBuf,
    /// Additional CLI arguments appended to the smoldyn invocation.
    /// Useful for `-w` (suppress warnings), `-q` (quiet), `-t` (text
    /// output mode), or `--define name=value` overrides.
    pub extra_args: Vec<String>,
}

impl SmoldynInput {
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
            .and_then(|v| v.get("smoldyn"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.smoldyn] section",
                    case_toml.display()
                ))
            })?;

        let config = block
            .get("config")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.smoldyn].config required")))?;
        if config.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.smoldyn].config must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.smoldyn].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.smoldyn].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            config: PathBuf::from(config),
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
        let d = tempdir("smoldyn-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "smoldyn.simulate"

[bio.smoldyn]
config = "system.txt"
"#,
        )
        .unwrap();
        let input = SmoldynInput::from_case_dir(&d).unwrap();
        assert_eq!(input.config, PathBuf::from("system.txt"));
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_extra_args() {
        let d = tempdir("smoldyn-extras");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "smoldyn.simulate"

[bio.smoldyn]
config     = "diffusion.txt"
extra_args = ["-w", "-q"]
"#,
        )
        .unwrap();
        let input = SmoldynInput::from_case_dir(&d).unwrap();
        assert_eq!(input.config, PathBuf::from("diffusion.txt"));
        assert_eq!(input.extra_args, vec!["-w".to_string(), "-q".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_config() {
        // An empty config means smoldyn has no system definition —
        // it would crash immediately on startup. Reject up front so
        // the failure is fast and obvious.
        let d = tempdir("smoldyn-noconf");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "smoldyn.simulate"

[bio.smoldyn]
config = ""
"#,
        )
        .unwrap();
        let err = SmoldynInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("config"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
