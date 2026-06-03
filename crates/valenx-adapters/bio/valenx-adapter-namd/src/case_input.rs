//! `[bio.namd]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "namd.simulate"
//!
//! [bio.namd]
//! config     = "production.namd"
//! processors = 8                       # optional, defaults to 1
//! extra_args = []                      # optional, defaults to []
//! ```
//!
//! NAMD takes a single configuration file (the `.namd` / `.conf`
//! Tcl-flavoured input deck) that pulls in the topology (`.psf`),
//! coordinates (`.pdb`), parameter files (`.prm` / `.par`), and
//! integration / output settings. Everything after the deck path on
//! the command line is forwarded to NAMD as-is via `extra_args`.
//!
//! The `processors` knob is forwarded as the Charm++ `+pN` argument
//! (a single OsString with no space — `+p4`, `+p16`, …). This is the
//! standard way to ask NAMD's bundled Charm++ runtime for shared-
//! memory parallelism. Pass `1` (the default) for the serial path.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct NamdInput {
    /// Path to the NAMD configuration file (`.namd` / `.conf`),
    /// relative to the case directory or absolute.
    pub config: PathBuf,
    /// Number of Charm++ worker threads — forwarded as `+pN` (single
    /// OsString arg, no space).
    pub processors: u32,
    /// Additional CLI arguments appended after the configuration
    /// file path.
    pub extra_args: Vec<String>,
}

impl NamdInput {
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
            .and_then(|v| v.get("namd"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.namd] section",
                    case_toml.display()
                ))
            })?;

        let config = block
            .get("config")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.namd].config required")))?;
        if config.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.namd].config must not be empty"
            )));
        }

        let processors = block
            .get("processors")
            .and_then(|v| v.as_integer())
            .map(|i| i.max(1) as u32)
            .unwrap_or(1);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.namd].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.namd].extra_args entries must be strings"
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
            processors,
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
        // Config-only form. Processors defaults to 1; extras empty.
        let d = tempdir("namd");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "namd.simulate"

[bio.namd]
config = "production.namd"
"#,
        )
        .unwrap();
        let input = NamdInput::from_case_dir(&d).unwrap();
        assert_eq!(input.config, PathBuf::from("production.namd"));
        assert_eq!(input.processors, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_processors_and_extras() {
        // Realistic shared-memory parallel run: 8 workers + a couple
        // of forwarded flags.
        let d = tempdir("namd");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "namd.simulate"

[bio.namd]
config     = "equilibrate.namd"
processors = 8
extra_args = ["+isomalloc_sync", "+setcpuaffinity"]
"#,
        )
        .unwrap();
        let input = NamdInput::from_case_dir(&d).unwrap();
        assert_eq!(input.config, PathBuf::from("equilibrate.namd"));
        assert_eq!(input.processors, 8);
        assert_eq!(
            input.extra_args,
            vec!["+isomalloc_sync".to_string(), "+setcpuaffinity".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_config_is_invalid() {
        // [bio.namd] present but no `config` — must reject loudly.
        let d = tempdir("namd");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "namd.simulate"

[bio.namd]
processors = 4
"#,
        )
        .unwrap();
        let err = NamdInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("config"),
            "error should reference missing config; got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
