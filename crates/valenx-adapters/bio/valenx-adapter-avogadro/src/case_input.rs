//! `[bio.avogadro]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "avogadro.script"
//!
//! [bio.avogadro]
//! script     = "render.py"
//! structure  = "molecule.cml"     # optional
//! headless   = true                # optional, defaults to true
//! extra_args = ["--debug"]         # optional, defaults to []
//! ```
//!
//! Avogadro 2's Python-scripted chemistry editor renders / edits
//! molecular structures via user-supplied scripts. The adapter stages
//! the script (and optional structure file) into the workdir and
//! invokes `avogadro2 --script <script>` with `--no-gui` when running
//! headlessly.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct AvogadroInput {
    pub script: PathBuf,
    pub structure: Option<PathBuf>,
    pub headless: bool,
    pub extra_args: Vec<String>,
}

impl AvogadroInput {
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
            .and_then(|v| v.get("avogadro"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.avogadro] section",
                    case_toml.display()
                ))
            })?;

        let script_str = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.avogadro].script required"))
            })?;
        if script_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.avogadro].script must not be empty"
            )));
        }

        let structure = block
            .get("structure")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);

        let headless = block
            .get("headless")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.avogadro].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.avogadro].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            script: PathBuf::from(script_str),
            structure,
            headless,
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
        let d = tempdir("avogadro");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "avogadro.script"

[bio.avogadro]
script = "render.py"
"#,
        )
        .unwrap();
        let input = AvogadroInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("render.py"));
        assert!(input.structure.is_none());
        // Default: headless. CI runs the happy path.
        assert!(input.headless);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_structure() {
        // Common case: load a CML structure, run a render script that
        // writes a PNG, exit. CML is Avogadro's native structure
        // format.
        let d = tempdir("avogadro");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "avogadro.script"

[bio.avogadro]
script     = "render.py"
structure  = "molecule.cml"
extra_args = ["--debug"]
"#,
        )
        .unwrap();
        let input = AvogadroInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("render.py"));
        assert_eq!(input.structure, Some(PathBuf::from("molecule.cml")));
        assert_eq!(input.extra_args, vec!["--debug".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn defaults_to_headless() {
        // Headless is the contract — CI / batch runs need it. The
        // explicit `headless = false` escape hatch supports the rare
        // interactive path (recording a demo session). Verify both
        // sides of that contract.
        let d = tempdir("avogadro");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "avogadro.script"

[bio.avogadro]
script   = "interactive.py"
headless = false
"#,
        )
        .unwrap();
        let input = AvogadroInput::from_case_dir(&d).unwrap();
        assert!(!input.headless);
        let _ = std::fs::remove_dir_all(&d);

        // Re-check that the absence of the key still defaults to true.
        let d2 = tempdir("avogadro");
        std::fs::write(
            d2.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "avogadro.script"

[bio.avogadro]
script = "render.py"
"#,
        )
        .unwrap();
        let input2 = AvogadroInput::from_case_dir(&d2).unwrap();
        assert!(input2.headless);
        let _ = std::fs::remove_dir_all(&d2);
    }
}
