//! `[bio.mcell]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "mcell.simulate"
//!
//! [bio.mcell]
//! mdl        = "system.mdl"
//! # seed     = 12345                      # optional
//! extra_args = []                         # optional, defaults to []
//! ```
//!
//! MCell is the Salk Institute (Stiles, Bartol) spatial stochastic
//! cell-scale simulator — Monte Carlo diffusion of individual
//! molecules through realistic 3D subcellular geometries with
//! reactions on surfaces and in volumes. The whole simulation
//! (geometry, species, reactions, diffusion coefficients, output
//! rules, observables) is described in MDL ("Model Description
//! Language") files (conventionally `*.mdl`) that MCell reads as the
//! sole positional argument.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct McellInput {
    /// Path to the MCell MDL file. MCell reads it as the sole
    /// positional argument: `mcell <mdl>`. Relative paths resolve
    /// against the case directory.
    pub mdl: PathBuf,
    /// Optional random-number seed. When `Some(N)`, the adapter
    /// emits `-seed N` (two separate args) so MCell uses the given
    /// seed instead of its default deterministic seed (1).
    pub seed: Option<u32>,
    /// Additional CLI arguments appended to the mcell invocation.
    /// Useful for `-quiet`, `-iterations N` (override the iteration
    /// count baked into the MDL), `-checkpoint_infile`, etc.
    pub extra_args: Vec<String>,
}

impl McellInput {
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
            .and_then(|v| v.get("mcell"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.mcell] section",
                    case_toml.display()
                ))
            })?;

        let mdl = block
            .get("mdl")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.mcell].mdl required")))?;
        if mdl.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.mcell].mdl must not be empty"
            )));
        }

        let seed = match block.get("seed") {
            Some(v) => {
                let n = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mcell].seed must be a non-negative integer"
                    ))
                })?;
                if !(0..=u32::MAX as i64).contains(&n) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.mcell].seed {n} out of range for u32"
                    )));
                }
                Some(n as u32)
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.mcell].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.mcell].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            mdl: PathBuf::from(mdl),
            seed,
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
        let d = tempdir("mcell-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mcell.simulate"

[bio.mcell]
mdl = "system.mdl"
"#,
        )
        .unwrap();
        let input = McellInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mdl, PathBuf::from("system.mdl"));
        assert_eq!(input.seed, None);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_seed_and_extras() {
        let d = tempdir("mcell-seed");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mcell.simulate"

[bio.mcell]
mdl        = "diffusion.mdl"
seed       = 12345
extra_args = ["-quiet", "-iterations", "1000"]
"#,
        )
        .unwrap();
        let input = McellInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mdl, PathBuf::from("diffusion.mdl"));
        assert_eq!(input.seed, Some(12345));
        assert_eq!(
            input.extra_args,
            vec![
                "-quiet".to_string(),
                "-iterations".to_string(),
                "1000".to_string()
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_mdl() {
        // An empty MDL means MCell has no model definition — it
        // would crash immediately on startup. Reject up front so
        // the failure is fast and obvious.
        let d = tempdir("mcell-noconf");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mcell.simulate"

[bio.mcell]
mdl = ""
"#,
        )
        .unwrap();
        let err = McellInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("mdl"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
