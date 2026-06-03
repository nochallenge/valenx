//! `[bio.bionetgen]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "bionetgen.simulate"
//!
//! [bio.bionetgen]
//! model           = "egfr.bngl"     # required, BNGL rule-based model
//! output_basename = "egfr"          # required, prefix for .net/.gdat/.cdat outputs
//! generate_only   = false           # optional, if true skip the simulate block (network only)
//! extra_args      = []              # optional, defaults to []
//! ```
//!
//! BioNetGen's canonical entry point is `BNG2.pl`, a Perl driver that
//! reads a `.bngl` model and runs the actions inside (`generate_network`,
//! `simulate`, `simulate_ssa`, etc.). `output_basename` becomes the
//! `-o` prefix every output file inherits — pin it so collect() finds
//! `<basename>.net`, `<basename>.gdat`, `<basename>.cdat` deterministically.
//! `generate_only = true` adds `--no-execute`, which skips simulation
//! actions and emits just the expanded reaction network (`.net`).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct BioNetGenInput {
    pub model: PathBuf,
    pub output_basename: String,
    pub generate_only: bool,
    pub extra_args: Vec<String>,
}

impl BioNetGenInput {
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
            .and_then(|v| v.get("bionetgen"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.bionetgen] section",
                    case_toml.display()
                ))
            })?;

        let model_str = block.get("model").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.bionetgen].model required (path to .bngl model)"
            ))
        })?;
        if model_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bionetgen].model must not be empty"
            )));
        }

        let output_basename_str = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.bionetgen].output_basename required (prefix for .net/.gdat/.cdat outputs)"
                ))
            })?;
        let output_basename = output_basename_str.trim().to_string();
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.bionetgen].output_basename must not be empty (after trim)"
            )));
        }

        let generate_only = block
            .get("generate_only")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.bionetgen].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.bionetgen].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            model: PathBuf::from(model_str),
            output_basename,
            generate_only,
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
        let d = tempdir("bionetgen");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bionetgen.simulate"

[bio.bionetgen]
model           = "egfr.bngl"
output_basename = "egfr"
"#,
        )
        .unwrap();
        let input = BioNetGenInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model, PathBuf::from("egfr.bngl"));
        assert_eq!(input.output_basename, "egfr");
        // Defaults: full simulate (not network-only), no extras.
        assert!(!input.generate_only);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_generate_only() {
        // Network-only mode: skip simulate actions, emit just the
        // expanded reaction network as `<basename>.net`.
        let d = tempdir("bionetgen");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bionetgen.simulate"

[bio.bionetgen]
model           = "fceri.bngl"
output_basename = "fceri"
generate_only   = true
extra_args      = ["--log"]
"#,
        )
        .unwrap();
        let input = BioNetGenInput::from_case_dir(&d).unwrap();
        assert!(input.generate_only);
        assert_eq!(input.extra_args, vec!["--log".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Whitespace-only counts as empty after trim — pin the check
        // so a typo doesn't silently produce files named ".net" etc.
        let d = tempdir("bionetgen");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bionetgen.simulate"

[bio.bionetgen]
model           = "egfr.bngl"
output_basename = "   "
"#,
        )
        .unwrap();
        let err = BioNetGenInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        assert!(msg.contains("empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
