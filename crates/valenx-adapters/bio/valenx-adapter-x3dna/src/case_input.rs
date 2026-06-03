//! `[bio.x3dna]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "x3dna.analyze"
//!
//! [bio.x3dna]
//! input_pdb       = "structure.pdb"
//! output_basename = "analysis"
//! extra_args      = []           # optional, defaults to []
//! ```
//!
//! X3DNA's `analyze` driver reads a PDB containing a nucleic-acid
//! structure, identifies base pairs, and emits a base-pair / step
//! parameter table (`*.par`) plus a per-run log (`*.out`). Output
//! filenames are derived by `analyze` from the input basename;
//! `output_basename` here is surfaced into the case so `collect()`
//! can label artefacts uniformly without scraping `analyze`'s
//! filename heuristics.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct X3dnaInput {
    /// Path to the input `.pdb` file (relative to the case
    /// directory, or absolute).
    pub input_pdb: PathBuf,
    /// Filename stem the user expects X3DNA to produce. Surfaced
    /// here so `collect()` can label artefacts uniformly.
    pub output_basename: String,
    /// Additional CLI arguments appended to the `analyze` call.
    pub extra_args: Vec<String>,
}

impl X3dnaInput {
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
            .and_then(|v| v.get("x3dna"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.x3dna] section",
                    case_toml.display()
                ))
            })?;

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.x3dna].input_pdb required"))
            })?;
        if input_pdb.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.x3dna].input_pdb must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.x3dna].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.x3dna].output_basename must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.x3dna].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.x3dna].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input_pdb: PathBuf::from(input_pdb),
            output_basename: output_basename.to_string(),
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
        let d = tempdir("x3dna-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "x3dna.analyze"

[bio.x3dna]
input_pdb       = "structure.pdb"
output_basename = "analysis"
"#,
        )
        .unwrap();
        let input = X3dnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input_pdb, PathBuf::from("structure.pdb"));
        assert_eq!(input.output_basename, "analysis");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input_pdb() {
        // The PDB is the entire input — empty string means analyze
        // has no structure to work on. Reject up front.
        let d = tempdir("x3dna-nopdb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "x3dna.analyze"

[bio.x3dna]
input_pdb       = ""
output_basename = "analysis"
"#,
        )
        .unwrap();
        let err = X3dnaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input_pdb"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename anchors collect()'s artefact labels;
        // empty string would leave the user with unlabelled
        // artefacts. Reject up front.
        let d = tempdir("x3dna-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "x3dna.analyze"

[bio.x3dna]
input_pdb       = "structure.pdb"
output_basename = ""
"#,
        )
        .unwrap();
        let err = X3dnaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
