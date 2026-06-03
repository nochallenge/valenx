//! `[bio.dssr]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "dssr.analyze"
//!
//! [bio.dssr]
//! input_pdb   = "structure.pdb"
//! output_json = "analysis.json"
//! extra_args  = ["--non-pair"]      # optional, defaults to []
//! ```
//!
//! DSSR reads a nucleic-acid PDB and emits a structured JSON
//! summary describing all detected structural features: base pairs,
//! multiplets, helices, stems, hairpin / internal / junction loops,
//! kissing loops, A-minor motifs, ribose zippers, pseudoknots, etc.
//! The JSON file is the canonical machine-readable output and the
//! anchor of downstream parsing.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct DssrInput {
    /// Path to the input `.pdb` file (relative to the case
    /// directory, or absolute).
    pub input_pdb: PathBuf,
    /// Path to the JSON output file DSSR will write (relative to
    /// the case directory, or absolute).
    pub output_json: PathBuf,
    /// Additional CLI arguments appended to the `x3dna-dssr`
    /// invocation.
    pub extra_args: Vec<String>,
}

impl DssrInput {
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
            .and_then(|v| v.get("dssr"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.dssr] section",
                    case_toml.display()
                ))
            })?;

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.dssr].input_pdb required")))?;
        if input_pdb.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.dssr].input_pdb must not be empty"
            )));
        }

        let output_json = block
            .get("output_json")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.dssr].output_json required"))
            })?;
        if output_json.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.dssr].output_json must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.dssr].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.dssr].extra_args entries must be strings"
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
            output_json: PathBuf::from(output_json),
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
        let d = tempdir("dssr-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dssr.analyze"

[bio.dssr]
input_pdb   = "structure.pdb"
output_json = "analysis.json"
"#,
        )
        .unwrap();
        let input = DssrInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input_pdb, PathBuf::from("structure.pdb"));
        assert_eq!(input.output_json, PathBuf::from("analysis.json"));
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input_pdb() {
        // The PDB is the entire input — empty string means DSSR
        // has no structure to work on. Reject up front.
        let d = tempdir("dssr-nopdb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dssr.analyze"

[bio.dssr]
input_pdb   = ""
output_json = "analysis.json"
"#,
        )
        .unwrap();
        let err = DssrInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input_pdb"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_output_json() {
        // Output JSON path anchors collect()'s artefact reporting;
        // empty string would leave the user with no machine-readable
        // output. Reject up front.
        let d = tempdir("dssr-nojson");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dssr.analyze"

[bio.dssr]
input_pdb   = "structure.pdb"
output_json = ""
"#,
        )
        .unwrap();
        let err = DssrInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_json"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
