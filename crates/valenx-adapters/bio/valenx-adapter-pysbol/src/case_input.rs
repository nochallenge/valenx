//! `[bio.pysbol]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "pysbol.compose"
//!
//! [bio.pysbol]
//! script          = "compose.py"
//! python          = "python3"             # optional, defaults to python3
//! input_sbol      = "starter.xml"         # optional — omit to compose from scratch
//! output_basename = "design"
//! ```
//!
//! pySBOL3 is the Python implementation of the Synthetic Biology
//! Open Language (SBOL) — the standard data model for genetic-
//! circuit designs. Components, sequences, interactions, and
//! constraints serialise to RDF/XML or JSON-LD and round-trip with
//! every SBOL-conformant tool (j5, Cello, SynBioHub, iBioSim, ...).
//!
//! `input_sbol` is optional: omit it to compose a brand-new SBOL
//! document from scratch, or supply an existing `.xml` that the
//! script reads as a starting point (e.g. extending a public-
//! database design with new variants).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PySbolInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute).
    pub script: PathBuf,
    /// Python interpreter to invoke. Defaults to `python3`.
    pub python: String,
    /// Optional path to an existing SBOL document the script reads
    /// as a starting point. `None` means "compose from scratch".
    pub input_sbol: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.xml` (SBOL document) and `<basename>*.json`
    /// (analysis logs) into the workdir.
    pub output_basename: String,
}

impl PySbolInput {
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
            .and_then(|v| v.get("pysbol"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.pysbol] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.pysbol].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pysbol].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_sbol = match block.get("input_sbol").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.pysbol].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pysbol].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_sbol,
            output_basename: output_basename.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("pysbol-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pysbol.compose"

[bio.pysbol]
script          = "compose.py"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = PySbolInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("compose.py"));
        assert_eq!(input.python, "python3");
        // No input_sbol — composing from scratch.
        assert_eq!(input.input_sbol, None);
        assert_eq!(input.output_basename, "design");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_sbol() {
        // Extending a public-database SBOL document with a pinned
        // conda interpreter.
        let d = tempdir("pysbol-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pysbol.compose"

[bio.pysbol]
script          = "extend.py"
python          = "/opt/conda/envs/sbol/bin/python"
input_sbol      = "starter.xml"
output_basename = "extended"
"#,
        )
        .unwrap();
        let input = PySbolInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/sbol/bin/python");
        assert_eq!(input.input_sbol, Some(PathBuf::from("starter.xml")));
        assert_eq!(input.output_basename, "extended");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename anchors collect()'s artefact filter; empty
        // string would surface every XML / JSON in the workdir,
        // including unrelated files. Reject up front.
        let d = tempdir("pysbol-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pysbol.compose"

[bio.pysbol]
script          = "compose.py"
output_basename = ""
"#,
        )
        .unwrap();
        let err = PySbolInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
