//! `[bio.tskit]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "tskit.analyze"
//!
//! [bio.tskit]
//! script          = "analyze.py"
//! python          = "python3"          # optional, defaults to python3
//! input_trees     = "sim.trees"
//! output_basename = "stats"
//! ```
//!
//! tskit is the canonical tree-sequence analysis library — the
//! downstream consumer of the `.trees` files msprime / SLiM emit.
//! The user authors an `analyze.py` driver that loads the input
//! tree sequence, computes statistics (`pi`, `Tajima's D`, `Fst`,
//! IBD, etc.), and writes per-site / per-window outputs as TSV /
//! CSV / PNG.
//!
//! `input_trees` is required — every tskit analysis consumes a
//! pre-computed tree sequence as input. The script reads the staged
//! filename from `valenx_params.json` so it can resolve it via a
//! relative path inside the workdir (no `cd`-shenanigans needed).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct TskitInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute).
    pub script: PathBuf,
    /// Python interpreter to invoke. Defaults to `python3`.
    pub python: String,
    /// Path to the input `.trees` tree-sequence file (relative to
    /// the case directory, or absolute). Required — every tskit
    /// analysis runs on a pre-computed tree sequence.
    pub input_trees: PathBuf,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.csv` / `<basename>*.tsv` (statistics tables)
    /// and `<basename>*.png` (plots) into the workdir.
    pub output_basename: String,
}

impl TskitInput {
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
            .and_then(|v| v.get("tskit"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.tskit] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.tskit].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.tskit].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_trees = block
            .get("input_trees")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.tskit].input_trees required (path to .trees file)"
                ))
            })?;
        if input_trees.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.tskit].input_trees must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.tskit].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.tskit].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_trees: PathBuf::from(input_trees),
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
        let d = tempdir("tskit-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "tskit.analyze"

[bio.tskit]
script          = "analyze.py"
input_trees     = "sim.trees"
output_basename = "stats"
"#,
        )
        .unwrap();
        let input = TskitInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("analyze.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_trees, PathBuf::from("sim.trees"));
        assert_eq!(input.output_basename, "stats");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input_trees() {
        // Empty input_trees would surface a script-side
        // FileNotFoundError after a long Python interpreter
        // spin-up. Reject up front.
        let d = tempdir("tskit-notrees");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "tskit.analyze"

[bio.tskit]
script          = "analyze.py"
input_trees     = ""
output_basename = "stats"
"#,
        )
        .unwrap();
        let err = TskitInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input_trees"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename anchors collect()'s artefact filter;
        // empty string would surface every CSV/TSV/PNG in the
        // workdir, including unrelated files. Reject up front.
        let d = tempdir("tskit-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "tskit.analyze"

[bio.tskit]
script          = "analyze.py"
input_trees     = "sim.trees"
output_basename = ""
"#,
        )
        .unwrap();
        let err = TskitInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
