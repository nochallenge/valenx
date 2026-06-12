//! `[bio.pksim]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "pksim.simulate"
//!
//! [bio.pksim]
//! project         = "vaccine_pbpk.pksim5"
//! output_basename = "results"
//! extra_args      = []                       # optional, defaults to []
//! ```
//!
//! PK-Sim is the Open Systems Pharmacology suite's whole-body
//! physiologically based pharmacokinetic (PBPK) modeling tool. A
//! `.pksim5` project file is an XML document describing compounds,
//! formulations, individuals, populations, and the simulation
//! protocol. The headless `pksim` (a.k.a. `PKSim.CLI`) reads the
//! project in place and writes simulation results into the working
//! directory under the user-chosen basename.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PkSimInput {
    /// Path to the `.pksim5` project file. Read in place by the
    /// PK-Sim CLI (`--project <path>`); the adapter does not stage
    /// or rewrite it.
    pub project: PathBuf,
    /// Basename for output files. PK-Sim writes
    /// `<basename>*.csv` (simulation tables) and
    /// `<basename>*.json` (metadata) under this stem.
    pub output_basename: String,
    /// Additional CLI arguments appended to the `pksim` invocation.
    /// Useful for `--population <file>`, `--individuals <list>`, or
    /// solver-tuning flags exposed by the OSP CLI.
    pub extra_args: Vec<String>,
}

impl PkSimInput {
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
            .and_then(|v| v.get("pksim"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.pksim] section",
                    case_toml.display()
                ))
            })?;

        let project = block
            .get("project")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.pksim].project required (path to .pksim5 project file)"
                ))
            })?;
        if project.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pksim].project must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.pksim].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pksim].output_basename must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.pksim].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.pksim].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            project: PathBuf::from(project),
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
        let d = tempdir("pksim-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pksim.simulate"

[bio.pksim]
project         = "vaccine_pbpk.pksim5"
output_basename = "results"
"#,
        )
        .unwrap();
        let input = PkSimInput::from_case_dir(&d).unwrap();
        assert_eq!(input.project, PathBuf::from("vaccine_pbpk.pksim5"));
        assert_eq!(input.output_basename, "results");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_extra_args() {
        let d = tempdir("pksim-extras");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pksim.simulate"

[bio.pksim]
project         = "model.pksim5"
output_basename = "out"
extra_args      = ["--population", "pop.csv"]
"#,
        )
        .unwrap();
        let input = PkSimInput::from_case_dir(&d).unwrap();
        assert_eq!(input.project, PathBuf::from("model.pksim5"));
        assert_eq!(input.output_basename, "out");
        assert_eq!(
            input.extra_args,
            vec!["--population".to_string(), "pop.csv".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_output_basename() {
        // An empty output basename leaves PK-Sim with no canonical
        // place to write its results; reject up front.
        let d = tempdir("pksim-noname");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pksim.simulate"

[bio.pksim]
project         = "model.pksim5"
output_basename = ""
"#,
        )
        .unwrap();
        let err = PkSimInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
