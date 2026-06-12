//! `[bio.copasi]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "copasi.simulate"
//!
//! [bio.copasi]
//! model      = "pathway.cps"           # required, COPASI native (.cps) or SBML (.xml/.sbml)
//! report     = "report.csv"            # optional, --save target for the run report
//! run_all    = false                   # optional, if true CopasiSE executes every defined task
//! extra_args = ["--nologo"]            # optional, defaults to []
//! ```
//!
//! `model` is the COPASI input file. `CopasiSE` reads either the
//! native `.cps` archive format or an SBML `.xml`. `report = Some` adds
//! `--save <report>` so the run output lands at a known path the
//! collect step can find without walking; when omitted the adapter
//! falls back to surfacing whatever `.csv` / `.txt` files appear at the
//! workdir top-level. `run_all` is forwarded as a flag the CLI
//! recognises in its standard mode (executes every task in the file).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CopasiInput {
    pub model: PathBuf,
    pub report: Option<PathBuf>,
    pub run_all: bool,
    pub extra_args: Vec<String>,
}

impl CopasiInput {
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
            .and_then(|v| v.get("copasi"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.copasi] section",
                    case_toml.display()
                ))
            })?;

        let model_str = block.get("model").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.copasi].model required (path to .cps or SBML .xml)"
            ))
        })?;
        if model_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.copasi].model must not be empty"
            )));
        }

        let report = match block.get("report") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.copasi].report must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.copasi].report must not be empty when present"
                    )));
                }
                Some(PathBuf::from(s))
            }
            None => None,
        };

        let run_all = block
            .get("run_all")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.copasi].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.copasi].extra_args entries must be strings"
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
            report,
            run_all,
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
        let d = tempdir("copasi");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "copasi.simulate"

[bio.copasi]
model = "pathway.cps"
"#,
        )
        .unwrap();
        let input = CopasiInput::from_case_dir(&d).unwrap();
        assert_eq!(input.model, PathBuf::from("pathway.cps"));
        // Defaults: no explicit report path, single primary task,
        // no extras.
        assert!(input.report.is_none());
        assert!(!input.run_all);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_report() {
        // Pinned report path lets the collect() phase find the run
        // output deterministically (no workdir walk).
        let d = tempdir("copasi");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "copasi.simulate"

[bio.copasi]
model      = "pathway.cps"
report     = "report.csv"
run_all    = true
extra_args = ["--nologo"]
"#,
        )
        .unwrap();
        let input = CopasiInput::from_case_dir(&d).unwrap();
        assert_eq!(input.report, Some(PathBuf::from("report.csv")));
        assert!(input.run_all);
        assert_eq!(input.extra_args, vec!["--nologo".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_model() {
        let d = tempdir("copasi");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "copasi.simulate"

[bio.copasi]
model = ""
"#,
        )
        .unwrap();
        let err = CopasiInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("model"), "msg: {msg}");
        assert!(msg.contains("empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
