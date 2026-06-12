//! `[bio.simrna]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "simrna.fold"
//!
//! [bio.simrna]
//! config          = "params.dat"
//! sequence        = "rna.seq"
//! output_basename = "tertiary"
//! n_replicas      = 4                       # optional, defaults to 1
//! extra_args      = []                      # optional, defaults to []
//! ```
//!
//! SimRNA is the Bujnicki-lab Monte Carlo engine for predicting
//! three-dimensional RNA tertiary structure from primary sequence.
//! The user supplies a configuration file (`params.dat`-style) plus a
//! `.seq` file containing the RNA primary sequence; SimRNA writes
//! candidate PDB models, replica trajectories (`*.trafl`), and an
//! energy log into the working directory under the chosen basename.
//! Replica-exchange Monte Carlo is enabled by raising `n_replicas`
//! above 1.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SimRnaInput {
    /// SimRNA configuration file (`-c` flag). Read in place by the
    /// SimRNA binary; not staged.
    pub config: PathBuf,
    /// SimRNA sequence file (`-s` flag). Read in place by the
    /// SimRNA binary; not staged.
    pub sequence: PathBuf,
    /// Output basename. SimRNA writes `<basename>*.pdb`,
    /// `<basename>*.trafl` and `<basename>*.txt` under this stem.
    pub output_basename: String,
    /// Number of replica-exchange replicas (`-R` flag). 1 disables
    /// replica exchange.
    pub n_replicas: u32,
    /// Additional CLI arguments appended to the SimRNA invocation.
    pub extra_args: Vec<String>,
}

impl SimRnaInput {
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
            .and_then(|v| v.get("simrna"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.simrna] section",
                    case_toml.display()
                ))
            })?;

        let config = block
            .get("config")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.simrna].config required (path to SimRNA config file)"
                ))
            })?;
        if config.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.simrna].config must not be empty"
            )));
        }

        let sequence = block
            .get("sequence")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.simrna].sequence required (path to .seq file)"
                ))
            })?;
        if sequence.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.simrna].sequence must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.simrna].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.simrna].output_basename must not be empty"
            )));
        }

        // Replica-exchange parameter: number of replicas. 1 disables
        // replica exchange (single-trajectory MC). Reject 0 — SimRNA
        // would have nothing to integrate.
        let n_replicas = match block.get("n_replicas") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.simrna].n_replicas must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.simrna].n_replicas must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.simrna].n_replicas must fit in u32, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.simrna].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.simrna].extra_args entries must be strings"
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
            sequence: PathBuf::from(sequence),
            output_basename: output_basename.to_string(),
            n_replicas,
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
        let d = tempdir("simrna-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "simrna.fold"

[bio.simrna]
config          = "params.dat"
sequence        = "rna.seq"
output_basename = "tertiary"
"#,
        )
        .unwrap();
        let input = SimRnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.config, PathBuf::from("params.dat"));
        assert_eq!(input.sequence, PathBuf::from("rna.seq"));
        assert_eq!(input.output_basename, "tertiary");
        assert_eq!(input.n_replicas, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_replicas_and_extras() {
        let d = tempdir("simrna-replicas");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "simrna.fold"

[bio.simrna]
config          = "params.dat"
sequence        = "rna.seq"
output_basename = "tertiary"
n_replicas      = 8
extra_args      = ["-T", "300"]
"#,
        )
        .unwrap();
        let input = SimRnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.n_replicas, 8);
        assert_eq!(input.extra_args, vec!["-T".to_string(), "300".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_replicas() {
        // n_replicas = 0 leaves SimRNA with no MC trajectory to
        // integrate. Reject up front so the failure is fast and
        // obvious.
        let d = tempdir("simrna-zero");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "simrna.fold"

[bio.simrna]
config          = "params.dat"
sequence        = "rna.seq"
output_basename = "tertiary"
n_replicas      = 0
"#,
        )
        .unwrap();
        let err = SimRnaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("n_replicas"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
