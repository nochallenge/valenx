//! `[bio.nwchem]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "nwchem.compute"
//!
//! [bio.nwchem]
//! input      = "input.nw"
//! output     = "output.out"
//! mpi_procs  = 8                # optional, defaults to 1 (serial)
//! extra_args = ["--bind-to", "core"]   # optional, defaults to []
//! ```
//!
//! NWChem reads its input from a single `.nw` file (its own DSL —
//! `geometry`, `basis`, `task` blocks). When `mpi_procs == 1` we
//! invoke `nwchem <input>` directly; for `mpi_procs > 1` the
//! invocation becomes `mpirun -n <N> nwchem <input>` (and the adapter
//! checks for `mpirun` on PATH at prepare time). Output is sent to
//! stdout in either case — the adapter redirects it to the named
//! output file.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct NwchemInput {
    pub input: PathBuf,
    pub output: PathBuf,
    pub mpi_procs: u32,
    pub extra_args: Vec<String>,
}

impl NwchemInput {
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
            .and_then(|v| v.get("nwchem"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.nwchem] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.nwchem].input required (path to NWChem .nw input file)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.nwchem].input must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.nwchem].output required (path for the NWChem output file)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.nwchem].output must not be empty"
            )));
        }

        let mpi_procs = match block.get("mpi_procs") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.nwchem].mpi_procs must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.nwchem].mpi_procs must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.nwchem].mpi_procs `{raw}` exceeds u32::MAX"
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
                        "[bio.nwchem].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.nwchem].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input: PathBuf::from(input_str),
            output: PathBuf::from(output_str),
            mpi_procs,
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
        let d = tempdir("nwchem");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nwchem.compute"

[bio.nwchem]
input  = "input.nw"
output = "output.out"
"#,
        )
        .unwrap();
        let input = NwchemInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.nw"));
        assert_eq!(input.output, PathBuf::from("output.out"));
        // Defaults: 1 MPI proc (serial), no extras.
        assert_eq!(input.mpi_procs, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_mpi() {
        // Parallel invocation with explicit MPI process count and an
        // mpirun-style placement extra (`--bind-to core` is the
        // OpenMPI canonical for one-process-per-core layouts).
        let d = tempdir("nwchem");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nwchem.compute"

[bio.nwchem]
input      = "h2o.nw"
output     = "h2o.out"
mpi_procs  = 16
extra_args = ["--bind-to", "core"]
"#,
        )
        .unwrap();
        let input = NwchemInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mpi_procs, 16);
        assert_eq!(
            input.extra_args,
            vec!["--bind-to".to_string(), "core".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_mpi() {
        let d = tempdir("nwchem");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nwchem.compute"

[bio.nwchem]
input     = "in.nw"
output    = "out.out"
mpi_procs = 0
"#,
        )
        .unwrap();
        let err = NwchemInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("mpi_procs"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
