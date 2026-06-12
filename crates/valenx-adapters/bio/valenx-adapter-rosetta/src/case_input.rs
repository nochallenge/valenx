//! `[bio.rosetta]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "rosetta.protocol"
//!
//! [bio.rosetta]
//! protocol        = "design.xml"
//! input_pdb       = "input.pdb"
//! output_basename = "design"
//! nstruct         = 10
//! database        = "/opt/rosetta/main/database"
//! extra_args      = ["-out:overwrite"]      # optional, defaults to []
//! ```
//!
//! Rosetta's `rosetta_scripts` binary takes an XML protocol that
//! drives the modeling pipeline (filters, movers, scorefunctions),
//! plus an input `.pdb` to operate on. `nstruct` is the number of
//! independent decoys to produce; the binary writes
//! `<output_basename>_<NNNN>.pdb` files plus a single `score.sc`
//! tab-separated scorefile summarising the run.
//!
//! `database` is required — every `rosetta_scripts` invocation needs
//! `-database <path>` pointing at the Rosetta data directory (energy
//! tables, fragment libraries, etc.). It's bundled with the source
//! distribution but isn't on PATH; the user must surface its absolute
//! path here.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RosettaInput {
    /// Path to the XML protocol script (relative to the case
    /// directory, or absolute).
    pub protocol: PathBuf,
    /// Path to the input `.pdb` Rosetta will operate on.
    pub input_pdb: PathBuf,
    /// Stem the binary uses to label output decoys
    /// (`<basename>_0001.pdb`, `<basename>_0002.pdb`, ...).
    pub output_basename: String,
    /// Number of independent decoys to generate. Rosetta runs
    /// `nstruct` independent trajectories; must be >= 1.
    pub nstruct: u32,
    /// Path to the Rosetta data `database/` directory. Required —
    /// `rosetta_scripts` won't start without `-database <path>`.
    pub database: PathBuf,
    /// Additional CLI arguments appended to the invocation.
    pub extra_args: Vec<String>,
}

impl RosettaInput {
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
            .and_then(|v| v.get("rosetta"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.rosetta] section",
                    case_toml.display()
                ))
            })?;

        let protocol = block
            .get("protocol")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rosetta].protocol required"))
            })?;
        if protocol.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rosetta].protocol must not be empty"
            )));
        }

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rosetta].input_pdb required"))
            })?;
        if input_pdb.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rosetta].input_pdb must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.rosetta].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rosetta].output_basename must not be empty"
            )));
        }

        let nstruct = match block.get("nstruct") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.rosetta].nstruct must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rosetta].nstruct must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.rosetta].nstruct `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.rosetta].nstruct required (number of decoys, >= 1)"
                )));
            }
        };

        let database = block
            .get("database")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.rosetta].database required (path to Rosetta data dir)"
                ))
            })?;
        if database.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.rosetta].database must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.rosetta].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.rosetta].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            protocol: PathBuf::from(protocol),
            input_pdb: PathBuf::from(input_pdb),
            output_basename: output_basename.to_string(),
            nstruct,
            database: PathBuf::from(database),
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
        let d = tempdir("rosetta-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosetta.protocol"

[bio.rosetta]
protocol        = "design.xml"
input_pdb       = "input.pdb"
output_basename = "design"
nstruct         = 10
database        = "/opt/rosetta/main/database"
"#,
        )
        .unwrap();
        let input = RosettaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.protocol, PathBuf::from("design.xml"));
        assert_eq!(input.input_pdb, PathBuf::from("input.pdb"));
        assert_eq!(input.output_basename, "design");
        assert_eq!(input.nstruct, 10);
        assert_eq!(input.database, PathBuf::from("/opt/rosetta/main/database"));
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Larger decoy run with an extras list — `-out:overwrite`
        // is the canonical "don't fail if a previous run already
        // wrote some decoys" flag from the Rosetta options manual.
        let d = tempdir("rosetta-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosetta.protocol"

[bio.rosetta]
protocol        = "fastdesign.xml"
input_pdb       = "scaffold.pdb"
output_basename = "fastdesign_run"
nstruct         = 1000
database        = "/usr/local/rosetta/database"
extra_args      = ["-out:overwrite", "-multithreading:total_threads", "8"]
"#,
        )
        .unwrap();
        let input = RosettaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.nstruct, 1000);
        assert_eq!(input.protocol, PathBuf::from("fastdesign.xml"));
        assert_eq!(input.output_basename, "fastdesign_run");
        assert_eq!(
            input.extra_args,
            vec![
                "-out:overwrite".to_string(),
                "-multithreading:total_threads".to_string(),
                "8".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_protocol() {
        // The XML protocol drives the entire pipeline — empty
        // string means rosetta_scripts has no work to do. Reject up
        // front so the user catches the typo at validation time.
        let d = tempdir("rosetta-noprot");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosetta.protocol"

[bio.rosetta]
protocol        = ""
input_pdb       = "input.pdb"
output_basename = "design"
nstruct         = 1
database        = "/opt/rosetta/main/database"
"#,
        )
        .unwrap();
        let err = RosettaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("protocol"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_nstruct() {
        // A zero-decoy run is meaningless — Rosetta would start,
        // load the protocol, and exit with no decoys produced.
        // Reject here so the failure is fast and obvious.
        let d = tempdir("rosetta-zero");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosetta.protocol"

[bio.rosetta]
protocol        = "design.xml"
input_pdb       = "input.pdb"
output_basename = "design"
nstruct         = 0
database        = "/opt/rosetta/main/database"
"#,
        )
        .unwrap();
        let err = RosettaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("nstruct"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_database() {
        // `-database <path>` is mandatory for rosetta_scripts —
        // empty string would crash on startup. Reject up front.
        let d = tempdir("rosetta-nodb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosetta.protocol"

[bio.rosetta]
protocol        = "design.xml"
input_pdb       = "input.pdb"
output_basename = "design"
nstruct         = 5
database        = ""
"#,
        )
        .unwrap();
        let err = RosettaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("database"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
