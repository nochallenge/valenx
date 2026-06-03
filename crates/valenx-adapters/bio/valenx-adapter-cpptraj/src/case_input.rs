//! `[bio.cpptraj]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "cpptraj.analyze"
//!
//! [bio.cpptraj]
//! script     = "analyse.in"
//! topology   = "system.prmtop"
//! extra_args = ["-tl"]                # optional, defaults to []
//! ```
//!
//! cpptraj is AmberTools' canonical trajectory analysis tool —
//! consumes Amber `.prmtop` / `.parm7` topologies plus
//! `.nc` / `.dcd` / `.mdcrd` trajectories, runs an analysis script
//! authored in cpptraj's domain language (`trajin ...`,
//! `rms ...`, `radgyr ...`, `hbond ...`), and writes results into
//! the workdir as `.dat` (one line per frame), `.agr` (XmGrace
//! plot data), or `.gnu` (gnuplot script + data).
//!
//! The `script` parameter points at the analysis script the user
//! wrote (typically a `.in` or `.cpptraj` file); `topology` is the
//! Amber topology cpptraj reads via `-p`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CpptrajInput {
    /// Path to the cpptraj analysis script. cpptraj reads it via
    /// `-i <script>`. Relative paths resolve against the case
    /// directory.
    pub script: PathBuf,
    /// Path to the Amber topology (`.prmtop` / `.parm7`). cpptraj
    /// reads it via `-p <topology>`.
    pub topology: PathBuf,
    /// Additional CLI arguments appended to the cpptraj invocation.
    /// Useful for `-tl` (toggle long output), `--mpi`, or extra
    /// `-y <traj>` topology / trajectory pairs the script references.
    pub extra_args: Vec<String>,
}

impl CpptrajInput {
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
            .and_then(|v| v.get("cpptraj"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.cpptraj] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.cpptraj].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cpptraj].script must not be empty"
            )));
        }

        let topology = block
            .get("topology")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cpptraj].topology required"))
            })?;
        if topology.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cpptraj].topology must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.cpptraj].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.cpptraj].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            script: PathBuf::from(script),
            topology: PathBuf::from(topology),
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
        let d = tempdir("cpptraj-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cpptraj.analyze"

[bio.cpptraj]
script   = "analyse.in"
topology = "system.prmtop"
"#,
        )
        .unwrap();
        let input = CpptrajInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("analyse.in"));
        assert_eq!(input.topology, PathBuf::from("system.prmtop"));
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_script() {
        // Empty script means cpptraj has no analysis instructions —
        // it would launch and exit immediately with no useful work.
        // Reject up front so the failure is fast and obvious.
        let d = tempdir("cpptraj-noscript");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cpptraj.analyze"

[bio.cpptraj]
script   = ""
topology = "system.prmtop"
"#,
        )
        .unwrap();
        let err = CpptrajInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("script"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_topology() {
        // cpptraj's `-p <topology>` is mandatory — empty string would
        // crash on startup (no atom definitions, no residue table).
        // Reject at validation time.
        let d = tempdir("cpptraj-notop");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cpptraj.analyze"

[bio.cpptraj]
script   = "analyse.in"
topology = ""
"#,
        )
        .unwrap();
        let err = CpptrajInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("topology"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
