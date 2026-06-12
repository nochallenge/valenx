//! `[fea.openradioss]` case-input parsing for the OpenRadioss adapter.
//!
//! The adapter operates on a pre-converted engine deck (typically
//! `<root>_0001.rad`) — the starter→engine conversion phase is left
//! to the user because it's a one-time step that runs on a
//! workstation, while the engine phase is what gets queued on the
//! cluster. This split matches how most OpenRadioss workflows are
//! actually organised.

use std::path::PathBuf;

use serde::Deserialize;

use valenx_core::AdapterError;

/// Parsed `[fea.openradioss]` block from a case.toml.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenRadiossInput {
    /// Path to the engine deck (`*_0001.rad`), relative to the case
    /// directory. The adapter stages this file plus any siblings
    /// (restart files, included decks) into the workdir before run.
    pub engine_input: PathBuf,
    /// Number of MPI ranks (`-nspmd`). Default 1.
    pub nspmd: u32,
    /// Number of OpenMP threads per rank (`-nthread`). Default 1.
    pub nthread: u32,
}

#[derive(Deserialize)]
struct CaseToml {
    case: Option<CaseHeader>,
    fea: Option<FeaTable>,
}

#[derive(Deserialize)]
struct CaseHeader {
    #[serde(default)]
    physics: String,
}

#[derive(Deserialize)]
struct FeaTable {
    openradioss: Option<OpenRadiossToml>,
}

#[derive(Deserialize)]
struct OpenRadiossToml {
    engine_input: String,
    #[serde(default = "default_nspmd")]
    nspmd: u32,
    #[serde(default = "default_nthread")]
    nthread: u32,
}

fn default_nspmd() -> u32 {
    1
}
fn default_nthread() -> u32 {
    1
}

impl OpenRadiossInput {
    /// Parse `case.toml` in `case_dir` and return the OpenRadioss
    /// section. Errors if the file is missing, unparseable, the
    /// physics tag isn't FEA, or the `[fea.openradioss]` section
    /// is absent.
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let toml_path = case_dir.join("case.toml");
        // Round-18 H1 (R17 sweep gap): cap the case.toml read at the
        // shared `MAX_PROJECT_FILE_BYTES`.
        let text = valenx_core::io_caps::read_capped_to_string(
            &toml_path,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", toml_path.display())))?;
        let parsed: CaseToml = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", toml_path.display()))
        })?;
        // Physics tag check — keep error messages helpful for users
        // who picked the wrong adapter.
        if let Some(ref hdr) = parsed.case {
            if !hdr.physics.is_empty()
                && !matches!(
                    hdr.physics.as_str(),
                    "fea" | "structural" | "crash" | "impact"
                )
            {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "case physics is `{}` — OpenRadioss handles fea / structural / crash / impact",
                    hdr.physics
                )));
            }
        }
        let block = parsed.fea.and_then(|f| f.openradioss).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} has no [fea.openradioss] section — add `engine_input = \"...\"`",
                toml_path.display()
            ))
        })?;
        Ok(OpenRadiossInput {
            engine_input: PathBuf::from(block.engine_input),
            nspmd: block.nspmd,
            nthread: block.nthread,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    fn write_case_toml(dir: &std::path::Path, content: &str) {
        std::fs::write(dir.join("case.toml"), content).unwrap();
    }

    #[test]
    fn parses_minimal_engine_input_with_defaults() {
        let d = tempdir("openradioss-minimal");
        write_case_toml(
            &d,
            r#"
[case]
physics = "fea"

[fea.openradioss]
engine_input = "model_0001.rad"
"#,
        );
        let parsed = OpenRadiossInput::from_case_dir(&d).expect("parse");
        assert_eq!(parsed.engine_input, PathBuf::from("model_0001.rad"));
        assert_eq!(parsed.nspmd, 1);
        assert_eq!(parsed.nthread, 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn picks_up_explicit_parallelism() {
        let d = tempdir("openradioss-parallel");
        write_case_toml(
            &d,
            r#"
[case]
physics = "crash"

[fea.openradioss]
engine_input = "drop_test_0001.rad"
nspmd = 4
nthread = 8
"#,
        );
        let parsed = OpenRadiossInput::from_case_dir(&d).expect("parse");
        assert_eq!(parsed.engine_input, PathBuf::from("drop_test_0001.rad"));
        assert_eq!(parsed.nspmd, 4);
        assert_eq!(parsed.nthread, 8);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_section_is_actionable_error() {
        let d = tempdir("openradioss-missing");
        write_case_toml(
            &d,
            r#"
[case]
physics = "fea"
"#,
        );
        let r = OpenRadiossInput::from_case_dir(&d);
        let msg = format!("{}", r.unwrap_err());
        assert!(
            msg.contains("[fea.openradioss]"),
            "expected actionable hint, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_wrong_physics_tag() {
        let d = tempdir("openradioss-wrong-physics");
        write_case_toml(
            &d,
            r#"
[case]
physics = "cfd"

[fea.openradioss]
engine_input = "x_0001.rad"
"#,
        );
        let r = OpenRadiossInput::from_case_dir(&d);
        let msg = format!("{}", r.unwrap_err());
        assert!(msg.contains("cfd"), "expected to mention `cfd`: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_case_toml_is_actionable() {
        let d = tempdir("openradioss-no-case-toml");
        let r = OpenRadiossInput::from_case_dir(&d);
        assert!(r.is_err());
        let _ = std::fs::remove_dir_all(&d);
    }
}
