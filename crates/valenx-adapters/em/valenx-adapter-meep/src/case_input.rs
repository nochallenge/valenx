//! `[em.meep]` case-input parsing for the Meep adapter.

use std::path::PathBuf;

use serde::Deserialize;

use valenx_core::AdapterError;

/// Parsed `[em.meep]` block.
#[derive(Clone, Debug, PartialEq)]
pub struct MeepInput {
    /// Path to the simulation script (relative to case dir). Either
    /// a Python file (preferred — `.py`) or legacy Scheme (`.ctl`).
    pub script: PathBuf,
    /// Whether `script` is Python (true) or Scheme/Meep-CTL (false).
    /// Defaults to true. Affects which interpreter prepare() picks.
    pub python: bool,
    /// Optional MPI ranks. `None` runs serial. Anything > 1 wraps
    /// the invocation in `mpirun -np N ...`.
    pub np: Option<u32>,
}

#[derive(Deserialize)]
struct CaseToml {
    case: Option<CaseHeader>,
    em: Option<EmTable>,
}

#[derive(Deserialize)]
struct CaseHeader {
    #[serde(default)]
    physics: String,
}

#[derive(Deserialize)]
struct EmTable {
    meep: Option<MeepToml>,
}

#[derive(Deserialize)]
struct MeepToml {
    script: String,
    #[serde(default = "default_python")]
    python: bool,
    #[serde(default)]
    np: Option<u32>,
}

fn default_python() -> bool {
    true
}

impl MeepInput {
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
        if let Some(ref hdr) = parsed.case {
            if !hdr.physics.is_empty()
                && !matches!(hdr.physics.as_str(), "em" | "fdtd" | "photonics" | "optics")
            {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "case physics is `{}` — Meep handles em / fdtd / photonics / optics",
                    hdr.physics
                )));
            }
        }
        let block = parsed.em.and_then(|e| e.meep).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} has no [em.meep] section — add `script = \"...\"`",
                toml_path.display()
            ))
        })?;
        if let Some(n) = block.np {
            if n == 0 {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "np must be > 0 (omit-or-1 for serial)"
                )));
            }
        }
        Ok(MeepInput {
            script: PathBuf::from(block.script),
            python: block.python,
            np: block.np,
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
    fn parses_minimal_python_script() {
        let d = tempdir("meep-python");
        write_case_toml(
            &d,
            r#"
[case]
physics = "em"

[em.meep]
script = "ring.py"
"#,
        );
        let input = MeepInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.script, PathBuf::from("ring.py"));
        assert!(input.python);
        assert!(input.np.is_none());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn picks_up_legacy_scheme_and_mpi() {
        let d = tempdir("meep-scheme-mpi");
        write_case_toml(
            &d,
            r#"
[case]
physics = "fdtd"

[em.meep]
script = "legacy.ctl"
python = false
np = 4
"#,
        );
        let input = MeepInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.script, PathBuf::from("legacy.ctl"));
        assert!(!input.python);
        assert_eq!(input.np, Some(4));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_section_actionable() {
        let d = tempdir("meep-missing");
        write_case_toml(&d, "[case]\nphysics = \"em\"\n");
        assert!(format!("{}", MeepInput::from_case_dir(&d).unwrap_err()).contains("[em.meep]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_np() {
        let d = tempdir("meep-zero-np");
        write_case_toml(&d, "[em.meep]\nscript = \"x.py\"\nnp = 0\n");
        assert!(MeepInput::from_case_dir(&d).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_wrong_physics() {
        let d = tempdir("meep-wrong");
        write_case_toml(
            &d,
            "[case]\nphysics = \"meshing\"\n[em.meep]\nscript = \"x.py\"\n",
        );
        assert!(format!("{}", MeepInput::from_case_dir(&d).unwrap_err()).contains("meshing"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
