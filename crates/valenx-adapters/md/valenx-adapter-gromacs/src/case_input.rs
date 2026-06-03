//! `[md.gromacs]` case-input parsing for the GROMACS adapter.
//!
//! Scope today: `gmx mdrun` on a pre-built `.tpr` file. The
//! `gmx grompp` preprocessing step is left to the user — it reads
//! topology + parameter files outside the scope a single adapter
//! can sanely manage, and most workflows keep the .tpr as the
//! "ready-to-run" artifact anyway.

use std::path::PathBuf;

use serde::Deserialize;

use valenx_core::AdapterError;

/// Parsed `[md.gromacs]` block.
#[derive(Clone, Debug, PartialEq)]
pub struct GromacsInput {
    /// Path to the `.tpr` file (relative to case dir).
    pub tpr: PathBuf,
    /// Output basename for `mdrun -deffnm`. Defaults to "md" so the
    /// trajectory lands at `<workdir>/md.trr` etc.
    pub deffnm: String,
    /// Optional thread count (`mdrun -nt`). `None` lets GROMACS
    /// auto-detect.
    pub nt: Option<u32>,
}

#[derive(Deserialize)]
struct CaseToml {
    case: Option<CaseHeader>,
    md: Option<MdTable>,
}

#[derive(Deserialize)]
struct CaseHeader {
    #[serde(default)]
    physics: String,
}

#[derive(Deserialize)]
struct MdTable {
    gromacs: Option<GromacsToml>,
}

#[derive(Deserialize)]
struct GromacsToml {
    tpr: String,
    #[serde(default = "default_deffnm")]
    deffnm: String,
    #[serde(default)]
    nt: Option<u32>,
}

fn default_deffnm() -> String {
    "md".to_string()
}

impl GromacsInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let toml_path = case_dir.join("case.toml");
        // Round-18 H1 (R17 sweep gap): cap the case.toml read at the
        // shared `MAX_PROJECT_FILE_BYTES`.
        let text = valenx_core::io_caps::read_capped_to_string(
            &toml_path,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("read {}: {e}", toml_path.display()))
        })?;
        let parsed: CaseToml = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", toml_path.display()))
        })?;
        if let Some(ref hdr) = parsed.case {
            if !hdr.physics.is_empty()
                && !matches!(
                    hdr.physics.as_str(),
                    "md" | "molecular_dynamics" | "molecular-dynamics" | "biomolecular"
                )
            {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "case physics is `{}` — GROMACS handles md / biomolecular",
                    hdr.physics
                )));
            }
        }
        let block = parsed.md.and_then(|m| m.gromacs).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} has no [md.gromacs] section — add `tpr = \"...\"`",
                toml_path.display()
            ))
        })?;
        if let Some(n) = block.nt {
            if n == 0 {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "nt must be > 0 (omit for GROMACS auto-detect)"
                )));
            }
        }
        Ok(GromacsInput {
            tpr: PathBuf::from(block.tpr),
            deffnm: block.deffnm,
            nt: block.nt,
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
    fn parses_minimal_tpr() {
        let d = tempdir("gromacs-min");
        write_case_toml(
            &d,
            r#"
[case]
physics = "md"

[md.gromacs]
tpr = "system.tpr"
"#,
        );
        let input = GromacsInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.tpr, PathBuf::from("system.tpr"));
        assert_eq!(input.deffnm, "md");
        assert!(input.nt.is_none());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn picks_up_deffnm_and_nt() {
        let d = tempdir("gromacs-full");
        write_case_toml(
            &d,
            r#"
[case]
physics = "biomolecular"

[md.gromacs]
tpr = "lysozyme.tpr"
deffnm = "prod"
nt = 8
"#,
        );
        let input = GromacsInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.tpr, PathBuf::from("lysozyme.tpr"));
        assert_eq!(input.deffnm, "prod");
        assert_eq!(input.nt, Some(8));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_section_actionable() {
        let d = tempdir("gromacs-missing");
        write_case_toml(&d, "[case]\nphysics = \"md\"\n");
        assert!(
            format!("{}", GromacsInput::from_case_dir(&d).unwrap_err()).contains("[md.gromacs]")
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_nt() {
        let d = tempdir("gromacs-zero-nt");
        write_case_toml(&d, "[md.gromacs]\ntpr = \"x.tpr\"\nnt = 0\n");
        assert!(GromacsInput::from_case_dir(&d).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_wrong_physics() {
        let d = tempdir("gromacs-wrong");
        write_case_toml(
            &d,
            "[case]\nphysics = \"cfd\"\n[md.gromacs]\ntpr = \"x.tpr\"\n",
        );
        assert!(format!("{}", GromacsInput::from_case_dir(&d).unwrap_err()).contains("cfd"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
