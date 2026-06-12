//! `[fea.code_aster]` case-input parsing for the Code_Aster adapter.
//!
//! Code_Aster's run model is driven by a `.export` file that
//! describes the comm + mmed + version + memory limits all in one
//! place. The user provides the .export and any companion files
//! it references; the adapter stages the lot into the workdir and
//! invokes `as_run`.

use std::path::PathBuf;

use serde::Deserialize;

use valenx_core::AdapterError;

/// Parsed `[fea.code_aster]` block.
#[derive(Clone, Debug, PartialEq)]
pub struct CodeAsterInput {
    /// Path to the `.export` file (relative to case dir).
    pub export: PathBuf,
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
    code_aster: Option<CodeAsterToml>,
}

#[derive(Deserialize)]
struct CodeAsterToml {
    export: String,
}

impl CodeAsterInput {
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
                && !matches!(
                    hdr.physics.as_str(),
                    "fea" | "structural" | "thermomechanical" | "multi-physics" | "multiphysics"
                )
            {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "case physics is `{}` — Code_Aster handles fea / thermomechanical / multi-physics",
                    hdr.physics
                )));
            }
        }
        let block = parsed.fea.and_then(|f| f.code_aster).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} has no [fea.code_aster] section — add `export = \"...\"`",
                toml_path.display()
            ))
        })?;
        Ok(CodeAsterInput {
            export: PathBuf::from(block.export),
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
    fn parses_minimal_export() {
        let d = tempdir("code-aster-min");
        write_case_toml(
            &d,
            r#"
[case]
physics = "fea"

[fea.code_aster]
export = "case.export"
"#,
        );
        let input = CodeAsterInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.export, PathBuf::from("case.export"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_section_actionable() {
        let d = tempdir("code-aster-missing");
        write_case_toml(&d, "[case]\nphysics = \"fea\"\n");
        assert!(
            format!("{}", CodeAsterInput::from_case_dir(&d).unwrap_err())
                .contains("[fea.code_aster]")
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_wrong_physics() {
        let d = tempdir("code-aster-wrong");
        write_case_toml(
            &d,
            "[case]\nphysics = \"cfd\"\n[fea.code_aster]\nexport = \"x.export\"\n",
        );
        assert!(format!("{}", CodeAsterInput::from_case_dir(&d).unwrap_err()).contains("cfd"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn accepts_thermomechanical_and_multiphysics() {
        for physics in ["thermomechanical", "multi-physics", "multiphysics"] {
            let d = tempdir("code-aster-variants");
            write_case_toml(
                &d,
                &format!(
                    "[case]\nphysics = \"{physics}\"\n[fea.code_aster]\nexport = \"x.export\"\n"
                ),
            );
            assert!(
                CodeAsterInput::from_case_dir(&d).is_ok(),
                "should accept physics={physics}"
            );
            let _ = std::fs::remove_dir_all(&d);
        }
    }
}
