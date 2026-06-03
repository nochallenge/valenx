//! `[meshing.netgen]` case-input parsing for the Netgen adapter.

use std::path::PathBuf;

use serde::Deserialize;

use valenx_core::AdapterError;

/// Parsed `[meshing.netgen]` block.
#[derive(Clone, Debug, PartialEq)]
pub struct NetgenInput {
    /// Path to the geometry source, relative to the case dir. Netgen
    /// accepts CSG (`.geo`, `.geo2d`) and BREP (`.step`, `.stp`,
    /// `.iges`, `.igs`, `.brep`) formats.
    pub geometry_file: PathBuf,
    /// Optional global mesh-size hint (passed to netgen via
    /// `-meshsize=`). `None` lets netgen use its built-in default.
    pub mesh_size: Option<f64>,
    /// Output mesh filename (relative to workdir). Default
    /// `mesh.vol`.
    pub output: String,
}

#[derive(Deserialize)]
struct CaseToml {
    case: Option<CaseHeader>,
    meshing: Option<MeshingTable>,
}

#[derive(Deserialize)]
struct CaseHeader {
    #[serde(default)]
    physics: String,
}

#[derive(Deserialize)]
struct MeshingTable {
    netgen: Option<NetgenToml>,
}

#[derive(Deserialize)]
struct NetgenToml {
    geometry_file: String,
    #[serde(default)]
    mesh_size: Option<f64>,
    #[serde(default = "default_output")]
    output: String,
}

fn default_output() -> String {
    "mesh.vol".to_string()
}

impl NetgenInput {
    /// Parse `case.toml` in `case_dir`. Errors if the file is
    /// missing, unparseable, the physics tag isn't meshing /
    /// geometry / cad, or `[meshing.netgen]` is absent.
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
                    "meshing" | "geometry" | "cad" | "mesh"
                )
            {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "case physics is `{}` — Netgen handles meshing / geometry",
                    hdr.physics
                )));
            }
        }
        let block = parsed.meshing.and_then(|m| m.netgen).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} has no [meshing.netgen] section — add `geometry_file = \"...\"`",
                toml_path.display()
            ))
        })?;
        // Sanity-check the mesh_size if provided.
        if let Some(ms) = block.mesh_size {
            if !(ms.is_finite() && ms > 0.0) {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "mesh_size must be finite and > 0; got {ms}"
                )));
            }
        }
        Ok(NetgenInput {
            geometry_file: PathBuf::from(block.geometry_file),
            mesh_size: block.mesh_size,
            output: block.output,
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
    fn parses_minimal_geometry_with_defaults() {
        let d = tempdir("netgen-min");
        write_case_toml(
            &d,
            r#"
[case]
physics = "meshing"

[meshing.netgen]
geometry_file = "shape.geo"
"#,
        );
        let input = NetgenInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.geometry_file, PathBuf::from("shape.geo"));
        assert_eq!(input.mesh_size, None);
        assert_eq!(input.output, "mesh.vol");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn picks_up_mesh_size_and_custom_output() {
        let d = tempdir("netgen-custom");
        write_case_toml(
            &d,
            r#"
[case]
physics = "meshing"

[meshing.netgen]
geometry_file = "wing.step"
mesh_size = 0.05
output = "wing_mesh.vol"
"#,
        );
        let input = NetgenInput::from_case_dir(&d).expect("parse");
        assert_eq!(input.geometry_file, PathBuf::from("wing.step"));
        assert_eq!(input.mesh_size, Some(0.05));
        assert_eq!(input.output, "wing_mesh.vol");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_negative_mesh_size() {
        let d = tempdir("netgen-neg");
        write_case_toml(
            &d,
            r#"
[meshing.netgen]
geometry_file = "x.geo"
mesh_size = -0.1
"#,
        );
        let r = NetgenInput::from_case_dir(&d);
        assert!(r.is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_section_is_actionable() {
        let d = tempdir("netgen-missing");
        write_case_toml(&d, "[case]\nphysics = \"meshing\"\n");
        let r = NetgenInput::from_case_dir(&d);
        let msg = format!("{}", r.unwrap_err());
        assert!(msg.contains("[meshing.netgen]"), "got: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_wrong_physics() {
        let d = tempdir("netgen-wrong");
        write_case_toml(
            &d,
            r#"
[case]
physics = "cfd"

[meshing.netgen]
geometry_file = "x.geo"
"#,
        );
        let r = NetgenInput::from_case_dir(&d);
        let msg = format!("{}", r.unwrap_err());
        assert!(msg.contains("cfd"), "got: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
