//! Phase 121 — extended Wavefront OBJ reader (materials + groups).
//!
//! ## What OCCT does
//!
//! `RWObj_Reader` parses an OBJ file and, when an `mtllib` directive
//! is present, fetches the named `.mtl` from the same directory.
//! Materials are exposed via the `XCAFDoc_VisMaterialTool` so
//! downstream renderers can bind them to face groups. Per-group
//! `usemtl` switches are recorded as `XCAFDoc_DocumentTool`
//! attribute entries.
//!
//! ## v1 status
//!
//! **Honest implementation** for the mesh geometry only —
//! delegates to [`valenx_mesh::format::obj::read_path`]. The
//! sidecar `.mtl` is parsed into a [`crate::obj_writer_extended::MaterialLib`]
//! when present, and the bundle is returned alongside the mesh.
//! Per-group material binding is deferred to Phase 121.5 (needs
//! a `Mesh` element-block material-index field).

use std::path::Path;

use valenx_mesh::Mesh;

use crate::error::OcctExchangeError;
use crate::obj_writer_extended::{MaterialLib, ObjMaterial};

/// Geometry + sidecar materials recovered from a `.obj` file.
#[derive(Clone, Debug)]
pub struct ObjImport {
    /// Triangle-surface mesh.
    pub mesh: Mesh,
    /// Materials from the sidecar `.mtl`, if any. Empty otherwise.
    pub materials: MaterialLib,
}

/// Read a `.obj` file plus its sidecar `.mtl` (if any).
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.obj`.
/// - [`OcctExchangeError::Parse`] for malformed OBJ.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn obj_reader_extended(path: &Path) -> Result<ObjImport, OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("obj") {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .obj",
        ));
    }
    let mesh = valenx_mesh::format::obj::read_path(path)
        .map_err(|e| OcctExchangeError::parse("obj file", format!("{e}")))?;

    // Scan the OBJ for `mtllib <name>.mtl` to know what sidecar to
    // load. v1 only honours the first directive.
    //
    // Round-23 sweep: bound both reads at MAX_OBJ_FILE_BYTES (1 GiB)
    // — OBJ + MTL for a million-vertex mesh in ASCII is in the low
    // hundreds of MiB; 1 GiB refuses the cat /dev/zero DoS.
    let obj_text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_OBJ_FILE_BYTES as usize,
    )?;
    let mut materials = MaterialLib::default();
    if let Some(mtl_name) = obj_text
        .lines()
        .find_map(|l| l.trim().strip_prefix("mtllib "))
    {
        let mtl_path = path.with_file_name(mtl_name.trim());
        if mtl_path.exists() {
            let mtl_text = valenx_core::io_caps::read_capped_to_string(
                &mtl_path,
                valenx_core::io_caps::MAX_OBJ_FILE_BYTES as usize,
            )?;
            materials = parse_mtl(&mtl_text);
        }
    }
    Ok(ObjImport { mesh, materials })
}

/// Pure parser for `.mtl` text — pulled out so tests don't hit
/// disk.
fn parse_mtl(text: &str) -> MaterialLib {
    let mut lib = MaterialLib::default();
    let mut current: Option<ObjMaterial> = None;
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.split_whitespace();
        let Some(tag) = parts.next() else { continue };
        match tag {
            "newmtl" => {
                if let Some(mat) = current.take() {
                    lib.materials.push(mat);
                }
                let name = parts.collect::<Vec<_>>().join(" ");
                current = Some(ObjMaterial {
                    name,
                    diffuse: [1.0, 1.0, 1.0],
                    ambient: [0.0, 0.0, 0.0],
                    specular: [0.0, 0.0, 0.0],
                    map_kd: None,
                });
            }
            "Kd" | "Ka" | "Ks" => {
                if let Some(mat) = current.as_mut() {
                    let rgb = parse_rgb(parts);
                    match tag {
                        "Kd" => mat.diffuse = rgb,
                        "Ka" => mat.ambient = rgb,
                        "Ks" => mat.specular = rgb,
                        _ => unreachable!(),
                    }
                }
            }
            "map_Kd" => {
                if let Some(mat) = current.as_mut() {
                    let s = parts.collect::<Vec<_>>().join(" ");
                    if !s.is_empty() {
                        mat.map_kd = Some(s);
                    }
                }
            }
            _ => {}
        }
    }
    if let Some(mat) = current.take() {
        lib.materials.push(mat);
    }
    lib
}

fn parse_rgb<'a>(parts: impl Iterator<Item = &'a str>) -> [f32; 3] {
    let vals: Vec<f32> = parts.take(3).filter_map(|s| s.parse::<f32>().ok()).collect();
    [
        *vals.first().unwrap_or(&0.0),
        *vals.get(1).unwrap_or(&0.0),
        *vals.get(2).unwrap_or(&0.0),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let err = obj_reader_extended(&PathBuf::from("a.ply")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn parses_simple_mtl() {
        let mtl = "newmtl Red\nKd 1 0 0\nKa 0.1 0.1 0.1\nKs 0.5 0.5 0.5\n";
        let lib = parse_mtl(mtl);
        assert_eq!(lib.materials.len(), 1);
        assert_eq!(lib.materials[0].name, "Red");
        assert_eq!(lib.materials[0].diffuse, [1.0, 0.0, 0.0]);
        assert_eq!(lib.materials[0].ambient, [0.1, 0.1, 0.1]);
        assert_eq!(lib.materials[0].specular, [0.5, 0.5, 0.5]);
    }

    #[test]
    fn parses_two_materials() {
        let mtl = "newmtl A\nKd 1 0 0\nnewmtl B\nKd 0 1 0\n";
        let lib = parse_mtl(mtl);
        assert_eq!(lib.materials.len(), 2);
        assert_eq!(lib.materials[0].name, "A");
        assert_eq!(lib.materials[1].name, "B");
    }

    #[test]
    fn parses_map_kd() {
        let mtl = "newmtl Tex\nKd 1 1 1\nmap_Kd albedo.png\n";
        let lib = parse_mtl(mtl);
        assert_eq!(lib.materials.len(), 1);
        assert_eq!(lib.materials[0].map_kd.as_deref(), Some("albedo.png"));
    }
}
