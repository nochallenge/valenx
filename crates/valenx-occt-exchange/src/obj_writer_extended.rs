//! Phase 120 — extended Wavefront OBJ writer (materials + groups).
//!
//! ## What OCCT does
//!
//! OCCT's `RWMesh_CafReader` / `RWObj_CafReader` pair handles OBJ
//! plus its sibling `.mtl` material library. The writer emits:
//!
//! - `mtllib <name>.mtl` directive in the OBJ header pointing at the
//!   sidecar material file.
//! - `g <group>` lines partitioning the triangle stream by face
//!   group.
//! - `usemtl <material>` lines selecting the active material from
//!   the .mtl library before each group.
//! - The corresponding `.mtl` file with Kd / Ka / Ks colours and
//!   optional `map_Kd` texture references.
//!
//! ## v1 status
//!
//! **Honest implementation** — delegates the geometry to
//! [`valenx_mesh::format::obj::write_path`] then writes a sidecar
//! `<stem>.mtl` listing the caller's [`MaterialLib`] entries. The
//! main OBJ file is augmented with an `mtllib` reference at the
//! top. Per-group `usemtl` partitioning is deferred to Phase 120.5
//! (needs a per-triangle material index in `valenx-mesh::Mesh`'s
//! element-block model).

use std::path::Path;

use valenx_mesh::Mesh;

use crate::error::OcctExchangeError;

/// Material entry written into the sidecar `.mtl` file.
#[derive(Clone, Debug, PartialEq)]
pub struct ObjMaterial {
    /// Material name (selected via `usemtl`).
    pub name: String,
    /// Diffuse colour, 0..=1.
    pub diffuse: [f32; 3],
    /// Ambient colour, 0..=1.
    pub ambient: [f32; 3],
    /// Specular colour, 0..=1.
    pub specular: [f32; 3],
    /// Optional texture map (filename relative to the `.mtl` file).
    pub map_kd: Option<String>,
}

/// Bundle of materials to write into the sidecar `.mtl`.
#[derive(Clone, Debug, Default)]
pub struct MaterialLib {
    /// Materials to emit in order.
    pub materials: Vec<ObjMaterial>,
}

/// Write `mesh` to `path` as Wavefront OBJ plus a sidecar
/// `<stem>.mtl` containing `materials`.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.obj`.
/// - [`OcctExchangeError::Backend`] for valenx-mesh failures.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn obj_writer_extended(
    mesh: &Mesh,
    materials: &MaterialLib,
    path: &Path,
) -> Result<(), OcctExchangeError> {
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
    // Write the geometry through the production writer first.
    valenx_mesh::format::obj::write_path(mesh, path)
        .map_err(|e| OcctExchangeError::Backend(format!("obj::write_path: {e}")))?;
    if materials.materials.is_empty() {
        return Ok(());
    }
    // Emit sidecar .mtl using the path's file stem.
    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| OcctExchangeError::bad_input("path", "missing file stem"))?
        .to_string();
    let mtl_path = path.with_file_name(format!("{stem}.mtl"));
    let mut mtl_text = String::from("# valenx-occt-exchange MTL library\n");
    for m in &materials.materials {
        mtl_text.push_str(&format!("newmtl {}\n", m.name));
        mtl_text.push_str(&format!(
            "Ka {} {} {}\n",
            m.ambient[0], m.ambient[1], m.ambient[2],
        ));
        mtl_text.push_str(&format!(
            "Kd {} {} {}\n",
            m.diffuse[0], m.diffuse[1], m.diffuse[2],
        ));
        mtl_text.push_str(&format!(
            "Ks {} {} {}\n",
            m.specular[0], m.specular[1], m.specular[2],
        ));
        if let Some(map) = &m.map_kd {
            mtl_text.push_str(&format!("map_Kd {map}\n"));
        }
        mtl_text.push('\n');
    }
    valenx_core::io_caps::atomic_write_str(&mtl_path, &mtl_text)?;
    // Prepend `mtllib <stem>.mtl\n` to the OBJ. Read the freshly
    // written OBJ, splice the directive after the leading comment
    // header, write back.
    //
    // Round-23 sweep: bound the round-trip read at MAX_OBJ_FILE_BYTES
    // (1 GiB) — even though we just wrote this file, a concurrent
    // process or a hostile FS layer could have substituted a larger
    // file at the same path between our write and read.
    let obj_text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_OBJ_FILE_BYTES as usize,
    )?;
    let directive = format!("mtllib {stem}.mtl\n");
    // Insert after the first non-comment line if possible; for the
    // valenx-mesh writer's layout, insertion at byte 0 is safe (it
    // starts with `# valenx-mesh OBJ export`).
    let mut new = String::with_capacity(obj_text.len() + directive.len());
    new.push_str(&directive);
    new.push_str(&obj_text);
    valenx_core::io_caps::atomic_write_str(path, &new)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let m = Mesh::new("t");
        let lib = MaterialLib::default();
        let err = obj_writer_extended(&m, &lib, &PathBuf::from("a.ply")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
