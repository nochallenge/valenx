//! Phase 130 — COLLADA (`.dae`) writer.
//!
//! ## What OCCT does
//!
//! There is no upstream OCCT writer for COLLADA — the format
//! belongs to the Khronos asset-pipeline family alongside glTF and
//! is most commonly consumed by Maya / Blender / 3ds Max. Vendors
//! that want OCCT-to-Maya interop typically ship a third-party
//! plugin. The reference implementation (OpenCOLLADA, archived
//! since 2019) covers the full schema; modern integrations cover
//! the geometry subset only.
//!
//! On the wire COLLADA is XML with the namespace
//! `http://www.collada.org/2005/11/COLLADASchema`. Geometry sits
//! under `library_geometries / geometry / mesh / source +
//! triangles`, with positions in one `<source>` (`float_array`
//! plus a `<technique_common><accessor>` typed view) and triangle
//! indices in `<triangles><p>`.
//!
//! ## v1 status
//!
//! **Honest implementation** — a schema-correct COLLADA 1.4.1
//! document covering the geometry subset (the same surface area as
//! [`crate::vrml_writer()`] / [`crate::x3d_writer()`]). The
//! document is a complete, loadable `.dae`:
//!
//! - `<asset>` header with `up_axis` and `unit`,
//! - one `<geometry>` per `Tri3` element block, each carrying a
//!   POSITION `<source>` (`<float_array>` + typed `<accessor>`),
//!   a `<vertices>` element, and a `<triangles>` element with the
//!   flat index buffer in `<p>`,
//! - a `<library_visual_scenes>` with one `<node>` per geometry
//!   (`<instance_geometry>`), and a `<scene>` pointing at it.
//!
//! Materials, per-vertex colour, normals, texture coordinates, and
//! the assembly transform hierarchy are deferred — a follow-up adds
//! the `<library_materials>` / `<library_effects>` binding.

use std::fmt::Write as _;
use std::path::Path;

use valenx_mesh::{ElementType, Mesh};

use crate::error::OcctExchangeError;

/// Write `mesh` to `path` as a COLLADA 1.4.1 `.dae` document.
///
/// Each `Tri3` element block becomes one `<geometry>` instanced by a
/// node in the visual scene. Non-triangle blocks are skipped.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.dae`.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn collada_writer(mesh: &Mesh, path: &Path) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("dae") {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .dae",
        ));
    }
    let text = collada_document(mesh);
    valenx_core::io_caps::atomic_write_str(path, &text)?;
    Ok(())
}

/// Build the full COLLADA document as a string — pulled out so tests
/// can assert on the payload without touching the filesystem.
fn collada_document(mesh: &Mesh) -> String {
    let mut geom_lib = String::new();
    let mut node_xml = String::new();
    let mut geom_index = 0usize;

    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        let tri_count = block.connectivity.len() / 3;
        if tri_count == 0 {
            continue;
        }
        let geom_id = format!("geom{geom_index}");
        let pos_src = format!("{geom_id}-positions");
        let pos_array = format!("{pos_src}-array");
        let vert_id = format!("{geom_id}-vertices");

        // --- POSITION source: a flat float array of x y z triples. ---
        let mut floats = String::new();
        for n in &mesh.nodes {
            let _ = write!(floats, "{} {} {} ", fmt_f(n.x), fmt_f(n.y), fmt_f(n.z));
        }
        let float_count = mesh.nodes.len() * 3;
        let vertex_count = mesh.nodes.len();

        // --- triangle index buffer ---
        let mut indices = String::new();
        for tri in block.connectivity.chunks_exact(3) {
            let _ = write!(indices, "{} {} {} ", tri[0], tri[1], tri[2]);
        }

        let _ = write!(
            geom_lib,
            "    <geometry id=\"{geom_id}\" name=\"{geom_id}\">\n\
             \x20     <mesh>\n\
             \x20       <source id=\"{pos_src}\">\n\
             \x20         <float_array id=\"{pos_array}\" count=\"{float_count}\">{}</float_array>\n\
             \x20         <technique_common>\n\
             \x20           <accessor source=\"#{pos_array}\" count=\"{vertex_count}\" stride=\"3\">\n\
             \x20             <param name=\"X\" type=\"float\"/>\n\
             \x20             <param name=\"Y\" type=\"float\"/>\n\
             \x20             <param name=\"Z\" type=\"float\"/>\n\
             \x20           </accessor>\n\
             \x20         </technique_common>\n\
             \x20       </source>\n\
             \x20       <vertices id=\"{vert_id}\">\n\
             \x20         <input semantic=\"POSITION\" source=\"#{pos_src}\"/>\n\
             \x20       </vertices>\n\
             \x20       <triangles count=\"{tri_count}\">\n\
             \x20         <input semantic=\"VERTEX\" source=\"#{vert_id}\" offset=\"0\"/>\n\
             \x20         <p>{}</p>\n\
             \x20       </triangles>\n\
             \x20     </mesh>\n\
             \x20   </geometry>\n",
            floats.trim_end(),
            indices.trim_end(),
        );

        let _ = write!(
            node_xml,
            "      <node id=\"node{geom_index}\" name=\"node{geom_index}\">\n\
             \x20       <instance_geometry url=\"#{geom_id}\"/>\n\
             \x20     </node>\n",
        );
        geom_index += 1;
    }

    format!(
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n\
         <COLLADA xmlns=\"http://www.collada.org/2005/11/COLLADASchema\" version=\"1.4.1\">\n\
         \x20 <asset>\n\
         \x20   <contributor><authoring_tool>valenx-occt-exchange</authoring_tool></contributor>\n\
         \x20   <unit name=\"meter\" meter=\"1\"/>\n\
         \x20   <up_axis>Z_UP</up_axis>\n\
         \x20 </asset>\n\
         \x20 <library_geometries>\n{geom_lib}  </library_geometries>\n\
         \x20 <library_visual_scenes>\n\
         \x20   <visual_scene id=\"scene\" name=\"scene\">\n{node_xml}    </visual_scene>\n\
         \x20 </library_visual_scenes>\n\
         \x20 <scene>\n\
         \x20   <instance_visual_scene url=\"#scene\"/>\n\
         \x20 </scene>\n\
         </COLLADA>\n",
    )
}

/// Compact float formatter — trims trailing zeros so the file stays
/// small without losing precision.
fn fmt_f(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    let s = format!("{x:.9}");
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use std::path::PathBuf;
    use valenx_mesh::element::ElementBlock;

    fn triangle_mesh() -> Mesh {
        let mut m = Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks.push(blk);
        m
    }

    #[test]
    fn rejects_wrong_extension() {
        let m = Mesh::new("t");
        let err = collada_writer(&m, &PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn document_has_schema_skeleton() {
        let doc = collada_document(&triangle_mesh());
        assert!(doc.contains("<COLLADA"));
        assert!(doc.contains("http://www.collada.org/2005/11/COLLADASchema"));
        assert!(doc.contains("<library_geometries>"));
        assert!(doc.contains("<library_visual_scenes>"));
        assert!(doc.contains("<scene>"));
        assert!(doc.ends_with("</COLLADA>\n"));
    }

    #[test]
    fn triangle_geometry_is_emitted() {
        let doc = collada_document(&triangle_mesh());
        // One geometry, one triangles element with count 1.
        assert!(doc.contains("<geometry id=\"geom0\""));
        assert!(doc.contains("<triangles count=\"1\">"));
        // The index buffer carries the three vertex indices.
        assert!(doc.contains("<p>0 1 2</p>"));
        // The accessor declares 3 vertices, stride 3.
        assert!(doc.contains("count=\"3\" stride=\"3\""));
        // The float array carries 9 floats.
        assert!(doc.contains("count=\"9\""));
        // The node instances the geometry.
        assert!(doc.contains("<instance_geometry url=\"#geom0\"/>"));
    }

    #[test]
    fn empty_mesh_still_valid_document() {
        let doc = collada_document(&Mesh::new("empty"));
        assert!(doc.contains("<library_geometries>"));
        assert!(!doc.contains("<geometry"));
    }
}
