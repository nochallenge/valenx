//! Phase 129 — X3D (Extensible 3D, ISO/IEC 19775) writer.
//!
//! ## What OCCT does
//!
//! `RWX3D_Writer` emits the XML-syntax X3D form (`.x3d`), the
//! modern successor to VRML 2.0. The element vocabulary mirrors
//! VRML's: `Shape > IndexedFaceSet > Coordinate / coordIndex`,
//! plus optional `Appearance > Material / DiffuseColor` for
//! per-shape colour.
//!
//! ## v1 status
//!
//! **Honest implementation** for the minimal `Shape /
//! IndexedFaceSet` form analogous to [`crate::vrml_writer()`].
//! Header carries the X3D 3.3 profile declaration and `Scene`
//! root. Appearance / per-vertex colour deferred to Phase 129.5.

use std::path::Path;

use valenx_mesh::{ElementType, Mesh};

use crate::error::OcctExchangeError;

/// Write `mesh` to `path` as X3D XML.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.x3d`.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn x3d_writer(mesh: &Mesh, path: &Path) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("x3d") {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .x3d",
        ));
    }
    let mut text = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.3//EN\" \"http://www.web3d.org/specifications/x3d-3.3.dtd\">\n\
<X3D profile=\"Interchange\" version=\"3.3\">\n\
  <head><meta name=\"generator\" content=\"valenx-occt-exchange\"/></head>\n\
  <Scene>\n",
    );
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        let mut coord = String::new();
        for n in &mesh.nodes {
            coord.push_str(&format!("{} {} {} ", n.x, n.y, n.z));
        }
        let mut coord_index = String::new();
        for tri in block.connectivity.chunks_exact(3) {
            coord_index.push_str(&format!("{} {} {} -1 ", tri[0], tri[1], tri[2]));
        }
        text.push_str("    <Shape>\n");
        text.push_str(&format!(
            "      <IndexedFaceSet coordIndex=\"{}\"><Coordinate point=\"{}\"/></IndexedFaceSet>\n",
            coord_index.trim_end(),
            coord.trim_end(),
        ));
        text.push_str("    </Shape>\n");
    }
    text.push_str("  </Scene>\n</X3D>\n");
    valenx_core::io_caps::atomic_write_str(path, &text)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let m = Mesh::new("t");
        let err = x3d_writer(&m, &PathBuf::from("a.wrl")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
