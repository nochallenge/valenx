//! Phase 124 — extended STL writer (binary + ASCII modes).
//!
//! ## What OCCT does
//!
//! `RWStl_Reader::WriteAscii` / `RWStl_Reader::WriteBinary` emit
//! the same triangle soup in either of the two STL flavours. ASCII
//! is human-readable and widely supported for tiny test inputs;
//! binary is ~4-5x more compact and faster to parse, and what every
//! production tool prefers. OCCT picks ASCII unless the caller
//! flips a flag.
//!
//! ## v1 status
//!
//! **Honest implementation** for both modes. Binary delegates to
//! [`valenx_mesh::stl_write::write_stl_binary`] (Phase 7
//! production writer). ASCII is a hand-rolled emitter at this
//! crate — the `valenx-mesh` crate doesn't currently expose ASCII
//! STL, and Phase 124 is the natural home for it given the
//! data-exchange theme.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use valenx_mesh::{ElementType, Mesh};

use crate::error::OcctExchangeError;

/// Output format selector for [`stl_writer_extended`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum StlMode {
    /// Binary STL (80-byte header + LE triangle records).
    Binary,
    /// ASCII STL (`solid` / `facet` text format).
    Ascii,
}

/// Write `mesh` to `path` as STL in the selected mode.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.stl`.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn stl_writer_extended(
    mesh: &Mesh,
    mode: StlMode,
    path: &Path,
) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("stl") {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .stl",
        ));
    }
    match mode {
        StlMode::Binary => valenx_mesh::stl_write::write_stl_binary(mesh, path)?,
        StlMode::Ascii => write_ascii(mesh, path)?,
    }
    Ok(())
}

fn write_ascii(mesh: &Mesh, path: &Path) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut sink = BufWriter::new(file);
    writeln!(sink, "solid valenx_mesh")?;
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let v0 = mesh.nodes[tri[0] as usize];
            let v1 = mesh.nodes[tri[1] as usize];
            let v2 = mesh.nodes[tri[2] as usize];
            let edge1 = v1 - v0;
            let edge2 = v2 - v0;
            let mut normal = edge1.cross(&edge2);
            let len = normal.norm();
            if len > 1e-20 {
                normal /= len;
            } else {
                normal[0] = 0.0;
                normal[1] = 0.0;
                normal[2] = 1.0;
            }
            writeln!(
                sink,
                "  facet normal {} {} {}",
                normal[0], normal[1], normal[2],
            )?;
            sink.write_all(b"    outer loop\n")?;
            for v in [v0, v1, v2] {
                writeln!(sink, "      vertex {} {} {}", v[0], v[1], v[2])?;
            }
            sink.write_all(b"    endloop\n")?;
            sink.write_all(b"  endfacet\n")?;
        }
    }
    writeln!(sink, "endsolid valenx_mesh")?;
    sink.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let m = Mesh::new("t");
        let err = stl_writer_extended(&m, StlMode::Binary, &PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
