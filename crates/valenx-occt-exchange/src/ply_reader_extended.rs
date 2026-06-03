//! Phase 123 — extended PLY reader (vertex colours + normals).
//!
//! ## What OCCT does
//!
//! `RWPly_Reader` parses a Stanford PLY header, identifies the
//! `red/green/blue` and `nx/ny/nz` properties (when present), and
//! exposes them alongside the geometry. The OCCT-side mesh stores
//! them as per-vertex attributes that downstream renderers consume.
//!
//! ## v1 status
//!
//! **Honest implementation** — the geometry path delegates to
//! [`valenx_mesh::format::ply::read_path`] (which tolerates the extra
//! vertex columns by reading + discarding them, and as of Phase 26.5
//! reads binary PLY too); the colour / normal pass re-scans the same
//! file with a minimal ASCII parser to recover the discarded values.
//! For binary PLY the geometry imports but the colour / normal
//! annotations come back empty — recovering them needs a typed binary
//! rescan, a small follow-up.

use std::path::Path;

use valenx_mesh::Mesh;

use crate::error::OcctExchangeError;
use crate::ply_writer_extended::PlyVertexAnnotations;

/// Geometry + colour / normal annotations recovered from a PLY file.
#[derive(Clone, Debug)]
pub struct PlyImport {
    /// Imported triangle-surface mesh.
    pub mesh: Mesh,
    /// Per-vertex colour / normal — empty fields when the file
    /// didn't carry them.
    pub annotations: PlyVertexAnnotations,
}

/// Read a `.ply` file plus any per-vertex colour / normal columns.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.ply`.
/// - [`OcctExchangeError::Parse`] for malformed PLY.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn ply_reader_extended(path: &Path) -> Result<PlyImport, OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("ply") {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .ply",
        ));
    }
    let mesh = valenx_mesh::format::ply::read_path(path)
        .map_err(|e| OcctExchangeError::parse("ply file", format!("{e}")))?;
    // The colour / normal rescan is ASCII-only. Binary PLB geometry
    // now reads via the underlying reader (Phase 26.5), but a binary
    // body is not UTF-8 — fall back to empty annotations rather than
    // surfacing an Io error for an otherwise-valid binary file.
    //
    // Round-21 M1: bound the annotation rescan at
    // MAX_CAD_INTERCHANGE_FILE_BYTES (256 MiB — sister to the STEP /
    // IGES cap). Pre-fix a hostile multi-GB ASCII `.ply` would have
    // slurped here before the annotation parser ran.
    let annotations = match valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES as usize,
    ) {
        Ok(text) => parse_annotations(&text),
        Err(_) => PlyVertexAnnotations::default(),
    };
    Ok(PlyImport { mesh, annotations })
}

/// Pure ASCII-PLY annotation extractor — pulled out so tests can
/// feed in synthetic text.
fn parse_annotations(text: &str) -> PlyVertexAnnotations {
    // Walk the header to discover where vertex columns sit and
    // whether colour / normal columns are present.
    let mut lines = text.lines();
    let mut vertex_count: usize = 0;
    let mut properties: Vec<String> = Vec::new();
    let mut in_vertex_block = false;
    for line in lines.by_ref() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("element ") {
            in_vertex_block = rest.starts_with("vertex ");
            if in_vertex_block {
                vertex_count = rest
                    .strip_prefix("vertex ")
                    .and_then(|c| c.trim().parse::<usize>().ok())
                    .unwrap_or(0);
            }
            continue;
        }
        if in_vertex_block {
            if let Some(rest) = trimmed.strip_prefix("property ") {
                // "property <type> <name>" — keep <name>.
                let mut parts = rest.split_whitespace();
                let _ty = parts.next();
                if let Some(name) = parts.next() {
                    properties.push(name.to_string());
                }
            }
        }
        if trimmed == "end_header" {
            break;
        }
    }
    let normal_idx = properties.iter().position(|p| p == "nx").and_then(|i| {
        if properties.get(i + 1).is_some_and(|p| p == "ny")
            && properties.get(i + 2).is_some_and(|p| p == "nz")
        {
            Some(i)
        } else {
            None
        }
    });
    let color_idx = properties.iter().position(|p| p == "red").and_then(|i| {
        if properties.get(i + 1).is_some_and(|p| p == "green")
            && properties.get(i + 2).is_some_and(|p| p == "blue")
        {
            Some(i)
        } else {
            None
        }
    });
    if normal_idx.is_none() && color_idx.is_none() {
        return PlyVertexAnnotations::default();
    }
    let mut colors = Vec::new();
    let mut normals = Vec::new();
    for (i, line) in lines.enumerate() {
        if i >= vertex_count {
            break;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        if let Some(ni) = normal_idx {
            let nx = cols.get(ni).and_then(|s| s.parse::<f32>().ok()).unwrap_or(0.0);
            let ny = cols
                .get(ni + 1)
                .and_then(|s| s.parse::<f32>().ok())
                .unwrap_or(0.0);
            let nz = cols
                .get(ni + 2)
                .and_then(|s| s.parse::<f32>().ok())
                .unwrap_or(0.0);
            normals.push([nx, ny, nz]);
        }
        if let Some(ci) = color_idx {
            let r = cols.get(ci).and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
            let g = cols.get(ci + 1).and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
            let b = cols.get(ci + 2).and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
            colors.push([r, g, b]);
        }
    }
    PlyVertexAnnotations { colors, normals }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let err = ply_reader_extended(&PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn parses_colors() {
        let ply = "\
ply
format ascii 1.0
element vertex 2
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 0
property list uchar int vertex_indices
end_header
0 0 0 255 0 0
1 0 0 0 255 0
";
        let ann = parse_annotations(ply);
        assert_eq!(ann.colors, vec![[255, 0, 0], [0, 255, 0]]);
        assert!(ann.normals.is_empty());
    }

    #[test]
    fn parses_normals() {
        let ply = "\
ply
format ascii 1.0
element vertex 1
property float x
property float y
property float z
property float nx
property float ny
property float nz
element face 0
property list uchar int vertex_indices
end_header
0 0 0 0 0 1
";
        let ann = parse_annotations(ply);
        assert_eq!(ann.normals, vec![[0.0, 0.0, 1.0]]);
        assert!(ann.colors.is_empty());
    }

    #[test]
    fn returns_empty_for_plain_ply() {
        let ply = "\
ply
format ascii 1.0
element vertex 1
property float x
property float y
property float z
element face 0
property list uchar int vertex_indices
end_header
0 0 0
";
        let ann = parse_annotations(ply);
        assert!(ann.colors.is_empty());
        assert!(ann.normals.is_empty());
    }
}
