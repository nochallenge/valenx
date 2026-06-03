//! Phase 125 — extended STL reader (binary + ASCII auto-detect).
//!
//! ## What OCCT does
//!
//! `RWStl_Reader::ReadFile` sniffs the first ~80 bytes: a leading
//! `b"solid "` token *plus* the absence of a NUL byte in the first
//! few KB means ASCII; otherwise it falls back to binary. Some
//! exotic binary STL writers put `solid ` in their header which
//! makes the sniff non-trivial — OCCT uses a hybrid check
//! (decimal-fraction count consistency) for robustness.
//!
//! ## v1 status
//!
//! **Honest implementation** for the geometry path. Both binary
//! and ASCII variants are parsed here directly (the `valenx-mesh`
//! crate doesn't yet ship an STL reader). The detector follows the
//! OCCT-style hybrid heuristic: if the file is exactly
//! `80 + 4 + n*50` bytes long for some `n` that matches the
//! declared triangle count, treat it as binary; otherwise try
//! ASCII.

use std::path::Path;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::OcctExchangeError;

/// Read an STL file from `path`, auto-detecting binary vs ASCII.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.stl`.
/// - [`OcctExchangeError::Parse`] for malformed input (wrong byte
///   count, malformed facet/loop text).
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn stl_reader_extended(path: &Path) -> Result<Mesh, OcctExchangeError> {
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
    // Round-21 M1: bound the read at MAX_CAD_INTERCHANGE_FILE_BYTES
    // (256 MiB — sister to the STEP / IGES / glTF caps). Pre-fix
    // this was a bare `fs::read(path)`; a hostile multi-GB STL
    // would have slurped before the binary-vs-ASCII sniff ran.
    let bytes = valenx_core::io_caps::read_capped_to_bytes(
        path,
        valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES,
    )?;
    let id = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("stl")
        .to_string();
    if is_binary(&bytes) {
        parse_binary(id, &bytes)
    } else {
        let text = std::str::from_utf8(&bytes)
            .map_err(|e| OcctExchangeError::parse("stl ascii", format!("not utf-8: {e}")))?;
        parse_ascii(id, text)
    }
}

/// Hybrid binary-vs-ASCII detector — see module header.
fn is_binary(bytes: &[u8]) -> bool {
    if bytes.len() < 84 {
        return false;
    }
    // The 4 bytes following the 80-byte header carry the declared
    // triangle count for binary STL. Predict the total file size.
    let tri_count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
    let expected = 80 + 4 + tri_count * 50;
    if bytes.len() == expected {
        return true;
    }
    // Some writers (notably Blender) emit a binary header that starts
    // with `solid ` even though the file is binary. The size match
    // above takes care of that case unambiguously.
    if bytes.starts_with(b"solid ") && !bytes[..bytes.len().min(2048)].contains(&0) {
        return false;
    }
    // Trust the size heuristic over the ASCII tag.
    bytes.len() == expected
}

fn parse_binary(id: String, bytes: &[u8]) -> Result<Mesh, OcctExchangeError> {
    if bytes.len() < 84 {
        return Err(OcctExchangeError::parse(
            "stl binary",
            "file shorter than 84 bytes (no header + count)",
        ));
    }
    let tri_count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
    let expected = 80 + 4 + tri_count * 50;
    if bytes.len() != expected {
        return Err(OcctExchangeError::parse(
            "stl binary",
            format!(
                "expected {expected} bytes for {tri_count} triangles; got {}",
                bytes.len(),
            ),
        ));
    }
    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(tri_count * 3);
    let mut conn: Vec<u32> = Vec::with_capacity(tri_count * 3);
    let mut cursor = 84;
    for _ in 0..tri_count {
        // Skip normal (12 bytes), read three vertices (36 bytes),
        // skip attribute (2 bytes).
        cursor += 12;
        for _ in 0..3 {
            let x = f32::from_le_bytes([
                bytes[cursor],
                bytes[cursor + 1],
                bytes[cursor + 2],
                bytes[cursor + 3],
            ]) as f64;
            let y = f32::from_le_bytes([
                bytes[cursor + 4],
                bytes[cursor + 5],
                bytes[cursor + 6],
                bytes[cursor + 7],
            ]) as f64;
            let z = f32::from_le_bytes([
                bytes[cursor + 8],
                bytes[cursor + 9],
                bytes[cursor + 10],
                bytes[cursor + 11],
            ]) as f64;
            conn.push(nodes.len() as u32);
            nodes.push(Vector3::new(x, y, z));
            cursor += 12;
        }
        cursor += 2;
    }
    let mut mesh = Mesh::new(id);
    mesh.nodes = nodes;
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(mesh)
}

fn parse_ascii(id: String, text: &str) -> Result<Mesh, OcctExchangeError> {
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut conn: Vec<u32> = Vec::new();
    let mut current_tri: Vec<u32> = Vec::with_capacity(3);
    for (i, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("vertex ") {
            let coords: Vec<f64> = rest
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if coords.len() != 3 {
                return Err(OcctExchangeError::parse(
                    format!("stl ascii line {}", i + 1),
                    format!("vertex needs 3 coords, got {}", coords.len()),
                ));
            }
            let idx = nodes.len() as u32;
            nodes.push(Vector3::new(coords[0], coords[1], coords[2]));
            current_tri.push(idx);
        } else if trimmed == "endloop" {
            if current_tri.len() != 3 {
                return Err(OcctExchangeError::parse(
                    format!("stl ascii line {}", i + 1),
                    format!("loop has {} vertices, expected 3", current_tri.len()),
                ));
            }
            conn.extend_from_slice(&current_tri);
            current_tri.clear();
        }
    }
    let mut mesh = Mesh::new(id);
    mesh.nodes = nodes;
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let err = stl_reader_extended(&PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn parses_minimal_ascii_triangle() {
        let stl = "\
solid test
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
      vertex 0 1 0
    endloop
  endfacet
endsolid test
";
        let mesh = parse_ascii("t".into(), stl).unwrap();
        assert_eq!(mesh.nodes.len(), 3);
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(mesh.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn rejects_truncated_ascii_loop() {
        let stl = "\
solid test
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
    endloop
  endfacet
endsolid test
";
        let err = parse_ascii("t".into(), stl).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.parse");
    }

    #[test]
    fn parses_minimal_binary_zero_triangles() {
        let mut bytes = vec![0u8; 84];
        // 0 triangles → 80 + 4 = 84 bytes total → valid binary.
        let mesh = parse_binary("t".into(), &bytes).unwrap();
        assert_eq!(mesh.nodes.len(), 0);
        // Setting wrong size triggers an error.
        bytes.push(0);
        let err = parse_binary("t".into(), &bytes).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.parse");
    }

    #[test]
    fn is_binary_recognises_size_match() {
        let mut bytes = vec![0u8; 84];
        // declared triangle count = 0 → expected size = 84
        assert!(is_binary(&bytes));
        // Now bump it by one byte; mismatch → not binary.
        bytes.push(0);
        assert!(!is_binary(&bytes));
    }
}
