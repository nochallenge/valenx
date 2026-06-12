//! Phase 126 — glTF 2.0 writer (JSON-only embedded form).
//!
//! ## What OCCT does
//!
//! `RWGltf_CafWriter` walks an `XCAFDoc_DocumentTool` document and
//! emits a glTF 2.0 file. The Khronos glTF 2.0 spec defines two
//! container shapes:
//!
//! - `.gltf` — the manifest is JSON with `bufferUri` pointing at an
//!   external `.bin` blob (or a data URI for embedded payload).
//! - `.glb` — a single-file binary container with a length-prefixed
//!   JSON chunk followed by a binary chunk.
//!
//! Either way, geometry sits in a `buffer` as packed little-endian
//! `float32 positions`, `float32 normals` (optional), and
//! `uint32 indices`. `accessor` records map ranges of the buffer to
//! `mesh.primitives`. `mesh` references `accessor`s, `node`s
//! reference `mesh`es, the top-level `scene` references `node`s.
//!
//! ## v1 status
//!
//! **Honest implementation** for the JSON-only embedded form
//! (`.gltf`). The geometry is base64-encoded into a single
//! `data:application/octet-stream;base64,...` URI in the `buffers`
//! array. This is the simplest glTF flavour and is what every glTF
//! 2.0 reader supports; the trade-off is ~33% size overhead from
//! base64. Phase 126.5 will add `.glb` support (binary chunking).

use std::path::Path;

use valenx_mesh::{ElementType, Mesh};

use crate::error::OcctExchangeError;

/// Write `mesh` to `path` as glTF 2.0 with the buffer payload
/// embedded as a base64 data URI.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.gltf`.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn gltf2_writer(mesh: &Mesh, path: &Path) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    if ext.as_deref() != Some("gltf") {
        return Err(OcctExchangeError::bad_input(
            "path",
            "extension must be .gltf (use .glb for binary, deferred to Phase 126.5)",
        ));
    }
    // Pack vertex positions into a contiguous float32 buffer.
    let mut buf: Vec<u8> = Vec::with_capacity(mesh.nodes.len() * 12);
    let mut min_xyz = [f32::INFINITY; 3];
    let mut max_xyz = [f32::NEG_INFINITY; 3];
    for n in &mesh.nodes {
        let x = n.x as f32;
        let y = n.y as f32;
        let z = n.z as f32;
        buf.extend_from_slice(&x.to_le_bytes());
        buf.extend_from_slice(&y.to_le_bytes());
        buf.extend_from_slice(&z.to_le_bytes());
        for (i, &v) in [x, y, z].iter().enumerate() {
            if v < min_xyz[i] {
                min_xyz[i] = v;
            }
            if v > max_xyz[i] {
                max_xyz[i] = v;
            }
        }
    }
    let pos_byte_length = buf.len();
    // Append uint32 indices for all Tri3 blocks.
    let mut index_count: u32 = 0;
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for &i in &block.connectivity {
            buf.extend_from_slice(&i.to_le_bytes());
            index_count += 1;
        }
    }
    let idx_byte_length = buf.len() - pos_byte_length;
    let total_byte_length = buf.len();
    let b64 = base64_encode(&buf);
    let json = format!(
        "{{\
\"asset\":{{\"version\":\"2.0\",\"generator\":\"valenx-occt-exchange\"}},\
\"scene\":0,\
\"scenes\":[{{\"nodes\":[0]}}],\
\"nodes\":[{{\"mesh\":0}}],\
\"meshes\":[{{\"primitives\":[{{\"attributes\":{{\"POSITION\":0}},\"indices\":1,\"mode\":4}}]}}],\
\"buffers\":[{{\"uri\":\"data:application/octet-stream;base64,{b64}\",\"byteLength\":{total_byte_length}}}],\
\"bufferViews\":[\
{{\"buffer\":0,\"byteOffset\":0,\"byteLength\":{pos_byte_length},\"target\":34962}},\
{{\"buffer\":0,\"byteOffset\":{pos_byte_length},\"byteLength\":{idx_byte_length},\"target\":34963}}\
],\
\"accessors\":[\
{{\"bufferView\":0,\"componentType\":5126,\"count\":{vertex_count},\"type\":\"VEC3\",\"min\":[{minx},{miny},{minz}],\"max\":[{maxx},{maxy},{maxz}]}},\
{{\"bufferView\":1,\"componentType\":5125,\"count\":{index_count},\"type\":\"SCALAR\"}}\
]\
}}",
        vertex_count = mesh.nodes.len(),
        minx = min_xyz[0],
        miny = min_xyz[1],
        minz = min_xyz[2],
        maxx = max_xyz[0],
        maxy = max_xyz[1],
        maxz = max_xyz[2],
    );
    valenx_core::io_caps::atomic_write_str(path, &json)?;
    Ok(())
}

/// Tiny base64 encoder — std doesn't ship one and we don't want a
/// new dependency just for this writer. Standard alphabet, no
/// padding shortcut (always pads to multiple of 4).
fn base64_encode(input: &[u8]) -> String {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(input.len().div_ceil(3) * 4);
    for chunk in input.chunks(3) {
        let b0 = chunk[0];
        let b1 = chunk.get(1).copied().unwrap_or(0);
        let b2 = chunk.get(2).copied().unwrap_or(0);
        let triple = ((b0 as u32) << 16) | ((b1 as u32) << 8) | b2 as u32;
        out.push(ALPHABET[((triple >> 18) & 0x3F) as usize] as char);
        out.push(ALPHABET[((triple >> 12) & 0x3F) as usize] as char);
        if chunk.len() >= 2 {
            out.push(ALPHABET[((triple >> 6) & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
        if chunk.len() >= 3 {
            out.push(ALPHABET[(triple & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let m = Mesh::new("t");
        let err = gltf2_writer(&m, &PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn base64_known_vectors() {
        // RFC 4648 test vectors.
        assert_eq!(base64_encode(b""), "");
        assert_eq!(base64_encode(b"f"), "Zg==");
        assert_eq!(base64_encode(b"fo"), "Zm8=");
        assert_eq!(base64_encode(b"foo"), "Zm9v");
        assert_eq!(base64_encode(b"foob"), "Zm9vYg==");
        assert_eq!(base64_encode(b"fooba"), "Zm9vYmE=");
        assert_eq!(base64_encode(b"foobar"), "Zm9vYmFy");
    }
}
