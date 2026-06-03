//! Binary STL writer for canonical [`Mesh`]es with `Tri3` blocks.
//!
//! ## Format
//!
//! Per the de-facto STL binary format:
//!
//! - 80-byte header (ASCII or NUL-padded — we write a NUL-padded
//!   marker so consumers that ASCII-sniff the first bytes correctly
//!   treat the file as binary).
//! - 4-byte little-endian `u32` triangle count.
//! - `N × 50` byte records of `{ normal: f32×3, v0: f32×3, v1: f32×3,
//!   v2: f32×3, attribute: u16 }`.
//!
//! ## Scope
//!
//! Only `ElementType::Tri3` element blocks contribute triangles —
//! STL is a triangle-soup format and doesn't represent volume
//! elements. Non-`Tri3` blocks are silently skipped, so a mixed
//! mesh exports just its surface.
//!
//! Normals are computed from each triangle's winding order
//! (right-hand rule) at write time — the canonical [`Mesh`] doesn't
//! cache per-face normals, and recomputing on write keeps the
//! output consistent with whatever transformations the mesh has
//! been through.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::element::ElementType;
use crate::mesh::Mesh;

/// Write a binary STL file of every `Tri3` block in `mesh` to
/// `path`. Errors propagate from the underlying file IO; callers
/// typically wrap a result-bearing `?` over the write.
///
/// The 80-byte header carries a `b"valenx-mesh binary STL"` marker
/// followed by NUL padding — distinctive enough to recognise our
/// output in a hex dump, short enough to fit without truncation.
pub fn write_stl_binary(mesh: &Mesh, path: impl AsRef<Path>) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Header — 80 bytes, NUL-padded. STL binary readers must ignore
    // header contents; only the 4-byte triangle count that follows
    // matters for parsing.
    let mut header = [0u8; 80];
    let tag = b"valenx-mesh binary STL";
    header[..tag.len()].copy_from_slice(tag);
    writer.write_all(&header)?;

    // Two-pass approach: first count Tri3 triangles so we can emit
    // the right declared count, then walk again writing each record.
    // Cheap — connectivity arrays are flat Vec<u32>.
    //
    // R34 S2 (defense-in-depth): both passes apply the SAME
    // `tri_in_range` predicate so the declared count and the records
    // written stay consistent. A triangle citing a vertex past
    // `nodes.len()` is dropped (graceful degrade) rather than panicking
    // `mesh.nodes[..]` — skipping a degenerate triangle is the correct
    // behaviour for a writer. The per-loader parse guards (OBJ/gmsh/
    // netgen/PLY) are the first line; this seal backs them so a future
    // un-hardened loader still produces a valid (if smaller) STL.
    let tri_in_range = |tri: &[u32]| {
        (tri[0] as usize) < mesh.nodes.len()
            && (tri[1] as usize) < mesh.nodes.len()
            && (tri[2] as usize) < mesh.nodes.len()
    };
    let total_tris: u32 = mesh
        .element_blocks
        .iter()
        .filter(|b| b.element_type == ElementType::Tri3)
        .map(|b| {
            b.connectivity
                .chunks_exact(3)
                .filter(|tri| tri_in_range(tri))
                .count() as u32
        })
        .sum();
    writer.write_all(&total_tris.to_le_bytes())?;

    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            // Same predicate as the count pass — keep them in lockstep.
            let (Some(&v0), Some(&v1), Some(&v2)) = (
                mesh.nodes.get(tri[0] as usize),
                mesh.nodes.get(tri[1] as usize),
                mesh.nodes.get(tri[2] as usize),
            ) else {
                continue;
            };
            let edge1 = v1 - v0;
            let edge2 = v2 - v0;
            let mut normal = edge1.cross(&edge2);
            let len = normal.norm();
            if len > 1e-20 {
                normal /= len;
            } else {
                // Degenerate triangle — STL allows any normal; emit
                // a default rather than a NaN-laden value.
                normal[0] = 0.0;
                normal[1] = 0.0;
                normal[2] = 1.0;
            }
            for &c in normal.as_slice() {
                writer.write_all(&(c as f32).to_le_bytes())?;
            }
            for v in [v0, v1, v2] {
                for &c in v.as_slice() {
                    writer.write_all(&(c as f32).to_le_bytes())?;
                }
            }
            writer.write_all(&0u16.to_le_bytes())?; // attribute
        }
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::{ElementBlock, ElementType};
    use nalgebra::Vector3;
    use std::io::Read;

    fn tri_mesh() -> Mesh {
        let mut m = Mesh::new("tri-out");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    #[test]
    fn writes_header_count_and_record() {
        let dir = tempdir();
        let path = dir.join("tri.stl");
        write_stl_binary(&tri_mesh(), &path).expect("write");
        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();
        // 80 + 4 + 50 = 134 bytes.
        assert_eq!(buf.len(), 134, "got {}", buf.len());
        // Header carries the marker.
        assert_eq!(&buf[..22], b"valenx-mesh binary STL");
        // Triangle count is 1.
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count, 1);
    }

    #[test]
    fn skips_non_tri3_blocks() {
        // Mesh with one Tri3 and one Hex8 — only the Tri3 should
        // contribute.
        let mut m = tri_mesh();
        let mut hex = ElementBlock::new(ElementType::Hex8);
        hex.connectivity = vec![0, 1, 2, 0, 0, 1, 2, 0];
        m.element_blocks.push(hex);
        let dir = tempdir();
        let path = dir.join("mixed.stl");
        write_stl_binary(&m, &path).expect("write");
        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count, 1, "Hex8 must not contribute");
    }

    /// R34 S2 (RED→GREEN): defense-in-depth sink seal. A mesh whose
    /// Tri3 connectivity cites a vertex past `nodes.len()` must NOT
    /// panic the writer. Pre-fix the write loop did
    /// `mesh.nodes[tri[0] as usize]` and panicked "index out of
    /// bounds". Post-fix the bad triangle is dropped, the declared
    /// count matches the records actually written, and the file stays
    /// a valid STL. We assert: one valid triangle survives, the
    /// out-of-range one is skipped, count == 1, file length == 134.
    #[test]
    fn out_of_range_triangle_is_skipped_count_consistent() {
        let mut m = Mesh::new("hostile");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // First tri valid; second cites vertex 9 (out of range).
        blk.connectivity = vec![0, 1, 2, 0, 1, 9];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        let dir = tempdir();
        let path = dir.join("hostile.stl");
        // Must not panic.
        write_stl_binary(&m, &path).expect("write");
        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();
        // Exactly the one valid triangle is declared and written:
        // 80 (header) + 4 (count) + 50 (one record) = 134 bytes.
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count, 1, "only the valid triangle should be declared");
        assert_eq!(
            buf.len(),
            134,
            "declared count must match the records written (no orphan/missing record)"
        );
    }

    #[test]
    fn empty_mesh_writes_just_header_and_zero_count() {
        let m = Mesh::new("empty");
        let dir = tempdir();
        let path = dir.join("empty.stl");
        write_stl_binary(&m, &path).expect("write");
        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();
        assert_eq!(buf.len(), 84);
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count, 0);
    }

    fn tempdir() -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        p.push(format!("valenx-mesh-stl-test-{nanos}"));
        std::fs::create_dir_all(&p).unwrap();
        p
    }
}
