//! Quantized-OBJ <-> valenx Mesh codec.

use crate::model::ModelProfile;
use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// Decode a (possibly partial) quantized-OBJ token string into a valenx `Mesh`.
///
/// Two-pass + tolerant: pass 1 collects every well-formed `v` line (dequantized
/// via `profile`); pass 2 fan-triangulates every `f` line, skipping any face
/// that references an out-of-range vertex. Malformed lines are skipped, never
/// fatal — the model output is untrusted text.
pub fn decode(text: &str, profile: &ModelProfile) -> Mesh {
    let mut mesh = Mesh::new("meshgen");

    // Pass 1: vertices.
    for line in text.lines() {
        let mut it = line.split_whitespace();
        if it.next() != Some("v") {
            continue;
        }
        let coords: Vec<f64> = it
            .take(3)
            .filter_map(|t| t.parse::<i64>().ok())
            .map(|q| profile.dequant(q))
            .collect();
        if coords.len() == 3 {
            mesh.nodes
                .push(Vector3::new(coords[0], coords[1], coords[2]));
        }
    }

    // Pass 2: faces (1-based OBJ indices -> 0-based, fan-triangulated).
    let n = mesh.nodes.len() as u32;
    let mut block = ElementBlock::new(ElementType::Tri3);
    for line in text.lines() {
        let mut it = line.split_whitespace();
        if it.next() != Some("f") {
            continue;
        }
        let verts: Vec<u32> = it
            // OBJ face tokens may be `v`, `v/vt`, `v//vn`; take the vertex part.
            .filter_map(|tok| tok.split('/').next()?.parse::<i64>().ok())
            // 1-based, positive only; convert to 0-based.
            .filter_map(|i| u32::try_from(i - 1).ok())
            .collect();
        if verts.len() < 3 {
            continue;
        }
        for k in 1..verts.len() - 1 {
            let (a, b, c) = (verts[0], verts[k], verts[k + 1]);
            if a < n && b < n && c < n {
                block.connectivity.extend_from_slice(&[a, b, c]);
            }
        }
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

/// Encode a `Mesh` to quantized-OBJ text using `profile`'s grid. Only Tri3
/// blocks are emitted (the mesh-LLM format is triangle soup); other element
/// types are skipped. Faces are written 1-based, per OBJ.
pub fn encode(mesh: &Mesh, profile: &ModelProfile) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    for v in &mesh.nodes {
        let _ = writeln!(
            out,
            "v {} {} {}",
            profile.quant(v.x),
            profile.quant(v.y),
            profile.quant(v.z)
        );
    }
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let _ = writeln!(out, "f {} {} {}", tri[0] + 1, tri[1] + 1, tri[2] + 1);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // A unit tetrahedron expressed in LLAMA_MESH bins (128-bin grid, -1..1):
    // bin 0 -> -1, bin 127 -> +1, bin 64 -> ~0.008. Faces are 1-based, fan ok.
    const TET: &str = "\
v 0 0 0
v 127 0 0
v 0 127 0
v 0 0 127
f 1 2 3
f 1 2 4
f 1 3 4
f 2 3 4
";

    #[test]
    fn decode_reads_vertices_and_triangulates_faces() {
        let m = decode(TET, &ModelProfile::LLAMA_MESH);
        assert_eq!(m.nodes.len(), 4, "4 vertices");
        let tris: usize = m
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len() / 3)
            .sum();
        assert_eq!(tris, 4, "4 triangles");
        // bin 0 -> coord_min (-1.0)
        assert!((m.nodes[0].x - (-1.0)).abs() < 1e-9);
        // bin 127 -> coord_max (+1.0)
        assert!((m.nodes[1].x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn decode_fan_triangulates_a_quad() {
        let quad = "v 0 0 0\nv 127 0 0\nv 127 127 0\nv 0 127 0\nf 1 2 3 4\n";
        let m = decode(quad, &ModelProfile::LLAMA_MESH);
        let tris: usize = m
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len() / 3)
            .sum();
        assert_eq!(tris, 2, "a quad fans into 2 triangles");
    }

    #[test]
    fn decode_skips_garbage_and_out_of_range_faces_without_panicking() {
        let junk = "\
hello world
v 0 0 0
v not a number
v 127 0 0
v 0 127 0
f 1 2 999
f 1 2 3
garbage line
f
";
        let m = decode(junk, &ModelProfile::LLAMA_MESH);
        // Only the 3 well-formed `v` lines parse.
        assert_eq!(m.nodes.len(), 3);
        // `f 1 2 999` is dropped (999 out of range); `f 1 2 3` is kept; `f` is dropped.
        let tris: usize = m
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len() / 3)
            .sum();
        assert_eq!(tris, 1);
    }

    #[test]
    fn decode_empty_input_is_an_empty_mesh() {
        let m = decode("", &ModelProfile::LLAMA_MESH);
        assert_eq!(m.nodes.len(), 0);
        assert_eq!(
            m.element_blocks
                .iter()
                .map(|b| b.connectivity.len())
                .sum::<usize>(),
            0
        );
    }

    #[test]
    fn encode_then_decode_preserves_topology() {
        let original = decode(TET, &ModelProfile::LLAMA_MESH);
        let text = encode(&original, &ModelProfile::LLAMA_MESH);
        let round = decode(&text, &ModelProfile::LLAMA_MESH);
        assert_eq!(round.nodes.len(), original.nodes.len());
        let t0: usize = original
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len())
            .sum();
        let t1: usize = round
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len())
            .sum();
        assert_eq!(t1, t0);
        // Coords survive a quant->dequant round-trip to within one bin width.
        let bin = (ModelProfile::LLAMA_MESH.coord_max - ModelProfile::LLAMA_MESH.coord_min)
            / (ModelProfile::LLAMA_MESH.quant_bins as f64 - 1.0);
        for (a, b) in original.nodes.iter().zip(round.nodes.iter()) {
            assert!((a - b).norm() <= bin * 2.0);
        }
    }

    #[test]
    fn encode_emits_v_and_f_lines() {
        let m = decode(TET, &ModelProfile::LLAMA_MESH);
        let text = encode(&m, &ModelProfile::LLAMA_MESH);
        assert_eq!(text.lines().filter(|l| l.starts_with("v ")).count(), 4);
        assert_eq!(text.lines().filter(|l| l.starts_with("f ")).count(), 4);
    }
}
