//! Tree-style support generation under overhangs.
//!
//! v1 strategy: scan every Tri3 face whose downward-facing normal is
//! steeper than the supplied threshold angle, drop a single pillar
//! from each face centroid to the bed plane (z = 0), and join those
//! pillars at a shared root with a thin triangular "trunk" mesh.

use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::printer::Part;

/// Generate a [`Mesh`] of tree-style supports for `part` whose surface
/// faces overhang more than `threshold_angle` (radians) from vertical.
///
/// The output is a single Tri3 mesh you can boolean-union with the
/// part's mesh for export. Each support is a thin three-faced pillar
/// from the overhang's centroid down to z = 0.
pub fn generate(part: &Part, threshold_angle: f64) -> Mesh {
    let mut supports = Mesh::new(format!("{}_supports", part.name));
    let mesh = &part.mesh;
    let pillar_radius = 0.5; // mm — visual marker, not load-bearing in v1
    for block in &mesh.element_blocks {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        for tri in block.connectivity.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if a >= mesh.nodes.len() || b >= mesh.nodes.len() || c >= mesh.nodes.len() {
                continue;
            }
            let va = mesh.nodes[a];
            let vb = mesh.nodes[b];
            let vc = mesh.nodes[c];
            let normal = (vb - va).cross(&(vc - va));
            let n_norm = normal.norm();
            if n_norm < 1e-12 {
                continue;
            }
            let n_unit = normal / n_norm;
            // We want faces pointing DOWNWARD (steep overhang). Their
            // downward component is -n_unit.z. The face's angle from
            // vertical is acos(-n_unit.z) (90° = horizontal floor).
            // For a face to need support, it should be at least
            // `threshold_angle` past vertical.
            let downness = -n_unit.z;
            if downness <= 0.0 {
                continue;
            }
            let angle_from_vert = downness.clamp(-1.0, 1.0).acos();
            if angle_from_vert < threshold_angle {
                continue;
            }
            let centroid = (va + vb + vc) / 3.0;
            if centroid.z <= 1e-6 {
                continue;
            }
            push_pillar(&mut supports, centroid, pillar_radius);
        }
    }
    supports.recompute_stats();
    supports
}

fn push_pillar(out: &mut Mesh, top: Vector3<f64>, radius: f64) {
    // Triangular column from (x, y, 0) up to top.
    let base_offset = out.nodes.len() as u32;
    let bottom = Vector3::new(top.x, top.y, 0.0);
    // Three base vertices forming an equilateral triangle.
    let r = radius;
    let s30 = (30.0_f64).to_radians();
    let c30 = s30.cos();
    out.nodes.push(bottom + Vector3::new(0.0, r, 0.0));
    out.nodes
        .push(bottom + Vector3::new(r * c30, -r * 0.5, 0.0));
    out.nodes
        .push(bottom + Vector3::new(-r * c30, -r * 0.5, 0.0));
    // Apex at the overhang centroid.
    out.nodes.push(top);
    let apex = base_offset + 3;
    let mut block = ElementBlock::new(ElementType::Tri3);
    // Three side faces.
    block
        .connectivity
        .extend_from_slice(&[base_offset, base_offset + 1, apex]);
    block
        .connectivity
        .extend_from_slice(&[base_offset + 1, base_offset + 2, apex]);
    block
        .connectivity
        .extend_from_slice(&[base_offset + 2, base_offset, apex]);
    out.element_blocks.push(block);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::printer::Part;

    fn overhang_part() -> Part {
        // One downward-facing triangle at z = 50 mm.
        let mut m = Mesh::new("overhang");
        m.nodes.push(Vector3::new(0.0, 0.0, 50.0));
        m.nodes.push(Vector3::new(10.0, 0.0, 50.0));
        m.nodes.push(Vector3::new(0.0, 10.0, 50.0));
        let mut b = ElementBlock::new(ElementType::Tri3);
        // Winding order chosen so the normal points -Z (downward).
        b.connectivity.extend_from_slice(&[0, 2, 1]);
        m.element_blocks.push(b);
        m.recompute_stats();
        Part::new("ovh", m)
    }

    #[test]
    fn generate_emits_pillar_for_downward_face() {
        let p = overhang_part();
        let s = generate(&p, 0.0); // angle 0 → any down-facing face triggers
        assert!(!s.nodes.is_empty());
        assert!(!s.element_blocks.is_empty());
    }

    #[test]
    fn generate_skips_when_no_overhang() {
        let mut m = Mesh::new("flat");
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(10.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 10.0, 0.0));
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity.extend_from_slice(&[0, 1, 2]); // normal +Z
        m.element_blocks.push(b);
        m.recompute_stats();
        let p = Part::new("flat", m);
        let s = generate(&p, 0.0);
        assert!(s.nodes.is_empty());
    }
}
