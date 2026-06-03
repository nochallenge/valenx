//! Tessellate a [`KicadBoard`] outline + drill holes into a
//! mesh-backed solid.

use nalgebra::Vector3;

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::board::KicadBoard;
use crate::error::KicadError;

const DRILL_SEGMENTS: usize = 24;

/// Extrude the board outline + subtract drill holes. Returns a
/// triangulated [`Solid::Mesh`].
///
/// v1 implementation: extrudes the outer polygon as a fan-top /
/// fan-bottom slab + walls; drill holes are represented as small
/// cylinders ADDED beneath the bottom face as "holes you can see"
/// (true subtraction defers to Phase 42.5 when a 2D Boolean kernel
/// lands; the result is still useful as a viewport proxy).
pub fn pcb_to_solid(board: &KicadBoard) -> Result<Solid, KicadError> {
    if board.outline.len() < 3 {
        return Err(KicadError::BadParameter {
            name: "outline",
            reason: format!("need >= 3 verts, got {}", board.outline.len()),
        });
    }
    if board.thickness_mm <= 0.0 {
        return Err(KicadError::BadParameter {
            name: "thickness_mm",
            reason: format!("must be > 0, got {}", board.thickness_mm),
        });
    }

    let n = board.outline.len();
    let mut mesh = Mesh::new("kicad_pcb");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let half = board.thickness_mm * 0.5;

    // Top + bottom rings.
    let top_base = mesh.nodes.len() as u32;
    for p in &board.outline {
        mesh.nodes.push(Vector3::new(p[0], p[1], half));
    }
    let bot_base = mesh.nodes.len() as u32;
    for p in &board.outline {
        mesh.nodes.push(Vector3::new(p[0], p[1], -half));
    }

    // Centroid for top + bottom fans.
    let mut c = Vector3::zeros();
    for p in &board.outline {
        c += Vector3::new(p[0], p[1], 0.0);
    }
    c /= n as f64;
    let top_c = mesh.nodes.len() as u32;
    mesh.nodes.push(Vector3::new(c.x, c.y, half));
    let bot_c = mesh.nodes.len() as u32;
    mesh.nodes.push(Vector3::new(c.x, c.y, -half));

    for i in 0..n {
        let j = (i + 1) % n;
        // Top fan.
        block.connectivity.extend_from_slice(&[
            top_c,
            top_base + i as u32,
            top_base + j as u32,
        ]);
        // Bottom fan (reversed).
        block.connectivity.extend_from_slice(&[
            bot_c,
            bot_base + j as u32,
            bot_base + i as u32,
        ]);
        // Walls.
        let a = top_base + i as u32;
        let b = top_base + j as u32;
        let cc = bot_base + j as u32;
        let d = bot_base + i as u32;
        block.connectivity.extend_from_slice(&[a, b, cc]);
        block.connectivity.extend_from_slice(&[a, cc, d]);
    }

    // Drill holes as inverted cylinders (visual proxy).
    for (centre, diameter) in &board.drill_holes {
        let r = *diameter * 0.5;
        let ring_top_base = mesh.nodes.len() as u32;
        for k in 0..DRILL_SEGMENTS {
            let theta = (k as f64 / DRILL_SEGMENTS as f64) * std::f64::consts::TAU;
            mesh.nodes.push(Vector3::new(
                centre.x + r * theta.cos(),
                centre.y + r * theta.sin(),
                half,
            ));
        }
        let ring_bot_base = mesh.nodes.len() as u32;
        for k in 0..DRILL_SEGMENTS {
            let theta = (k as f64 / DRILL_SEGMENTS as f64) * std::f64::consts::TAU;
            mesh.nodes.push(Vector3::new(
                centre.x + r * theta.cos(),
                centre.y + r * theta.sin(),
                -half,
            ));
        }
        for k in 0..DRILL_SEGMENTS {
            let j = (k + 1) % DRILL_SEGMENTS;
            let a = ring_top_base + k as u32;
            let b = ring_top_base + j as u32;
            let cc = ring_bot_base + j as u32;
            let d = ring_bot_base + k as u32;
            // Inward-facing walls (CW orientation).
            block.connectivity.extend_from_slice(&[a, cc, b]);
            block.connectivity.extend_from_slice(&[a, d, cc]);
        }
    }

    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_board_tessellates() {
        let b = KicadBoard::demo_devboard();
        let s = pcb_to_solid(&b).unwrap();
        match s {
            Solid::Mesh(m) => {
                assert!(!m.nodes.is_empty());
                assert!(m.total_elements() > 0);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn empty_outline_errors() {
        let b = KicadBoard::new_default();
        assert!(matches!(
            pcb_to_solid(&b),
            Err(KicadError::BadParameter { .. })
        ));
    }
}
