//! Opening cut helpers — produce a wall mesh with window / door
//! cut-outs applied.
//!
//! ## Strategy
//!
//! True mesh-domain boolean subtraction across triangles is fragile
//! (we'd need a CSG kernel). v1 takes a simpler approach: instead of
//! cutting the wall mesh, we emit a sequence of axis-aligned
//! sub-rectangles that tile the wall faces around each opening. The
//! result is visually equivalent for axis-aligned openings on
//! straight walls — every pixel inside an opening is left blank and
//! the rest is covered by triangles.
//!
//! Specifically, for each wall side face (+perp and -perp), we tile
//! the face into a grid of sub-rectangles by the union of all
//! opening x-ranges and z-ranges, then emit a triangle for every
//! sub-rect that's NOT inside any opening. Top, bottom, and end caps
//! are emitted untouched (openings don't pierce them by construction —
//! windows have a sill above the bottom, doors stop short of the
//! top).
//!
//! Limitations:
//! - Openings must lie entirely within the wall (no partial overlap
//!   off the end). We simply clip them silently.
//! - Two openings that overlap each other get merged.
//! - Openings that pierce the top of the wall (header height ≥ wall
//!   height) are clipped to the wall height — so a door of height
//!   = wall height still produces a closed top edge.

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::door::DoorParams;
use crate::error::ArchError;
use crate::wall::WallParams;
use crate::window::WindowParams;

/// One opening's footprint along the wall axis + z, in the wall's
/// local frame.
#[derive(Copy, Clone, Debug)]
struct Opening {
    /// Start position along the wall axis (clipped to [0, length]).
    x0: f64,
    /// End position along the wall axis.
    x1: f64,
    /// Bottom z relative to wall.start.z.
    z0: f64,
    /// Top z relative to wall.start.z.
    z1: f64,
}

impl Opening {
    fn from_window(w: &WindowParams, wall: &WallParams) -> Self {
        let half = w.width * 0.5;
        let x0 = (w.position_along_wall - half).clamp(0.0, wall.length());
        let x1 = (w.position_along_wall + half).clamp(0.0, wall.length());
        let z0 = w.position_height.clamp(0.0, wall.height);
        let z1 = (w.position_height + w.height).clamp(0.0, wall.height);
        Self { x0, x1, z0, z1 }
    }

    fn from_door(d: &DoorParams, wall: &WallParams) -> Self {
        let half = d.width * 0.5;
        let x0 = (d.position_along_wall - half).clamp(0.0, wall.length());
        let x1 = (d.position_along_wall + half).clamp(0.0, wall.length());
        let z0 = 0.0;
        let z1 = d.height.clamp(0.0, wall.height);
        Self { x0, x1, z0, z1 }
    }

    fn is_inside(&self, x: f64, z: f64) -> bool {
        x > self.x0 && x < self.x1 && z > self.z0 && z < self.z1
    }
}

/// Tessellate a wall with any number of window / door openings cut
/// from its side faces.
///
/// The returned mesh has:
/// - top, bottom, and start/end cap faces from the closed wall, and
/// - the two side faces (`+perp` and `-perp`) tessellated as a grid
///   of axis-aligned sub-rects with the opening cells skipped.
///
/// Pass empty `windows` + `doors` slices to get the same triangulation
/// as [`WallParams::tessellate_mesh`] (modulo grid splitting on the
/// side faces; the topology is equivalent under any z-buffer
/// renderer).
pub fn wall_with_openings(
    wall: &WallParams,
    windows: &[&WindowParams],
    doors: &[&DoorParams],
) -> Result<Mesh, ArchError> {
    wall.validate()?;
    let mut openings: Vec<Opening> = Vec::new();
    for w in windows {
        w.validate()?;
        openings.push(Opening::from_window(w, wall));
    }
    for d in doors {
        d.validate()?;
        openings.push(Opening::from_door(d, wall));
    }

    // Build the unique sorted x and z grid lines from the opening
    // edges plus the wall's outer extents.
    let mut xs: Vec<f64> = vec![0.0, wall.length()];
    let mut zs: Vec<f64> = vec![0.0, wall.height];
    for o in &openings {
        xs.push(o.x0);
        xs.push(o.x1);
        zs.push(o.z0);
        zs.push(o.z1);
    }
    sort_dedup(&mut xs);
    sort_dedup(&mut zs);

    let mut mesh = Mesh::new("wall_open");
    let mut block = ElementBlock::new(ElementType::Tri3);

    // Bottom, top, and the two end caps — emit from the wall's full
    // 8 corners.
    let c = wall.corners();
    for v in &c {
        mesh.nodes.push(*v);
    }
    // bottom 0-1-2-3.
    block.connectivity.extend_from_slice(&[0, 1, 2]);
    block.connectivity.extend_from_slice(&[0, 2, 3]);
    // top 4-7-6-5.
    block.connectivity.extend_from_slice(&[4, 7, 6]);
    block.connectivity.extend_from_slice(&[4, 6, 5]);
    // start cap 0-4-5-1.
    block.connectivity.extend_from_slice(&[0, 4, 5]);
    block.connectivity.extend_from_slice(&[0, 5, 1]);
    // end cap 3-2-6-7.
    block.connectivity.extend_from_slice(&[3, 2, 6]);
    block.connectivity.extend_from_slice(&[3, 6, 7]);

    // Side faces — tile (xs × zs) grid, skip cells inside an opening.
    let axis = wall.axis_xy();
    let perp = wall.perp_xy();
    let half_t = wall.thickness * 0.5;
    let z_base = wall.start.z;

    emit_side(
        &mut mesh,
        &mut block,
        &xs,
        &zs,
        &openings,
        wall.start,
        axis,
        perp * half_t,
        z_base,
        /* outward_normal_sign = */ 1.0,
    );
    emit_side(
        &mut mesh,
        &mut block,
        &xs,
        &zs,
        &openings,
        wall.start,
        axis,
        -perp * half_t,
        z_base,
        -1.0,
    );

    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(mesh)
}

fn sort_dedup(v: &mut Vec<f64>) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
}

#[allow(clippy::too_many_arguments)]
fn emit_side(
    mesh: &mut Mesh,
    block: &mut ElementBlock,
    xs: &[f64],
    zs: &[f64],
    openings: &[Opening],
    origin: Vector3<f64>,
    axis: Vector3<f64>,
    perp_offset: Vector3<f64>,
    z_base: f64,
    sign: f64,
) {
    // For each cell in (xs × zs) check whether its centre lies
    // inside an opening; if not, emit two triangles for it.
    for i in 0..(xs.len() - 1) {
        for k in 0..(zs.len() - 1) {
            let x0 = xs[i];
            let x1 = xs[i + 1];
            let z0 = zs[k];
            let z1 = zs[k + 1];
            let cx = 0.5 * (x0 + x1);
            let cz = 0.5 * (z0 + z1);
            let in_opening = openings.iter().any(|o| o.is_inside(cx, cz));
            if in_opening {
                continue;
            }
            // 4 corners of the sub-rect on this side face.
            let p_bl =
                origin + axis * x0 + perp_offset + Vector3::new(0.0, 0.0, z_base + z0 - origin.z);
            let p_br =
                origin + axis * x1 + perp_offset + Vector3::new(0.0, 0.0, z_base + z0 - origin.z);
            let p_tr =
                origin + axis * x1 + perp_offset + Vector3::new(0.0, 0.0, z_base + z1 - origin.z);
            let p_tl =
                origin + axis * x0 + perp_offset + Vector3::new(0.0, 0.0, z_base + z1 - origin.z);
            let n0 = mesh.nodes.len() as u32;
            mesh.nodes.push(p_bl);
            mesh.nodes.push(p_br);
            mesh.nodes.push(p_tr);
            mesh.nodes.push(p_tl);
            // Outward normal: for +sign (perp + side), CCW order is
            // bl, br, tr, tl; for -sign flip the winding.
            if sign > 0.0 {
                block.connectivity.extend_from_slice(&[n0, n0 + 1, n0 + 2]);
                block.connectivity.extend_from_slice(&[n0, n0 + 2, n0 + 3]);
            } else {
                block.connectivity.extend_from_slice(&[n0, n0 + 3, n0 + 2]);
                block.connectivity.extend_from_slice(&[n0, n0 + 2, n0 + 1]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::door::{DoorStyle, Side};
    use crate::window::WindowStyle;

    fn wall_5m() -> WallParams {
        WallParams {
            start: Vector3::zeros(),
            end: Vector3::new(5.0, 0.0, 0.0),
            height: 2.7,
            thickness: 0.2,
            material: "Brick".into(),
        }
    }

    #[test]
    fn no_openings_matches_closed_box_triangle_count_minimum() {
        // Without any openings, the side-face grid has 1 cell per
        // side → 2 tris × 2 sides = 4. Plus bottom (2) + top (2) +
        // 2 caps (4) = 12 — same as the closed wall.
        let m = wall_with_openings(&wall_5m(), &[], &[]).unwrap();
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn one_door_removes_at_least_one_cell() {
        let d = DoorParams {
            host: 1,
            position_along_wall: 2.5,
            width: 0.9,
            height: 2.1,
            style: DoorStyle::Single,
            hinge_side: Side::Left,
        };
        let m_closed = wall_with_openings(&wall_5m(), &[], &[]).unwrap();
        let m_open = wall_with_openings(&wall_5m(), &[], &[&d]).unwrap();
        // With a door, the side-face grid now has multiple cells; the
        // total triangle count includes the unchanged top/bottom/caps
        // but the side cells skip the door's footprint. With one
        // door, each side has (3 × 2) = 6 cells minus the (1 × 1)
        // door cell = 5 cells × 2 tris = 10 tris/side × 2 sides + 8
        // top/bot/caps = 28.
        assert!(m_open.total_elements() > m_closed.total_elements());
        // And it should still be valid.
        assert!(!m_open.nodes.is_empty());
    }

    #[test]
    fn window_and_door_together() {
        let w = WindowParams {
            host: 1,
            position_along_wall: 1.0,
            position_height: 1.0,
            width: 0.8,
            height: 1.0,
            frame_thickness: 0.05,
            style: WindowStyle::Casement,
        };
        let d = DoorParams {
            host: 1,
            position_along_wall: 3.5,
            width: 0.9,
            height: 2.1,
            style: DoorStyle::Single,
            hinge_side: Side::Left,
        };
        let m = wall_with_openings(&wall_5m(), &[&w], &[&d]).unwrap();
        assert!(m.total_elements() > 0);
    }
}
