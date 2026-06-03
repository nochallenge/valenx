//! Pocket operation — hollow out the cross-section interior at each
//! Z step-down level.
//!
//! ## Algorithm
//!
//! 1. For each Z level (top of stock, decreasing by `step_down`
//!    until reaching `top_z() - depth`):
//!    1. Cut the source mesh by a horizontal plane at that Z and
//!       stitch the resulting line segments into a closed polygon
//!       (same path as [`crate::op::profile`]).
//!    2. Offset the polygon inward by `tool.radius_mm()` so the
//!       cutter centre stays inside the part boundary.
//!    3. Fill the offset polygon by [`PocketStrategy`]:
//!       - `ZigZag` — alternating-direction raster lines.
//!       - `Parallel` — one-way raster (lift between passes).
//!       - `Spiral` — concentric inward rings.
//!    4. Emit rapid → plunge → cut → rapid sequences.
//!
//! ## v1 simplifications
//!
//! - **Single closed polygon per Z level** — see [`crate::op::profile`].
//! - **Step-over enforcement** — warns (but does not reject) if
//!   `step_over > tool.diameter * 0.5`. v1 trusts the caller.
//! - **No lead-in / lead-out arcs.** Plunges go straight in at full
//!   plunge feed.

use nalgebra::Vector3;
use valenx_mesh::cut::intersect_plane;
use valenx_mesh::Mesh;

use crate::error::CamError;
use crate::offset;
use crate::op::profile::stitch_segments;
use crate::operation::{PocketParams, PocketStrategy};
use crate::raster;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Generate a pocket-clear toolpath. See module docs for the algorithm
/// and v1 simplifications.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &PocketParams,
    tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params, tool)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));

    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "pocket")?;
    let mut had_polygon = false;
    for k in 1..=n_passes {
        let depth_below_top = (params.step_down * k as f64).min(params.depth);
        let z = stock.top_z() - depth_below_top;
        let segments = intersect_plane(
            source,
            Vector3::new(0.0, 0.0, z),
            Vector3::new(0.0, 0.0, 1.0),
        );
        if segments.is_empty() {
            continue;
        }
        let polygon = match stitch_segments(&segments) {
            Some(p) => p,
            None => continue,
        };
        // Inward offset by tool radius — pocket centre must stay
        // inside the cavity.
        let offsets = offset::polygon(&polygon, -tool.radius_mm());
        if offsets.is_empty() {
            // Pocket too small for tool — skip this level.
            continue;
        }
        let mut ring = offsets.into_iter().next().unwrap();
        for v in &mut ring {
            v.z = z;
        }
        match params.strategy {
            PocketStrategy::ZigZag => {
                let mut lines = raster::zigzag(&ring, params.step_over, params.raster_angle_deg);
                if !params.climb {
                    for line in &mut lines {
                        line.reverse();
                    }
                }
                emit_zigzag(&mut tp, &lines, z, safe_z, params);
            }
            PocketStrategy::Parallel => {
                let lines = raster::parallel(&ring, params.step_over, params.raster_angle_deg);
                emit_parallel(&mut tp, &lines, z, safe_z, params);
            }
            PocketStrategy::Spiral => {
                let spiral = raster::spiral(&ring, params.step_over);
                emit_spiral(&mut tp, &spiral, z, safe_z, params);
            }
        }
        had_polygon = true;
    }
    if !had_polygon {
        return Err(CamError::BadOperation {
            name: "pocket".into(),
            reason: "no pocket polygon could be extracted from the source mesh".into(),
        });
    }
    Ok(tp)
}

fn emit_zigzag(
    tp: &mut Toolpath,
    lines: &[Vec<Vector3<f64>>],
    z: f64,
    safe_z: f64,
    params: &PocketParams,
) {
    if lines.is_empty() {
        return;
    }
    let first = lines[0][0];
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(first.x, first.y, safe_z),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Plunge,
        Vector3::new(first.x, first.y, z),
        params.plunge_feed,
    ));
    // First line — cut through it.
    for v in lines[0].iter().skip(1) {
        tp.push(Move::new(MoveKind::Cut, *v, params.feed_mm_per_min));
    }
    // Each subsequent line: cut directly to its start (the zigzag
    // alternation already reverses every other line so consecutive
    // line endpoints are adjacent).
    for line in lines.iter().skip(1) {
        for v in line {
            tp.push(Move::new(MoveKind::Cut, *v, params.feed_mm_per_min));
        }
    }
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
}

fn emit_parallel(
    tp: &mut Toolpath,
    lines: &[Vec<Vector3<f64>>],
    z: f64,
    safe_z: f64,
    params: &PocketParams,
) {
    for line in lines {
        if line.len() < 2 {
            continue;
        }
        let start = line[0];
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(start.x, start.y, safe_z),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(start.x, start.y, z),
            params.plunge_feed,
        ));
        for v in line.iter().skip(1) {
            tp.push(Move::new(MoveKind::Cut, *v, params.feed_mm_per_min));
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(v_last_or(line).x, v_last_or(line).y, safe_z),
            0.0,
        ));
    }
}

fn emit_spiral(
    tp: &mut Toolpath,
    spiral: &[Vector3<f64>],
    z: f64,
    safe_z: f64,
    params: &PocketParams,
) {
    if spiral.is_empty() {
        return;
    }
    let first = spiral[0];
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(first.x, first.y, safe_z),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Plunge,
        Vector3::new(first.x, first.y, z),
        params.plunge_feed,
    ));
    for v in spiral.iter().skip(1) {
        let mut p = *v;
        p.z = z;
        tp.push(Move::new(MoveKind::Cut, p, params.feed_mm_per_min));
    }
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
}

fn v_last_or(line: &[Vector3<f64>]) -> Vector3<f64> {
    *line.last().unwrap_or(&Vector3::zeros())
}

fn validate(params: &PocketParams, tool: &Tool) -> Result<(), CamError> {
    if !(params.step_down > 0.0) {
        return Err(CamError::BadOperation {
            name: "pocket".into(),
            reason: format!("step_down must be > 0 (got {})", params.step_down),
        });
    }
    if !(params.step_over > 0.0) {
        return Err(CamError::BadOperation {
            name: "pocket".into(),
            reason: format!("step_over must be > 0 (got {})", params.step_over),
        });
    }
    if !(params.depth > 0.0) {
        return Err(CamError::BadOperation {
            name: "pocket".into(),
            reason: format!("depth must be > 0 (got {})", params.depth),
        });
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(CamError::BadOperation {
            name: "pocket".into(),
            reason: format!("feed must be > 0 (got {})", params.feed_mm_per_min),
        });
    }
    // Tool engagement check — warn (do not reject) on aggressive
    // step-over. v1 trusts the caller.
    if params.step_over > tool.diameter_mm * 0.5 {
        tracing::warn!(
            "pocket step_over {:.3} > tool.diameter/2 = {:.3}; tool engagement may exceed safe limit",
            params.step_over,
            tool.diameter_mm * 0.5
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operation::{PocketParams, PocketStrategy};
    use crate::tool::{Tool, ToolKind};

    fn make_cube() -> valenx_mesh::Mesh {
        // Reuse the same cube mesh factory via path import.
        use valenx_mesh::element::{ElementBlock, ElementType};
        use valenx_mesh::Mesh;
        let s = 5.0_f64;
        let nodes = vec![
            Vector3::new(-s, -s, -s),
            Vector3::new(s, -s, -s),
            Vector3::new(s, s, -s),
            Vector3::new(-s, s, -s),
            Vector3::new(-s, -s, s),
            Vector3::new(s, -s, s),
            Vector3::new(s, s, s),
            Vector3::new(-s, s, s),
        ];
        let conn: Vec<u32> = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut mesh = Mesh::new("cube");
        mesh.nodes = nodes;
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = conn;
        mesh.element_blocks.push(block);
        mesh
    }

    #[test]
    fn pocket_zigzag_produces_raster_pattern() {
        let mesh = make_cube(); // 10x10x10 cube — we pretend it's the pocket boundary
        let stock = Stock::new(
            Vector3::new(-6.0, -6.0, -5.0),
            Vector3::new(12.0, 12.0, 10.0),
            "wood",
        )
        .unwrap();
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "carbide").unwrap();
        let params = PocketParams {
            step_over: 1.0,
            step_down: 2.0,
            depth: 4.0,
            strategy: PocketStrategy::ZigZag,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).expect("should generate");
        assert!(
            tp.len() > 10,
            "expected substantial toolpath, got {}",
            tp.len()
        );
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert!(plunges >= 2, "expected ≥2 plunges, got {plunges}");
    }

    #[test]
    fn pocket_spiral_strategy() {
        let mesh = make_cube();
        let stock = Stock::new(
            Vector3::new(-6.0, -6.0, -5.0),
            Vector3::new(12.0, 12.0, 10.0),
            "wood",
        )
        .unwrap();
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "carbide").unwrap();
        let params = PocketParams {
            step_over: 1.0,
            step_down: 2.0,
            depth: 2.0,
            strategy: PocketStrategy::Spiral,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).expect("should generate");
        assert!(!tp.is_empty());
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts > 5, "expected substantial spiral cuts, got {cuts}");
    }

    #[test]
    fn pocket_validate_rejects_zero_step_over() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let bad = PocketParams {
            step_over: 0.0,
            ..Default::default()
        };
        assert!(validate(&bad, &tool).is_err());
    }
}
