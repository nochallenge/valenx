//! Waterline 3D — Z-level finish on a 3D solid surface.
//!
//! Slices the source mesh at constant Z levels and emits a profile
//! cut at each level. Differs from [`crate::op::profile`] in that
//! every level is independent (no carry-through of the boundary) and
//! the source is expected to be a fully-3D solid (not a 2.5D
//! extrusion).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::cut::intersect_plane;
use valenx_mesh::Mesh;

use crate::error::CamError;
use crate::op::profile::stitch_segments;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Waterline 3D parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Waterline3DParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Z-level spacing (mm).
    pub step_down: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for Waterline3DParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 800.0,
            plunge_feed: 200.0,
            spindle_rpm: 12_000.0,
            step_down: 0.5,
            depth: 5.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a waterline-3D toolpath.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &Waterline3DParams,
    _tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    let n_levels = crate::op::compute_n_passes(params.depth, params.step_down, "waterline_3d")?;
    let mut emitted = false;
    for k in 1..=n_levels {
        let depth = (params.step_down * k as f64).min(params.depth);
        let z = stock.top_z() - depth;
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
        if polygon.is_empty() {
            continue;
        }
        let start = polygon[0];
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
        for v in polygon.iter().skip(1) {
            tp.push(Move::new(
                MoveKind::Cut,
                Vector3::new(v.x, v.y, z),
                params.feed_mm_per_min,
            ));
        }
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(start.x, start.y, z),
            params.feed_mm_per_min,
        ));
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(start.x, start.y, safe_z),
            0.0,
        ));
        emitted = true;
    }
    if !emitted {
        return Err(CamError::BadOperation {
            name: "waterline_3d".into(),
            reason: "no slice yielded a usable contour".into(),
        });
    }
    Ok(tp)
}

fn validate(params: &Waterline3DParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "waterline_3d".into(),
        reason,
    };
    if !(params.step_down > 0.0) {
        return Err(mk(format!(
            "step_down must be > 0 (got {})",
            params.step_down
        )));
    }
    if !(params.depth > 0.0) {
        return Err(mk(format!("depth must be > 0 (got {})", params.depth)));
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(mk("feed must be > 0".into()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;
    use valenx_mesh::element::{ElementBlock, ElementType};
    use valenx_mesh::Mesh;

    fn cube(size: f64) -> Mesh {
        let s = size * 0.5;
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
    fn waterline_emits_per_level() {
        let mesh = cube(10.0);
        let stock = Stock::new(
            Vector3::new(-6.0, -6.0, -5.0),
            Vector3::new(12.0, 12.0, 10.0),
            "alu",
        )
        .unwrap();
        let tool = Tool::new(1, "BM3", ToolKind::BallMill, 3.0, 25.0, 2, "").unwrap();
        let params = Waterline3DParams {
            step_down: 1.0,
            depth: 4.0,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).unwrap();
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert!(
            plunges >= 4,
            "expected >=4 plunges (one per level), got {plunges}"
        );
    }
}
