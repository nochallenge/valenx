//! Spiral pocket — pocket cleared by an Archimedean spiral instead
//! of a zig-zag raster.
//!
//! Reuses the existing pocket polygon extraction; difference vs.
//! [`crate::op::pocket`] with `PocketStrategy::Spiral` is that the
//! spiral is generated *analytically* from `step_over` around a
//! single centre, rather than chained inward offsets.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::cut::intersect_plane;
use valenx_mesh::Mesh;

use crate::error::CamError;
use crate::op::profile::stitch_segments;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Spiral-pocket parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpiralPocketParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Step-over between spiral turns (mm).
    pub step_over: f64,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Segments per spiral revolution (≥3).
    pub n_steps_per_rev: u32,
}

impl Default for SpiralPocketParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 500.0,
            plunge_feed: 150.0,
            spindle_rpm: 12_000.0,
            step_over: 1.5,
            step_down: 1.0,
            depth: 3.0,
            safe_z_clearance: 5.0,
            n_steps_per_rev: 36,
        }
    }
}

fn polygon_centroid(polygon: &[Vector3<f64>]) -> Vector3<f64> {
    let n = polygon.len() as f64;
    let mut s = Vector3::zeros();
    for v in polygon {
        s += *v;
    }
    s / n
}

fn polygon_max_radius_from(centre: Vector3<f64>, polygon: &[Vector3<f64>]) -> f64 {
    let mut max = 0.0_f64;
    for v in polygon {
        let dx = v.x - centre.x;
        let dy = v.y - centre.y;
        let r = (dx * dx + dy * dy).sqrt();
        if r > max {
            max = r;
        }
    }
    max
}

/// Generate a spiral-pocket toolpath.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &SpiralPocketParams,
    tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "spiral_pocket")?;
    let mut emitted = false;
    for k in 1..=n_passes {
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
        let centre = polygon_centroid(&polygon);
        let r_max = polygon_max_radius_from(centre, &polygon) - tool.radius_mm();
        if r_max <= 0.0 {
            continue;
        }
        // Build an Archimedean spiral from centre outward.
        let n = params.n_steps_per_rev.max(3) as usize;
        let dtheta = std::f64::consts::TAU / n as f64;
        let dr = params.step_over / n as f64;
        let mut theta = 0.0_f64;
        let mut r = 0.0_f64;
        // Plunge at the centre.
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(centre.x, centre.y, safe_z),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(centre.x, centre.y, z),
            params.plunge_feed,
        ));
        while r <= r_max {
            theta += dtheta;
            r += dr;
            let x = centre.x + r * theta.cos();
            let y = centre.y + r * theta.sin();
            tp.push(Move::new(
                MoveKind::Cut,
                Vector3::new(x, y, z),
                params.feed_mm_per_min,
            ));
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(centre.x, centre.y, safe_z),
            0.0,
        ));
        emitted = true;
    }
    if !emitted {
        return Err(CamError::BadOperation {
            name: "spiral_pocket".into(),
            reason: "no slice yielded a usable pocket polygon".into(),
        });
    }
    Ok(tp)
}

fn validate(params: &SpiralPocketParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "spiral_pocket".into(),
        reason,
    };
    if !(params.step_over > 0.0) {
        return Err(mk(format!(
            "step_over must be > 0 (got {})",
            params.step_over
        )));
    }
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
    if params.n_steps_per_rev < 3 {
        return Err(mk("n_steps_per_rev must be >= 3".into()));
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
    fn spiral_pocket_emits_path() {
        let mesh = cube(20.0);
        let stock = Stock::new(
            Vector3::new(-12.0, -12.0, -10.0),
            Vector3::new(24.0, 24.0, 20.0),
            "alu",
        )
        .unwrap();
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "").unwrap();
        let params = SpiralPocketParams {
            step_down: 2.0,
            depth: 2.0,
            step_over: 1.0,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts > 4, "expected spiral cuts, got {cuts}");
    }
}
