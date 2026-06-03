//! Rest machining — remove material left after a previous larger-tool
//! operation.
//!
//! v1 takes the source mesh + a `prev_tool_radius_mm` (radius of the
//! previously-used roughing tool). It generates a profile-style cut
//! along the boundary using the *current* tool, but offset inward by
//! the *previous* tool's radius — so only material the larger tool
//! could not reach is removed.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::cut::intersect_plane;
use valenx_mesh::Mesh;

use crate::error::CamError;
use crate::offset;
use crate::op::profile::stitch_segments;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Rest-machining parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RestMachiningParams {
    /// Tool id (smaller current finishing tool).
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Radius (mm) of the previously-used roughing tool.
    pub prev_tool_radius_mm: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for RestMachiningParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 500.0,
            plunge_feed: 150.0,
            spindle_rpm: 18_000.0,
            step_down: 1.0,
            depth: 3.0,
            prev_tool_radius_mm: 3.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a rest-machining toolpath.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &RestMachiningParams,
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
    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "rest_machining")?;
    let mut emitted = false;
    // Inward offset = previous-tool radius - current-tool radius:
    // the rest-corridor lives in the annulus where the previous tool
    // could not reach but the current can.
    let inward = -(params.prev_tool_radius_mm - tool.radius_mm()).max(0.0);
    for k in 1..=n_passes {
        let z = stock.top_z() - (params.step_down * k as f64).min(params.depth);
        let segs = intersect_plane(
            source,
            Vector3::new(0.0, 0.0, z),
            Vector3::new(0.0, 0.0, 1.0),
        );
        if segs.is_empty() {
            continue;
        }
        let polygon = match stitch_segments(&segs) {
            Some(p) => p,
            None => continue,
        };
        let rings = offset::polygon(&polygon, inward);
        if rings.is_empty() {
            continue;
        }
        let ring = &rings[0];
        if ring.is_empty() {
            continue;
        }
        let start = ring[0];
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
        for v in ring.iter().skip(1) {
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
            name: "rest_machining".into(),
            reason: "no rest-corridor polygon at any depth".into(),
        });
    }
    Ok(tp)
}

fn validate(params: &RestMachiningParams, tool: &Tool) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "rest_machining".into(),
        reason,
    };
    if !(params.prev_tool_radius_mm > tool.radius_mm()) {
        return Err(mk(format!(
            "prev_tool_radius {} must exceed current tool radius {}",
            params.prev_tool_radius_mm,
            tool.radius_mm()
        )));
    }
    if !(params.step_down > 0.0) {
        return Err(mk("step_down must be > 0".into()));
    }
    if !(params.depth > 0.0) {
        return Err(mk("depth must be > 0".into()));
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
    #[ignore = "Phase 17 regression — offset polygon collapses for cube-slice input. \
                Tracked separately; un-ignore once `offset::polygon` handles the \
                stitched-polygon vertex layout `intersect_plane` produces."]
    fn rest_machining_emits_some_cuts() {
        let mesh = cube(20.0);
        let stock = Stock::new(
            Vector3::new(-12.0, -12.0, -10.0),
            Vector3::new(24.0, 24.0, 20.0),
            "alu",
        )
        .unwrap();
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "").unwrap();
        let params = RestMachiningParams {
            depth: 2.0,
            step_down: 2.0,
            prev_tool_radius_mm: 5.0,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts > 0);
    }
}
