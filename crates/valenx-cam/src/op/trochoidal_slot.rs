//! Trochoidal slot — slot widening via trochoidal (loop) motion.
//!
//! Used to cut slots wider than the tool diameter while keeping
//! constant tool engagement. The tool follows the slot centreline,
//! superimposing circular loops with `loop_radius`. Each loop
//! advances the centre by `step_over` along the centreline.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Trochoidal-slot parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TrochoidalSlotParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Slot centreline endpoints (XY).
    pub centreline: Vec<Vector3<f64>>,
    /// Finished slot width (mm) — must be > tool diameter.
    pub slot_width: f64,
    /// Trochoidal loop radius (mm) — typically (slot_width - tool_diameter) / 2.
    pub loop_radius: f64,
    /// Step-over per loop along the centreline (mm).
    pub step_over: f64,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Segments per trochoidal loop (≥8).
    pub n_steps_per_loop: u32,
}

impl Default for TrochoidalSlotParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 1500.0,
            plunge_feed: 300.0,
            spindle_rpm: 18_000.0,
            centreline: Vec::new(),
            slot_width: 8.0,
            loop_radius: 1.0,
            step_over: 0.5,
            step_down: 4.0,
            depth: 4.0,
            safe_z_clearance: 5.0,
            n_steps_per_loop: 16,
        }
    }
}

/// Generate a trochoidal-slot toolpath.
pub fn generate(
    stock: &Stock,
    params: &TrochoidalSlotParams,
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
    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "trochoidal_slot")?;
    let n = params.n_steps_per_loop.max(8) as usize;
    let two_pi = std::f64::consts::TAU;
    for k in 1..=n_passes {
        let z = stock.top_z() - (params.step_down * k as f64).min(params.depth);
        // Walk the centreline.
        let segments = params.centreline.windows(2);
        let mut first = true;
        for win in segments {
            let p0 = win[0];
            let p1 = win[1];
            let dir = Vector3::new(p1.x - p0.x, p1.y - p0.y, 0.0);
            let total = dir.norm();
            if !(total > 0.0) {
                continue;
            }
            let dir_u = dir / total;
            let perp = Vector3::new(-dir_u.y, dir_u.x, 0.0);
            let n_loops = ((total / params.step_over).floor() as usize).max(1);
            for i in 0..=n_loops {
                let t = (i as f64) * params.step_over;
                let t = t.min(total);
                let centre = Vector3::new(p0.x + dir_u.x * t, p0.y + dir_u.y * t, z);
                if first {
                    // Approach centre at safe Z, then plunge.
                    tp.push(Move::new(
                        MoveKind::Rapid,
                        Vector3::new(
                            centre.x + params.loop_radius * perp.x,
                            centre.y + params.loop_radius * perp.y,
                            safe_z,
                        ),
                        0.0,
                    ));
                    tp.push(Move::new(
                        MoveKind::Plunge,
                        Vector3::new(
                            centre.x + params.loop_radius * perp.x,
                            centre.y + params.loop_radius * perp.y,
                            z,
                        ),
                        params.plunge_feed,
                    ));
                    first = false;
                }
                for s in 0..=n {
                    let theta = two_pi * (s as f64 / n as f64);
                    let x = centre.x
                        + params.loop_radius * (perp.x * theta.cos() + dir_u.x * theta.sin());
                    let y = centre.y
                        + params.loop_radius * (perp.y * theta.cos() + dir_u.y * theta.sin());
                    tp.push(Move::new(
                        MoveKind::Cut,
                        Vector3::new(x, y, z),
                        params.feed_mm_per_min,
                    ));
                }
            }
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(0.0, 0.0, safe_z),
            0.0,
        ));
    }
    Ok(tp)
}

fn validate(params: &TrochoidalSlotParams, tool: &Tool) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "trochoidal_slot".into(),
        reason,
    };
    if params.centreline.len() < 2 {
        return Err(mk("centreline needs >= 2 points".into()));
    }
    if !(params.slot_width > tool.diameter_mm) {
        return Err(mk(format!(
            "slot_width {} must exceed tool diameter {}",
            params.slot_width, tool.diameter_mm
        )));
    }
    if !(params.loop_radius > 0.0) {
        return Err(mk(format!(
            "loop_radius must be > 0 (got {})",
            params.loop_radius
        )));
    }
    if !(params.step_over > 0.0) {
        return Err(mk("step_over must be > 0".into()));
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

    #[test]
    fn trochoidal_slot_short_segment() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM4", ToolKind::EndMill, 4.0, 25.0, 2, "").unwrap();
        let params = TrochoidalSlotParams {
            centreline: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(20.0, 0.0, 0.0)],
            slot_width: 8.0,
            loop_radius: 1.0,
            step_over: 1.0,
            step_down: 4.0,
            depth: 4.0,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts > 20, "expected many trochoidal cuts, got {cuts}");
    }
}
