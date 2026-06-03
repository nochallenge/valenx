//! Slot — straight slot cut.
//!
//! Plunges to depth at one end of the slot, cuts to the other end,
//! then retracts. v1 emits a simple back-and-forth pass at each Z
//! step-down (no widening — tool diameter equals the slot width).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Slot parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SlotParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Slot start XY.
    pub start: Vector3<f64>,
    /// Slot end XY.
    pub end: Vector3<f64>,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for SlotParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 600.0,
            plunge_feed: 200.0,
            spindle_rpm: 12_000.0,
            start: Vector3::zeros(),
            end: Vector3::new(20.0, 0.0, 0.0),
            step_down: 1.0,
            depth: 3.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a slot toolpath.
pub fn generate(stock: &Stock, params: &SlotParams, _tool: &Tool) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "slot")?;
    let mut from = params.start;
    let mut to = params.end;
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(from.x, from.y, safe_z),
        0.0,
    ));
    for k in 1..=n_passes {
        let z = stock.top_z() - (params.step_down * k as f64).min(params.depth);
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(from.x, from.y, z),
            params.plunge_feed,
        ));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(to.x, to.y, z),
            params.feed_mm_per_min,
        ));
        // Reverse for the next pass so we don't rapid back.
        std::mem::swap(&mut from, &mut to);
    }
    let end_xy = from;
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(end_xy.x, end_xy.y, safe_z),
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &SlotParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "slot".into(),
        reason,
    };
    let len = (params.end - params.start).norm();
    if !(len > 0.0) {
        return Err(mk("slot start and end coincide".into()));
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

    #[test]
    fn slot_emits_one_plunge_per_pass() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let params = SlotParams {
            start: Vector3::new(0.0, 0.0, 0.0),
            end: Vector3::new(30.0, 0.0, 0.0),
            step_down: 1.0,
            depth: 3.0,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert_eq!(plunges, 3);
    }
}
