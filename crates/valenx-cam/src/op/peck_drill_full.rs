//! Peck-drill full — extended drill peck cycle with dwell and chip
//! clearance retract.
//!
//! Differs from [`crate::op::drill`] by adding a configurable dwell
//! at the bottom of each peck (to break the chip) and an explicit
//! full retract every N pecks for chip clearance.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Peck-drill-full parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PeckDrillFullParams {
    /// Tool id (Drill).
    pub tool_id: u32,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Per-peck depth (mm).
    pub peck_depth: f64,
    /// Total depth below stock top (mm).
    pub total_depth: f64,
    /// Retract Z clearance (mm above stock top) between pecks.
    pub retract_clearance: f64,
    /// Safe-Z clearance for inter-hole rapids (mm).
    pub safe_z_clearance: f64,
    /// Hole XY positions (Z ignored).
    pub hole_positions: Vec<Vector3<f64>>,
    /// Dwell time at bottom of each peck (seconds).
    pub dwell_at_bottom_s: f64,
    /// Number of pecks between full retracts (0 = no full retracts).
    pub full_retract_every_n_pecks: u32,
}

impl Default for PeckDrillFullParams {
    fn default() -> Self {
        Self {
            tool_id: 2,
            plunge_feed: 100.0,
            spindle_rpm: 1500.0,
            peck_depth: 1.0,
            total_depth: 5.0,
            retract_clearance: 1.0,
            safe_z_clearance: 5.0,
            hole_positions: Vec::new(),
            dwell_at_bottom_s: 0.1,
            full_retract_every_n_pecks: 3,
        }
    }
}

/// Generate a peck-drill-full toolpath. The dwell at the bottom of
/// each peck is recorded as a zero-distance `Cut` move at feed 0;
/// postprocessors can recognise this as a G4 dwell. For v1 we encode
/// the dwell as `feed = -seconds` (negative feed = dwell signal).
/// This stays compatible with simulation (zero-length segment) and
/// the Fanuc-family postprocessors can detect the marker.
pub fn generate(
    stock: &Stock,
    params: &PeckDrillFullParams,
    _tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let retract_z = stock.top_z() + params.retract_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    for hole in &params.hole_positions {
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(hole.x, hole.y, safe_z),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(hole.x, hole.y, retract_z),
            0.0,
        ));
        let mut depth = 0.0_f64;
        let mut peck_idx = 0_u32;
        while depth < params.total_depth - 1e-9 {
            depth = (depth + params.peck_depth).min(params.total_depth);
            let z = stock.top_z() - depth;
            tp.push(Move::new(
                MoveKind::Plunge,
                Vector3::new(hole.x, hole.y, z),
                params.plunge_feed,
            ));
            if params.dwell_at_bottom_s > 0.0 {
                // Encode dwell as a zero-length "cut" with negative feed.
                tp.push(Move::new(
                    MoveKind::Cut,
                    Vector3::new(hole.x, hole.y, z),
                    -params.dwell_at_bottom_s,
                ));
            }
            peck_idx += 1;
            if params.full_retract_every_n_pecks > 0
                && peck_idx % params.full_retract_every_n_pecks == 0
            {
                tp.push(Move::new(
                    MoveKind::Rapid,
                    Vector3::new(hole.x, hole.y, safe_z),
                    0.0,
                ));
                tp.push(Move::new(
                    MoveKind::Rapid,
                    Vector3::new(hole.x, hole.y, retract_z),
                    0.0,
                ));
            } else {
                tp.push(Move::new(
                    MoveKind::Rapid,
                    Vector3::new(hole.x, hole.y, retract_z),
                    0.0,
                ));
            }
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(hole.x, hole.y, safe_z),
            0.0,
        ));
    }
    Ok(tp)
}

fn validate(params: &PeckDrillFullParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "peck_drill_full".into(),
        reason,
    };
    if !(params.peck_depth > 0.0) {
        return Err(mk("peck_depth must be > 0".into()));
    }
    if !(params.total_depth > 0.0) {
        return Err(mk("total_depth must be > 0".into()));
    }
    if !(params.plunge_feed > 0.0) {
        return Err(mk("plunge_feed must be > 0".into()));
    }
    if params.hole_positions.is_empty() {
        return Err(mk("hole_positions is empty".into()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;

    #[test]
    fn peck_drill_full_emits_pecks() {
        let stock = Stock::default();
        let tool = Tool::new(2, "D5", ToolKind::Drill, 5.0, 30.0, 2, "HSS").unwrap();
        let params = PeckDrillFullParams {
            peck_depth: 1.0,
            total_depth: 3.0,
            hole_positions: vec![Vector3::new(0.0, 0.0, 0.0)],
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
