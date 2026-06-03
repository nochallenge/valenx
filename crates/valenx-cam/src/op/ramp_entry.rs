//! Ramp entry — angled entry into material instead of a vertical
//! plunge.
//!
//! Generates a back-and-forth ramp at `ramp_angle_deg` along a
//! supplied 2D entry segment, descending by the configured `depth`.
//! Used to safely engage materials that cannot tolerate full plunges
//! (e.g. aluminum with non-centre-cutting end-mills).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Ramp-entry parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RampEntryParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Ramp segment start XY (Z is ignored).
    pub ramp_start: Vector3<f64>,
    /// Ramp segment end XY (Z is ignored).
    pub ramp_end: Vector3<f64>,
    /// Ramp angle in degrees (1°–10° typical).
    pub ramp_angle_deg: f64,
    /// Total depth below stock top to reach (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for RampEntryParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 600.0,
            spindle_rpm: 12_000.0,
            ramp_start: Vector3::zeros(),
            ramp_end: Vector3::new(10.0, 0.0, 0.0),
            ramp_angle_deg: 3.0,
            depth: 3.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a ramp-entry toolpath.
pub fn generate(
    stock: &Stock,
    params: &RampEntryParams,
    _tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let dir = Vector3::new(
        params.ramp_end.x - params.ramp_start.x,
        params.ramp_end.y - params.ramp_start.y,
        0.0,
    );
    let len = dir.norm();
    if !(len > 0.0) {
        return Err(CamError::BadOperation {
            name: "ramp_entry".into(),
            reason: "ramp segment has zero length".into(),
        });
    }
    let dz_per_pass = len * params.ramp_angle_deg.to_radians().tan();
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(params.ramp_start.x, params.ramp_start.y, safe_z),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(params.ramp_start.x, params.ramp_start.y, stock.top_z()),
        0.0,
    ));
    let mut current_z = stock.top_z();
    let target_z = stock.top_z() - params.depth;
    let mut at_start = true;
    while current_z > target_z + 1e-9 {
        let new_z = (current_z - dz_per_pass).max(target_z);
        let dest = if at_start {
            params.ramp_end
        } else {
            params.ramp_start
        };
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(dest.x, dest.y, new_z),
            params.feed_mm_per_min,
        ));
        current_z = new_z;
        at_start = !at_start;
    }
    let last_xy = if at_start {
        params.ramp_start
    } else {
        params.ramp_end
    };
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(last_xy.x, last_xy.y, safe_z),
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &RampEntryParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "ramp_entry".into(),
        reason,
    };
    if !(params.ramp_angle_deg > 0.0 && params.ramp_angle_deg <= 45.0) {
        return Err(mk(format!(
            "ramp_angle_deg must be in (0, 45] (got {})",
            params.ramp_angle_deg
        )));
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
    fn ramp_entry_descends() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let params = RampEntryParams {
            ramp_start: Vector3::new(0.0, 0.0, 0.0),
            ramp_end: Vector3::new(20.0, 0.0, 0.0),
            ramp_angle_deg: 5.0,
            depth: 3.0,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts >= 1);
    }
}
