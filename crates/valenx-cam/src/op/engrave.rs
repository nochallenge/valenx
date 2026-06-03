//! Engrave — V-bit engraving along a 2D curve.
//!
//! Plunge depth = `chord_width / 2 / tan(half_angle)`. v1 uses the
//! configured `engrave_width_mm` to size the chord and computes the
//! corresponding plunge depth from the tool's `half_angle_deg`.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Engrave parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EngraveParams {
    /// Tool id — caller is responsible for picking a V-bit tool.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Width of the engraved line (mm) — drives the plunge depth.
    pub engrave_width_mm: f64,
    /// V-bit half-angle (degrees). Common: 30°, 45°, 60°.
    pub v_bit_half_angle_deg: f64,
    /// 2D curve.
    pub curve: Vec<Vector3<f64>>,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for EngraveParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 300.0,
            plunge_feed: 100.0,
            spindle_rpm: 18_000.0,
            engrave_width_mm: 0.5,
            v_bit_half_angle_deg: 30.0,
            curve: Vec::new(),
            safe_z_clearance: 5.0,
        }
    }
}

/// Compute the engrave depth for a given chord width + half-angle.
pub fn engrave_depth(width_mm: f64, half_angle_deg: f64) -> f64 {
    let half_rad = half_angle_deg.to_radians();
    (width_mm * 0.5) / half_rad.tan()
}

/// Generate an engrave toolpath.
pub fn generate(stock: &Stock, params: &EngraveParams, _tool: &Tool) -> Result<Toolpath, CamError> {
    validate(params)?;
    let depth = engrave_depth(params.engrave_width_mm, params.v_bit_half_angle_deg);
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let z = stock.top_z() - depth;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    let start = params.curve[0];
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
    for v in params.curve.iter().skip(1) {
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(v.x, v.y, z),
            params.feed_mm_per_min,
        ));
    }
    let last = *params.curve.last().unwrap();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(last.x, last.y, safe_z),
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &EngraveParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "engrave".into(),
        reason,
    };
    if params.curve.len() < 2 {
        return Err(mk("curve needs >= 2 points".into()));
    }
    if !(params.engrave_width_mm > 0.0) {
        return Err(mk(format!(
            "engrave_width_mm must be > 0 (got {})",
            params.engrave_width_mm
        )));
    }
    if !(params.v_bit_half_angle_deg > 0.0 && params.v_bit_half_angle_deg < 90.0) {
        return Err(mk(format!(
            "v_bit_half_angle_deg must be in (0, 90) (got {})",
            params.v_bit_half_angle_deg
        )));
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
    fn depth_for_30deg_v_bit() {
        // width=1mm, half_angle=30° -> depth = 0.5 / tan(30°) = 0.5 / 0.577 = 0.866
        let d = engrave_depth(1.0, 30.0);
        assert!((d - 0.866).abs() < 0.01, "got {d}");
    }

    #[test]
    fn engrave_short_curve() {
        let stock = Stock::default();
        let tool = Tool::new(1, "V30", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let params = EngraveParams {
            curve: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(10.0, 0.0, 0.0)],
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        assert!(tp.moves.iter().any(|m| m.kind == MoveKind::Cut));
    }
}
