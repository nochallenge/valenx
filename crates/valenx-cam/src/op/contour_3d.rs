//! Contour 3D — follow a 3D edge curve at its actual Z coordinates.
//!
//! Unlike [`crate::op::contour_2d`], the curve is followed in 3D
//! (the per-vertex Z is honoured). Used for ridge / fillet finishing
//! and rest-machining.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Contour 3D parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Contour3DParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// 3D source curve.
    pub curve: Vec<Vector3<f64>>,
    /// `true` to close the curve.
    pub closed: bool,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for Contour3DParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 600.0,
            plunge_feed: 200.0,
            spindle_rpm: 12_000.0,
            curve: Vec::new(),
            closed: false,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a 3D-contour toolpath.
pub fn generate(
    stock: &Stock,
    params: &Contour3DParams,
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
    let start = params.curve[0];
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(start.x, start.y, safe_z),
        0.0,
    ));
    tp.push(Move::new(MoveKind::Plunge, start, params.plunge_feed));
    for v in params.curve.iter().skip(1) {
        tp.push(Move::new(MoveKind::Cut, *v, params.feed_mm_per_min));
    }
    if params.closed {
        tp.push(Move::new(MoveKind::Cut, start, params.feed_mm_per_min));
    }
    let last = *params.curve.last().unwrap();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(last.x, last.y, safe_z),
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &Contour3DParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "contour_3d".into(),
        reason,
    };
    if params.curve.len() < 2 {
        return Err(mk(format!(
            "curve needs >= 2 points (got {})",
            params.curve.len()
        )));
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(mk(format!(
            "feed must be > 0 (got {})",
            params.feed_mm_per_min
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;

    #[test]
    fn contour_3d_descending_curve() {
        let stock = Stock::default();
        let tool = Tool::new(1, "BM6", ToolKind::BallMill, 6.0, 25.0, 2, "").unwrap();
        let params = Contour3DParams {
            curve: vec![
                Vector3::new(0.0, 0.0, 5.0),
                Vector3::new(10.0, 0.0, 0.0),
                Vector3::new(20.0, 0.0, -5.0),
            ],
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts >= 2);
    }
}
