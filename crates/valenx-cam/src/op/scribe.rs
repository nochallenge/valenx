//! Scribe — single-pass shallow line cut along a 2D curve.
//!
//! Differs from [`crate::op::engrave`] in that the depth is a fixed
//! constant (not derived from V-bit geometry). One pass at full
//! cut feed.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Scribe parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScribeParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Curve to scribe.
    pub curve: Vec<Vector3<f64>>,
    /// Depth below stock top (mm) — typically 0.05–0.2 mm.
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for ScribeParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 400.0,
            plunge_feed: 100.0,
            spindle_rpm: 12_000.0,
            curve: Vec::new(),
            depth: 0.1,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a scribe toolpath.
pub fn generate(stock: &Stock, params: &ScribeParams, _tool: &Tool) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let z = stock.top_z() - params.depth;
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

fn validate(params: &ScribeParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "scribe".into(),
        reason,
    };
    if params.curve.len() < 2 {
        return Err(mk("curve needs >= 2 points".into()));
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

    #[test]
    fn scribe_short_line() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM1", ToolKind::EndMill, 1.0, 10.0, 2, "").unwrap();
        let params = ScribeParams {
            curve: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0)],
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert_eq!(cuts, 1);
    }
}
