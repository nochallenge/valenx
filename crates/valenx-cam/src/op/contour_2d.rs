//! Contour 2D — follow a 2D curve (XY polyline) at constant Z.
//!
//! Generalises [`crate::op::profile`] for cases where the boundary
//! comes from a user-supplied curve rather than being extracted from
//! a source mesh. Used by the Engrave / Scribe / ThreadMill ops as a
//! lower-level primitive.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Contour 2D parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Contour2DParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Curve points (XY only — Z values are overridden per pass).
    pub curve: Vec<Vector3<f64>>,
    /// `true` to interpret the curve as closed (cut joins last back
    /// to first).
    pub closed: bool,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for Contour2DParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 600.0,
            plunge_feed: 200.0,
            spindle_rpm: 12_000.0,
            curve: Vec::new(),
            closed: false,
            step_down: 1.0,
            depth: 2.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a contour-2D toolpath following `curve` at each Z step.
pub fn generate(
    stock: &Stock,
    params: &Contour2DParams,
    _tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "contour_2d")?;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    for k in 1..=n_passes {
        let depth_below_top = (params.step_down * k as f64).min(params.depth);
        let z = stock.top_z() - depth_below_top;
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
        if params.closed {
            tp.push(Move::new(
                MoveKind::Cut,
                Vector3::new(start.x, start.y, z),
                params.feed_mm_per_min,
            ));
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(start.x, start.y, safe_z),
            0.0,
        ));
    }
    Ok(tp)
}

fn validate(params: &Contour2DParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "contour_2d".into(),
        reason,
    };
    if params.curve.len() < 2 {
        return Err(mk(format!(
            "curve needs >= 2 points (got {})",
            params.curve.len()
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
    fn contour_2d_open_curve() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let params = Contour2DParams {
            curve: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(10.0, 0.0, 0.0),
                Vector3::new(10.0, 10.0, 0.0),
            ],
            step_down: 1.0,
            depth: 2.0,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert_eq!(plunges, 2, "2 passes -> 2 plunges");
    }
}
