//! Plunge rough — rough a pocket by discrete full-depth plunges.
//!
//! For each pre-computed XY plunge position, the cutter rapids to
//! safe Z, plunges to depth, and retracts. No pecking — a single
//! straight plunge per position. Useful for soft materials or with
//! drill-style end-mills.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Plunge-rough parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlungeRoughParams {
    /// Tool id.
    pub tool_id: u32,
    /// Plunge feed (mm/min) — used for every Z descent.
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// XY plunge positions. Z is ignored (taken from the stock top).
    pub plunge_positions: Vec<Vector3<f64>>,
    /// Plunge depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for PlungeRoughParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            plunge_feed: 150.0,
            spindle_rpm: 12_000.0,
            plunge_positions: Vec::new(),
            depth: 5.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate a plunge-rough toolpath.
pub fn generate(
    stock: &Stock,
    params: &PlungeRoughParams,
    _tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let bottom_z = stock.top_z() - params.depth;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));
    for pos in &params.plunge_positions {
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(pos.x, pos.y, safe_z),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(pos.x, pos.y, bottom_z),
            params.plunge_feed,
        ));
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(pos.x, pos.y, safe_z),
            0.0,
        ));
    }
    Ok(tp)
}

fn validate(params: &PlungeRoughParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "plunge_rough".into(),
        reason,
    };
    if params.plunge_positions.is_empty() {
        return Err(mk("plunge_positions is empty".into()));
    }
    if !(params.depth > 0.0) {
        return Err(mk(format!("depth must be > 0 (got {})", params.depth)));
    }
    if !(params.plunge_feed > 0.0) {
        return Err(mk(format!(
            "plunge_feed must be > 0 (got {})",
            params.plunge_feed
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;

    #[test]
    fn plunge_rough_three_positions() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let params = PlungeRoughParams {
            plunge_positions: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(10.0, 0.0, 0.0),
                Vector3::new(20.0, 0.0, 0.0),
            ],
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

    #[test]
    fn validate_rejects_empty_positions() {
        let bad = PlungeRoughParams::default();
        assert!(validate(&bad).is_err());
    }
}
