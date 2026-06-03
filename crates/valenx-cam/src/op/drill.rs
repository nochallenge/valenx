//! Drill operation — peck-cycle vertical holes at given XY positions.
//!
//! ## Algorithm
//!
//! For each hole position:
//!
//! 1. Rapid up to safe-Z, rapid traverse to hole XY.
//! 2. Rapid down to retract clearance above stock top.
//! 3. Plunge down by `peck_depth`, rapid retract to retract clearance,
//!    plunge again to previous depth + `peck_depth`, retract — repeat
//!    until reaching `total_depth` below stock top.
//! 4. Final rapid up to safe-Z.

use nalgebra::Vector3;

use crate::error::CamError;
use crate::operation::DrillParams;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Generate a drill toolpath with a peck cycle at each hole position.
pub fn generate(stock: &Stock, params: &DrillParams, _tool: &Tool) -> Result<Toolpath, CamError> {
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
        // Rapid to hole XY at safe Z.
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(hole.x, hole.y, safe_z),
            0.0,
        ));
        // Rapid down to retract Z (just above stock).
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(hole.x, hole.y, retract_z),
            0.0,
        ));
        // Peck cycle.
        let mut cur_depth = 0.0;
        while cur_depth < params.total_depth - 1e-9 {
            cur_depth = (cur_depth + params.peck_depth).min(params.total_depth);
            let z = stock.top_z() - cur_depth;
            tp.push(Move::new(
                MoveKind::Plunge,
                Vector3::new(hole.x, hole.y, z),
                params.plunge_feed,
            ));
            // Retract for chip clearance.
            tp.push(Move::new(
                MoveKind::Rapid,
                Vector3::new(hole.x, hole.y, retract_z),
                0.0,
            ));
        }
        // Final rapid to safe Z so the next traversal can begin safely.
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(hole.x, hole.y, safe_z),
            0.0,
        ));
    }
    Ok(tp)
}

fn validate(params: &DrillParams) -> Result<(), CamError> {
    if !(params.peck_depth > 0.0) {
        return Err(CamError::BadOperation {
            name: "drill".into(),
            reason: format!("peck_depth must be > 0 (got {})", params.peck_depth),
        });
    }
    if !(params.total_depth > 0.0) {
        return Err(CamError::BadOperation {
            name: "drill".into(),
            reason: format!("total_depth must be > 0 (got {})", params.total_depth),
        });
    }
    if !(params.plunge_feed > 0.0) {
        return Err(CamError::BadOperation {
            name: "drill".into(),
            reason: format!("plunge_feed must be > 0 (got {})", params.plunge_feed),
        });
    }
    if params.hole_positions.is_empty() {
        return Err(CamError::BadOperation {
            name: "drill".into(),
            reason: "hole_positions is empty".into(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::{Tool, ToolKind};

    #[test]
    fn drill_three_holes_produces_three_peck_sequences() {
        let stock = Stock::default();
        let tool = Tool::new(2, "Drill5", ToolKind::Drill, 5.0, 30.0, 2, "HSS").unwrap();
        let params = DrillParams {
            tool_id: 2,
            peck_depth: 1.0,
            total_depth: 3.0,
            plunge_feed: 100.0,
            spindle_rpm: 1500.0,
            retract_clearance: 1.0,
            safe_z_clearance: 5.0,
            hole_positions: vec![
                Vector3::new(10.0, 10.0, 0.0),
                Vector3::new(20.0, 10.0, 0.0),
                Vector3::new(30.0, 10.0, 0.0),
            ],
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        // Per hole: 1 traversal + 1 down-to-retract + 3 plunges + 3 rapids + 1 final = 9 moves
        // Plus initial safe-Z move.
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert_eq!(plunges, 9, "3 holes × 3 pecks = 9 plunges");
    }

    #[test]
    fn drill_validate_rejects_empty_positions() {
        let p = DrillParams {
            hole_positions: vec![],
            ..Default::default()
        };
        assert!(validate(&p).is_err());
    }

    #[test]
    fn drill_validate_rejects_zero_peck() {
        let p = DrillParams {
            peck_depth: 0.0,
            hole_positions: vec![Vector3::new(0.0, 0.0, 0.0)],
            ..Default::default()
        };
        assert!(validate(&p).is_err());
    }

    #[test]
    fn drill_pecks_partial_final_depth() {
        // total_depth=3, peck_depth=2: pecks should go 2.0, 3.0.
        let stock = Stock::default();
        let tool = Tool::new(2, "Drill5", ToolKind::Drill, 5.0, 30.0, 2, "HSS").unwrap();
        let params = DrillParams {
            tool_id: 2,
            peck_depth: 2.0,
            total_depth: 3.0,
            plunge_feed: 100.0,
            spindle_rpm: 1500.0,
            retract_clearance: 1.0,
            safe_z_clearance: 5.0,
            hole_positions: vec![Vector3::new(10.0, 10.0, 0.0)],
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let plunges: Vec<_> = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .collect();
        assert_eq!(plunges.len(), 2);
        assert!((plunges[0].position.z - 8.0).abs() < 1e-9); // top=10, depth=2
        assert!((plunges[1].position.z - 7.0).abs() < 1e-9); // top=10, depth=3
    }
}
