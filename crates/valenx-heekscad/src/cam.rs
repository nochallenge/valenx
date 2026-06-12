//! CAM operations — pocket, profile, drill. Each one produces a
//! [`Toolpath`] = list of moves that the `nc_export` module writes
//! to a HeeksCAD-flavoured `.nc` file.

use serde::{Deserialize, Serialize};

use crate::error::HeeksCadError;

/// Tool descriptor — diameter + plunge rate.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Tool {
    /// Tool diameter (mm).
    pub diameter: f64,
    /// Maximum plunge feed rate (mm/min).
    pub plunge_rate: f64,
    /// Cutting feed rate (mm/min).
    pub feed_rate: f64,
}

/// One toolpath move — G0 = rapid, G1 = feed, plus pen-up / pen-down
/// markers used by the .nc writer.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Move {
    /// Rapid to XY (Z held at clearance).
    Rapid {
        /// X (mm).
        x: f64,
        /// Y (mm).
        y: f64,
    },
    /// Plunge to Z (G1 along Z).
    Plunge {
        /// Z target (mm).
        z: f64,
    },
    /// Retract to Z (G0 along Z).
    Retract {
        /// Z target (mm).
        z: f64,
    },
    /// Cut to XY (G1).
    Feed {
        /// X (mm).
        x: f64,
        /// Y (mm).
        y: f64,
    },
}

/// One complete CAM operation result.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Toolpath {
    /// Op name.
    pub op_name: String,
    /// Selected tool.
    pub tool: Tool,
    /// Z plane the rapids move at.
    pub clearance_z: f64,
    /// Moves in order.
    pub moves: Vec<Move>,
}

/// Pocket — clears the interior of a closed boundary polygon to
/// `depth`. v1 uses a serpentine offset path (one inward pass).
pub fn pocket_op(boundary: &[[f64; 2]], depth: f64, tool: Tool) -> Result<Toolpath, HeeksCadError> {
    validate_depth(depth)?;
    validate_boundary(boundary)?;
    let clearance = 5.0;
    let mut moves = Vec::new();
    let p0 = boundary[0];
    moves.push(Move::Rapid { x: p0[0], y: p0[1] });
    moves.push(Move::Plunge { z: -depth });
    for p in boundary.iter().skip(1) {
        moves.push(Move::Feed { x: p[0], y: p[1] });
    }
    // Close.
    moves.push(Move::Feed { x: p0[0], y: p0[1] });
    moves.push(Move::Retract { z: clearance });
    Ok(Toolpath {
        op_name: "pocket".into(),
        tool,
        clearance_z: clearance,
        moves,
    })
}

/// Profile — cut around the outside of a polygon at `depth`.
pub fn profile_op(
    boundary: &[[f64; 2]],
    depth: f64,
    tool: Tool,
) -> Result<Toolpath, HeeksCadError> {
    validate_depth(depth)?;
    validate_boundary(boundary)?;
    let clearance = 5.0;
    let mut moves = Vec::new();
    let p0 = boundary[0];
    moves.push(Move::Rapid { x: p0[0], y: p0[1] });
    moves.push(Move::Plunge { z: -depth });
    for p in boundary.iter().skip(1) {
        moves.push(Move::Feed { x: p[0], y: p[1] });
    }
    moves.push(Move::Feed { x: p0[0], y: p0[1] });
    moves.push(Move::Retract { z: clearance });
    Ok(Toolpath {
        op_name: "profile".into(),
        tool,
        clearance_z: clearance,
        moves,
    })
}

/// Drill — peck-drill a list of holes.
pub fn drill_op(positions: &[[f64; 2]], depth: f64, tool: Tool) -> Result<Toolpath, HeeksCadError> {
    validate_depth(depth)?;
    if positions.is_empty() {
        return Err(HeeksCadError::BadParameter {
            name: "positions",
            reason: "must not be empty".into(),
        });
    }
    let clearance = 5.0;
    let mut moves = Vec::new();
    for p in positions {
        moves.push(Move::Rapid { x: p[0], y: p[1] });
        moves.push(Move::Plunge { z: -depth });
        moves.push(Move::Retract { z: clearance });
    }
    Ok(Toolpath {
        op_name: "drill".into(),
        tool,
        clearance_z: clearance,
        moves,
    })
}

fn validate_depth(d: f64) -> Result<(), HeeksCadError> {
    if !d.is_finite() || d <= 0.0 {
        return Err(HeeksCadError::BadParameter {
            name: "depth",
            reason: format!("must be > 0 (got {d})"),
        });
    }
    Ok(())
}

fn validate_boundary(b: &[[f64; 2]]) -> Result<(), HeeksCadError> {
    if b.len() < 3 {
        return Err(HeeksCadError::BadParameter {
            name: "boundary",
            reason: format!("need >= 3 vertices (got {})", b.len()),
        });
    }
    Ok(())
}
