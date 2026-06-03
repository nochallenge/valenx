//! True 5-axis continuous tool-centre-point contour.
//!
//! v1 simplification: the path is supplied as a polyline + per-point
//! tool-orientation vector (unit vector pointing along the tool's
//! spindle axis). The op computes A/B axis values that align the
//! spindle with each orientation vector.
//!
//! A robotics-grade implementation would solve for redundant axes,
//! sing-cone avoidance, and TCP-mode kinematic compensation; those
//! are out of scope for v1.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::axis::{Move5Ax, Toolpath5Ax};
use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::MoveKind;

/// Parameters for a 5-axis TCP contour.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tcp5AxContourParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// 3D path points.
    pub path: Vec<Vector3<f64>>,
    /// Per-point tool axis directions (must equal `path.len()`).
    /// Unit vectors pointing from the part surface up into the tool.
    pub tool_axes: Vec<Vector3<f64>>,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for Tcp5AxContourParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 500.0,
            plunge_feed: 200.0,
            spindle_rpm: 12_000.0,
            path: Vec::new(),
            tool_axes: Vec::new(),
            safe_z_clearance: 5.0,
        }
    }
}

/// Convert a unit tool axis direction to A/B (degrees). Convention:
///
/// - `A` rotates the X axis toward Y (around X).
/// - `B` rotates Z toward X (around Y).
///
/// The straight-down axis (0, 0, 1) maps to A=0, B=0.
pub fn axis_to_ab(axis: Vector3<f64>) -> (f64, f64) {
    let n = axis.norm();
    if n < 1e-9 {
        return (0.0, 0.0);
    }
    let u = axis / n;
    // B = arctan2(ux, uz) — tilt about Y.
    let b = u.x.atan2(u.z);
    // A = arcsin(-uy) — tilt about X (project onto YZ first).
    let denom = (u.x * u.x + u.z * u.z).sqrt().max(1e-12);
    let a = (-u.y / denom).atan2(1.0);
    (a.to_degrees(), b.to_degrees())
}

/// Generate a 5-axis TCP contour toolpath.
pub fn generate(
    stock: &Stock,
    params: &Tcp5AxContourParams,
    _tool: &Tool,
) -> Result<Toolpath5Ax, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath5Ax::new();
    // Initial rapid to safe Z above first point.
    let first = params.path[0];
    tp.push(Move5Ax::new(
        MoveKind::Rapid,
        Vector3::new(first.x, first.y, safe_z),
        0.0,
        0.0,
        0.0,
    ));
    let (a0, b0) = axis_to_ab(params.tool_axes[0]);
    tp.push(Move5Ax::new(
        MoveKind::Plunge,
        first,
        a0,
        b0,
        params.plunge_feed,
    ));
    for (p, ax) in params.path.iter().zip(params.tool_axes.iter()).skip(1) {
        let (a, b) = axis_to_ab(*ax);
        tp.push(Move5Ax::new(
            MoveKind::Cut,
            *p,
            a,
            b,
            params.feed_mm_per_min,
        ));
    }
    // Retract.
    let last = *params.path.last().unwrap();
    tp.push(Move5Ax::new(
        MoveKind::Rapid,
        Vector3::new(last.x, last.y, safe_z),
        0.0,
        0.0,
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &Tcp5AxContourParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "tcp_5ax_contour".into(),
        reason,
    };
    if params.path.len() < 2 {
        return Err(mk("path needs >= 2 points".into()));
    }
    if params.tool_axes.len() != params.path.len() {
        return Err(mk(format!(
            "tool_axes length {} ≠ path length {}",
            params.tool_axes.len(),
            params.path.len()
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
    fn straight_down_axis_maps_to_zero_ab() {
        let (a, b) = axis_to_ab(Vector3::new(0.0, 0.0, 1.0));
        assert!(a.abs() < 1e-6, "got A={a}");
        assert!(b.abs() < 1e-6, "got B={b}");
    }

    #[test]
    fn tilted_path_maps_to_nonzero_ab() {
        let stock = Stock::default();
        let tool = Tool::new(1, "BM6", ToolKind::BallMill, 6.0, 25.0, 2, "").unwrap();
        let params = Tcp5AxContourParams {
            path: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(10.0, 0.0, 0.0)],
            tool_axes: vec![
                Vector3::new(1.0, 0.0, 1.0).normalize(),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        assert!(tp.len() >= 3);
        // First plunge should be at non-zero B.
        let plunge = tp
            .moves
            .iter()
            .find(|m| m.kind == MoveKind::Plunge)
            .unwrap();
        assert!(plunge.b_deg.abs() > 1e-6);
    }
}
