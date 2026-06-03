//! Thread mill — helical thread-milling toolpath inside a pre-bored
//! hole.
//!
//! Generates a single helical pass at radius
//! `(thread_diameter / 2 - tool_radius)` with `pitch_mm` per
//! revolution. v1 cuts only one revolution per thread (no multi-pass
//! depth control); production thread-milling uses multiple radial
//! step-overs for deep threads.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Thread-mill parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ThreadMillParams {
    /// Tool id — caller picks a thread mill.
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// XY centre of the hole.
    pub centre: Vector3<f64>,
    /// Finished thread (major) diameter (mm).
    pub thread_diameter: f64,
    /// Thread pitch (mm per revolution).
    pub pitch_mm: f64,
    /// Total thread depth below stock top (mm).
    pub thread_depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Segments per revolution (≥8).
    pub n_steps_per_rev: u32,
    /// `true` to climb-mill the thread (CCW seen from above for a
    /// right-hand thread cut top-down).
    pub climb: bool,
}

impl Default for ThreadMillParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 300.0,
            plunge_feed: 100.0,
            spindle_rpm: 4000.0,
            centre: Vector3::zeros(),
            thread_diameter: 8.0,
            pitch_mm: 1.25,
            thread_depth: 10.0,
            safe_z_clearance: 5.0,
            n_steps_per_rev: 36,
            climb: true,
        }
    }
}

/// Generate a thread-mill toolpath.
pub fn generate(
    stock: &Stock,
    params: &ThreadMillParams,
    tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params, tool)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let r = params.thread_diameter * 0.5 - tool.radius_mm();
    let n = params.n_steps_per_rev.max(8) as usize;
    let mut tp = Toolpath::new();
    let cx = params.centre.x;
    let cy = params.centre.y;
    // Rapid to centre at safe Z, then to start radius at safe Z.
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(cx, cy, safe_z),
        0.0,
    ));
    let start_xy = Vector3::new(cx + r, cy, safe_z);
    tp.push(Move::new(MoveKind::Rapid, start_xy, 0.0));
    tp.push(Move::new(
        MoveKind::Plunge,
        Vector3::new(cx + r, cy, stock.top_z()),
        params.plunge_feed,
    ));
    // Helical descent: total revolutions = depth / pitch.
    let revs = params.thread_depth / params.pitch_mm;
    let total_steps = (revs * n as f64).ceil() as usize;
    let sign = if params.climb { 1.0 } else { -1.0 };
    let two_pi = std::f64::consts::TAU;
    for step in 1..=total_steps {
        let frac = step as f64 / total_steps as f64;
        let theta = sign * two_pi * revs * frac;
        let z = stock.top_z() - params.thread_depth * frac;
        let x = cx + r * theta.cos();
        let y = cy + r * theta.sin();
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(x, y, z),
            params.feed_mm_per_min,
        ));
    }
    // Retract centre then safe Z.
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(cx, cy, safe_z),
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &ThreadMillParams, tool: &Tool) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "thread_mill".into(),
        reason,
    };
    if !(params.thread_diameter > tool.diameter_mm) {
        return Err(mk(format!(
            "thread_diameter {} must exceed tool diameter {}",
            params.thread_diameter, tool.diameter_mm
        )));
    }
    if !(params.pitch_mm > 0.0) {
        return Err(mk("pitch_mm must be > 0".into()));
    }
    if !(params.thread_depth > 0.0) {
        return Err(mk("thread_depth must be > 0".into()));
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(mk("feed must be > 0".into()));
    }
    if params.n_steps_per_rev < 8 {
        return Err(mk("n_steps_per_rev must be >= 8".into()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;

    #[test]
    fn thread_mill_helix() {
        let stock = Stock::default();
        let tool = Tool::new(1, "TM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let params = ThreadMillParams {
            thread_diameter: 10.0,
            pitch_mm: 1.5,
            thread_depth: 6.0,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts > 50, "expected many helix cuts, got {cuts}");
    }
}
