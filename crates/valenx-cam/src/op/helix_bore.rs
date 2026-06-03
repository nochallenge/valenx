//! Helix bore — helical descent for boring large holes.
//!
//! Generates a spiral-down toolpath that descends `pitch` per
//! revolution while orbiting at radius `(bore_radius - tool_radius)`.
//! For each revolution the cutter advances `pitch` mm downward; the
//! orbit is approximated as a polyline with `n_steps_per_rev` segments.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Helix bore parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HelicalBoreParams {
    /// Tool id (typically EndMill).
    pub tool_id: u32,
    /// Cut feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min) — used for the initial Z descent.
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// XY centre of the bore.
    pub centre: Vector3<f64>,
    /// Finished bore radius (mm) — must be > tool radius.
    pub bore_radius: f64,
    /// Per-revolution downward pitch (mm).
    pub pitch: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Polyline segments per orbit (≥3).
    pub n_steps_per_rev: u32,
}

impl Default for HelicalBoreParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 500.0,
            plunge_feed: 200.0,
            spindle_rpm: 12_000.0,
            centre: Vector3::zeros(),
            bore_radius: 5.0,
            pitch: 0.5,
            depth: 10.0,
            safe_z_clearance: 5.0,
            n_steps_per_rev: 32,
        }
    }
}

/// Generate a helical-bore toolpath.
pub fn generate(
    stock: &Stock,
    params: &HelicalBoreParams,
    tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params, tool)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let r = params.bore_radius - tool.radius_mm();
    let n = params.n_steps_per_rev.max(3) as usize;
    let mut tp = Toolpath::new();
    let cx = params.centre.x;
    let cy = params.centre.y;
    // Rapid to safe Z above centre, then plunge to top of stock at orbit start.
    let start_xy = Vector3::new(cx + r, cy, safe_z);
    tp.push(Move::new(MoveKind::Rapid, start_xy, 0.0));
    tp.push(Move::new(
        MoveKind::Plunge,
        Vector3::new(cx + r, cy, stock.top_z()),
        params.plunge_feed,
    ));
    // Spiral down by `pitch` per revolution until reaching depth.
    //
    // Bounded loop: once `z` reaches `z_end` (clamped by `.max`)
    // the previous `while z > z_end - 1e-9` form would never exit
    // because z stays exactly at z_end forever and that bound is
    // always true. We instead compute the exact step count up front
    // and iterate that many times — `n` steps per revolution times
    // `depth / pitch` revolutions, plus a final partial-revolution
    // remainder so the spiral lands cleanly on z_end.
    let two_pi = std::f64::consts::TAU;
    let dz_per_step = -params.pitch / n as f64;
    let z_end = stock.top_z() - params.depth;
    // Total descent / per-step descent, rounded up so we always
    // reach (or just-past then clamp to) z_end.
    let n_steps = ((params.depth / params.pitch) * n as f64).ceil() as usize;
    let mut z = stock.top_z();
    for step in 1..=n_steps {
        let theta = (step as f64) * two_pi / n as f64;
        z = (z + dz_per_step).max(z_end);
        let x = cx + r * theta.cos();
        let y = cy + r * theta.sin();
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(x, y, z),
            params.feed_mm_per_min,
        ));
        if (z - z_end).abs() < 1e-9 {
            break;
        }
    }
    // Final flat orbit at the bottom to clean up.
    for k in 0..=n {
        let theta = k as f64 * two_pi / n as f64;
        let x = cx + r * theta.cos();
        let y = cy + r * theta.sin();
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(x, y, z_end),
            params.feed_mm_per_min,
        ));
    }
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(cx, cy, safe_z),
        0.0,
    ));
    Ok(tp)
}

fn validate(params: &HelicalBoreParams, tool: &Tool) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "helix_bore".into(),
        reason,
    };
    if !(params.bore_radius > tool.radius_mm()) {
        return Err(mk(format!(
            "bore_radius {} must exceed tool radius {}",
            params.bore_radius,
            tool.radius_mm()
        )));
    }
    if !(params.pitch > 0.0) {
        return Err(mk(format!("pitch must be > 0 (got {})", params.pitch)));
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
    if params.n_steps_per_rev < 3 {
        return Err(mk(format!(
            "n_steps_per_rev must be >= 3 (got {})",
            params.n_steps_per_rev
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;

    #[test]
    fn helix_bore_descends() {
        let stock = Stock::default();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let params = HelicalBoreParams {
            bore_radius: 10.0,
            pitch: 1.0,
            depth: 4.0,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts > 4, "expected helical descent, got {cuts} cuts");
    }

    #[test]
    fn validate_rejects_small_bore() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let bad = HelicalBoreParams {
            bore_radius: 1.0,
            ..Default::default()
        };
        assert!(validate(&bad, &tool).is_err());
    }
}
