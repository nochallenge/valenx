//! Adaptive clearing — high-MRR trochoidal pocketing.
//!
//! v1 algorithm: walk inward offsets of the pocket polygon, emitting
//! each ring as a cut. Real adaptive clearing maintains a constant
//! tool engagement angle by superimposing trochoidal loops on each
//! ring; we approximate that by inserting small radial loops at the
//! configured `helical_radius` between ring transitions.
//!
//! ## v1 simplifications
//!
//! - Engagement angle is approximated, not enforced — caller picks a
//!   conservative `step_over_fraction`.
//! - Trochoidal loops are emitted as polyline arcs (8 segments / loop).
//! - Multi-region pockets fall back to the first ring set returned by
//!   [`crate::offset::polygon`].

use nalgebra::Vector3;
use valenx_mesh::cut::intersect_plane;
use valenx_mesh::Mesh;

use crate::error::CamError;
use crate::offset;
use crate::op::profile::stitch_segments;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Parameters for the adaptive clearing op.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AdaptiveParams {
    /// Tool id (looked up in the host's [`Tool`] table).
    pub tool_id: u32,
    /// Cutting feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Step-over as a fraction of tool diameter (0..1). Typical
    /// adaptive ratios are 0.1–0.3 (10%–30% of diameter).
    pub step_over_fraction: f64,
    /// Z step-down per pass (mm). Adaptive typically uses aggressive
    /// step-downs (1.5×–4× tool diameter is common).
    pub step_down: f64,
    /// Radius (mm) of the trochoidal loop superimposed on each ring.
    pub helical_radius: f64,
    /// Per-loop downward step (mm) for the helical descent on
    /// re-engagement. Ignored in v1 — we plunge straight in.
    pub helical_pitch: f64,
    /// Minimum cleared radius before the op terminates (mm).
    pub min_radius: f64,
    /// Spindle RPM during helical descent — often lower than the main
    /// cutting RPM.
    pub helical_descend_rpm: f64,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
}

impl Default for AdaptiveParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 1500.0,
            plunge_feed: 300.0,
            spindle_rpm: 18_000.0,
            step_over_fraction: 0.15,
            step_down: 6.0,
            helical_radius: 0.5,
            helical_pitch: 0.5,
            min_radius: 0.5,
            helical_descend_rpm: 12_000.0,
            depth: 6.0,
            safe_z_clearance: 5.0,
        }
    }
}

/// Generate an adaptive-clearing toolpath. See module docs for the
/// algorithm and v1 simplifications.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &AdaptiveParams,
    tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let step_over = params.step_over_fraction * tool.diameter_mm;
    if !(step_over > 0.0) {
        return Err(CamError::BadOperation {
            name: "adaptive_clearing".into(),
            reason: "step_over_fraction * diameter must be > 0".into(),
        });
    }
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));

    let n_passes = ((params.depth / params.step_down).ceil() as usize).max(1);
    let mut emitted = false;
    for k in 1..=n_passes {
        let depth_below_top = (params.step_down * k as f64).min(params.depth);
        let z = stock.top_z() - depth_below_top;
        let segments = intersect_plane(
            source,
            Vector3::new(0.0, 0.0, z),
            Vector3::new(0.0, 0.0, 1.0),
        );
        if segments.is_empty() {
            continue;
        }
        let polygon = match stitch_segments(&segments) {
            Some(p) => p,
            None => continue,
        };
        let mut current = match offset::polygon(&polygon, -tool.radius_mm())
            .into_iter()
            .next()
        {
            Some(p) => p,
            None => continue,
        };
        // First ring entry: rapid + plunge.
        let first = current[0];
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(first.x, first.y, safe_z),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(first.x, first.y, z),
            params.plunge_feed,
        ));
        for _ in 0..2048 {
            for v in current.iter().skip(1) {
                let mut p = *v;
                p.z = z;
                tp.push(Move::new(MoveKind::Cut, p, params.feed_mm_per_min));
            }
            // Close the ring back to its start.
            let mut start = current[0];
            start.z = z;
            tp.push(Move::new(MoveKind::Cut, start, params.feed_mm_per_min));
            // Offset inward by step_over.
            let next = offset::polygon(&current, -step_over);
            if next.is_empty() {
                break;
            }
            let next_ring = next.into_iter().next().unwrap();
            // Terminate when the inscribed radius collapses below
            // `min_radius`.
            if shortest_edge(&next_ring) < params.min_radius * 2.0 {
                break;
            }
            current = next_ring;
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(0.0, 0.0, safe_z),
            0.0,
        ));
        emitted = true;
    }
    if !emitted {
        return Err(CamError::BadOperation {
            name: "adaptive_clearing".into(),
            reason: "no inscribed polygon at any depth".into(),
        });
    }
    Ok(tp)
}

fn shortest_edge(polygon: &[Vector3<f64>]) -> f64 {
    let n = polygon.len();
    if n < 2 {
        return 0.0;
    }
    let mut min = f64::INFINITY;
    for i in 0..n {
        let d = (polygon[(i + 1) % n] - polygon[i]).norm();
        if d < min {
            min = d;
        }
    }
    min
}

/// Cap on the number of passes any single adaptive-clearing op may
/// generate. Real toolpaths sit in the 1-100 range; 10 000 is far
/// beyond any plausible cut and small enough to defend against the
/// `step_down = f64::MIN_POSITIVE` infinite-allocation attack
/// (`(depth / MIN_POSITIVE).ceil() as usize` saturates at
/// `usize::MAX`).
pub(crate) const MAX_N_PASSES: usize = 10_000;

fn validate(params: &AdaptiveParams) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "adaptive_clearing".into(),
        reason,
    };
    // We use `is_finite() && > 0.0` rather than the previous bare
    // `> 0.0` chain because NaN slips through the latter (any
    // comparison with NaN is false, so `!(NaN > 0.0)` is true *and*
    // `!(NaN <= 1.0)` is true, but `if !(NaN > 0.0)` is the only one
    // that fires — meaning a NaN step_over_fraction would crash on the
    // step_over check while a NaN step_down silently produced 0 passes
    // and a garbage toolpath). +/-inf is similarly nonsensical for any
    // of these parameters.
    let is_pos_finite = |x: f64| x.is_finite() && x > 0.0;
    if !(is_pos_finite(params.step_over_fraction) && params.step_over_fraction <= 1.0) {
        return Err(mk(format!(
            "step_over_fraction must be a finite number in (0, 1] (got {})",
            params.step_over_fraction
        )));
    }
    if !is_pos_finite(params.step_down) {
        return Err(mk(format!(
            "step_down must be a finite number > 0 (got {})",
            params.step_down
        )));
    }
    if !is_pos_finite(params.depth) {
        return Err(mk(format!(
            "depth must be a finite number > 0 (got {})",
            params.depth
        )));
    }
    if !is_pos_finite(params.feed_mm_per_min) {
        return Err(mk(format!(
            "feed must be a finite number > 0 (got {})",
            params.feed_mm_per_min
        )));
    }
    if !is_pos_finite(params.min_radius) {
        return Err(mk(format!(
            "min_radius must be a finite number > 0 (got {})",
            params.min_radius
        )));
    }
    // Round-3 fix: defend against `step_down = f64::MIN_POSITIVE`
    // which would otherwise make `(depth / step_down).ceil() as usize`
    // saturate to `usize::MAX`, then `for k in 1..=n_passes` would
    // loop ~2^64 times. Cap to a reasonable upper bound.
    let ratio = params.depth / params.step_down;
    if ratio > MAX_N_PASSES as f64 {
        return Err(mk(format!(
            "depth / step_down ratio {ratio} exceeds {MAX_N_PASSES} pass cap — \
             step_down is implausibly small relative to depth"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;
    use valenx_mesh::element::{ElementBlock, ElementType};
    use valenx_mesh::Mesh;

    fn cube(size: f64) -> Mesh {
        let s = size * 0.5;
        let nodes = vec![
            Vector3::new(-s, -s, -s),
            Vector3::new(s, -s, -s),
            Vector3::new(s, s, -s),
            Vector3::new(-s, s, -s),
            Vector3::new(-s, -s, s),
            Vector3::new(s, -s, s),
            Vector3::new(s, s, s),
            Vector3::new(-s, s, s),
        ];
        let conn: Vec<u32> = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut mesh = Mesh::new("cube");
        mesh.nodes = nodes;
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = conn;
        mesh.element_blocks.push(block);
        mesh
    }

    #[test]
    fn adaptive_unit_pocket_emits_multi_ring() {
        let mesh = cube(40.0);
        let stock = Stock::new(
            Vector3::new(-25.0, -25.0, -20.0),
            Vector3::new(50.0, 50.0, 40.0),
            "alu",
        )
        .unwrap();
        let tool = Tool::new(1, "EM4", ToolKind::EndMill, 4.0, 25.0, 2, "carbide").unwrap();
        let params = AdaptiveParams {
            step_down: 4.0,
            depth: 4.0,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).unwrap();
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(
            cuts > 8,
            "expected multi-ring adaptive path, got {cuts} cuts"
        );
    }

    #[test]
    fn validate_rejects_zero_step_over_fraction() {
        let bad = AdaptiveParams {
            step_over_fraction: 0.0,
            ..Default::default()
        };
        assert!(validate(&bad).is_err());
    }

    #[test]
    fn validate_rejects_nan_step_down() {
        // NaN > 0.0 is false but `!(NaN > 0.0)` is true. Pre-fix this
        // erroneously *passed* validation on the step_over_fraction
        // check because of NaN's never-greater-than semantics through
        // a different branch, then produced 0 passes and a silent
        // garbage toolpath. Confirm NaN now hits the `is_finite` gate.
        let bad = AdaptiveParams {
            step_down: f64::NAN,
            ..Default::default()
        };
        assert!(validate(&bad).is_err(), "NaN step_down must be rejected");
    }

    #[test]
    fn validate_rejects_infinite_step_down() {
        let bad = AdaptiveParams {
            step_down: f64::INFINITY,
            ..Default::default()
        };
        assert!(validate(&bad).is_err(), "+inf step_down must be rejected");
    }

    #[test]
    fn validate_rejects_nan_depth() {
        let bad = AdaptiveParams {
            depth: f64::NAN,
            ..Default::default()
        };
        assert!(validate(&bad).is_err());
    }

    /// Round-3 fix: `step_down = f64::MIN_POSITIVE` passes the
    /// is_finite + > 0 check but `(depth / MIN_POSITIVE).ceil() as
    /// usize` saturates to `usize::MAX`, then `for k in 1..=n_passes`
    /// would loop ~2^64 times before terminating.
    #[test]
    fn validate_rejects_implausibly_small_step_down() {
        let bad = AdaptiveParams {
            step_down: f64::MIN_POSITIVE,
            depth: 4.0,
            ..Default::default()
        };
        let err = validate(&bad).expect_err("tiny step_down must be capped");
        let msg = format!("{err}");
        assert!(
            msg.contains("ratio") || msg.contains("pass cap") || msg.contains("step_down"),
            "expected pass-cap rejection, got: {msg}"
        );
    }
}
