//! Face operation — level the top of the stock with parallel raster
//! passes.
//!
//! ## Algorithm
//!
//! For each Z step-down (from `top_z()` down to `top_z() - depth`):
//!
//! 1. Build a rectangular polygon over the stock's XY extent
//!    (slightly inset by `tool.radius_mm()` so the cutter stays
//!    over material).
//! 2. Generate parallel raster lines at `raster_angle_deg` spaced
//!    `step_over` apart.
//! 3. Reverse every other line if `climb` is true (zig-zag fill),
//!    or keep all in the same direction for conventional facing.
//! 4. Emit rapid → plunge → cut sequences (one plunge per pass).

use nalgebra::Vector3;

use crate::error::CamError;
use crate::operation::FaceParams;
use crate::raster;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Generate a face-mill toolpath spanning the stock's XY extent.
pub fn generate(stock: &Stock, params: &FaceParams, tool: &Tool) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));

    let (min, max) = stock.aabb();
    // Inset by tool radius so the cutter stays over material on the
    // first / last raster line.
    let r = tool.radius_mm();
    let polygon = vec![
        Vector3::new(min.x + r, min.y + r, 0.0),
        Vector3::new(max.x - r, min.y + r, 0.0),
        Vector3::new(max.x - r, max.y - r, 0.0),
        Vector3::new(min.x + r, max.y - r, 0.0),
    ];

    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "face")?;
    for k in 1..=n_passes {
        let depth_below_top = (params.step_down * k as f64).min(params.depth);
        let z = stock.top_z() - depth_below_top;
        // Climb = zig-zag (alternating); conventional = parallel
        // (one-way).
        let lines = if params.climb {
            raster::zigzag(&polygon, params.step_over, params.raster_angle_deg)
        } else {
            raster::parallel(&polygon, params.step_over, params.raster_angle_deg)
        };
        emit_raster(&mut tp, &lines, z, safe_z, params);
    }
    Ok(tp)
}

fn emit_raster(
    tp: &mut Toolpath,
    lines: &[Vec<Vector3<f64>>],
    z: f64,
    safe_z: f64,
    params: &FaceParams,
) {
    if lines.is_empty() {
        return;
    }
    // Climb / zigzag → cut directly between line endpoints (no
    // retract). Conventional / parallel → retract between passes.
    if params.climb {
        let first = lines[0][0];
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
        for v in lines[0].iter().skip(1) {
            tp.push(Move::new(
                MoveKind::Cut,
                Vector3::new(v.x, v.y, z),
                params.feed_mm_per_min,
            ));
        }
        for line in lines.iter().skip(1) {
            for v in line {
                tp.push(Move::new(
                    MoveKind::Cut,
                    Vector3::new(v.x, v.y, z),
                    params.feed_mm_per_min,
                ));
            }
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(0.0, 0.0, safe_z),
            0.0,
        ));
    } else {
        for line in lines {
            if line.len() < 2 {
                continue;
            }
            let start = line[0];
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
            for v in line.iter().skip(1) {
                tp.push(Move::new(
                    MoveKind::Cut,
                    Vector3::new(v.x, v.y, z),
                    params.feed_mm_per_min,
                ));
            }
            let end = *line.last().unwrap();
            tp.push(Move::new(
                MoveKind::Rapid,
                Vector3::new(end.x, end.y, safe_z),
                0.0,
            ));
        }
    }
}

fn validate(params: &FaceParams) -> Result<(), CamError> {
    if !(params.step_down > 0.0) {
        return Err(CamError::BadOperation {
            name: "face".into(),
            reason: format!("step_down must be > 0 (got {})", params.step_down),
        });
    }
    if !(params.step_over > 0.0) {
        return Err(CamError::BadOperation {
            name: "face".into(),
            reason: format!("step_over must be > 0 (got {})", params.step_over),
        });
    }
    if !(params.depth > 0.0) {
        return Err(CamError::BadOperation {
            name: "face".into(),
            reason: format!("depth must be > 0 (got {})", params.depth),
        });
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(CamError::BadOperation {
            name: "face".into(),
            reason: format!("feed must be > 0 (got {})", params.feed_mm_per_min),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::{Tool, ToolKind};

    #[test]
    fn face_rectangular_stock_parallel_raster() {
        let stock = Stock::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(50.0, 30.0, 10.0),
            "wood",
        )
        .unwrap();
        let tool = Tool::new(1, "FM10", ToolKind::FaceMill, 10.0, 30.0, 4, "carbide").unwrap();
        let params = FaceParams {
            step_over: 4.0,
            step_down: 0.5,
            depth: 0.5,
            climb: false, // parallel — retract between passes
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        // Stock Y extent = 30, tool radius inset = 5 each side, so
        // usable = 20. step_over = 4 → ~5 raster lines → 5 plunges.
        assert!(
            (4..=7).contains(&plunges),
            "expected 4–7 plunges, got {plunges}"
        );
    }

    #[test]
    fn face_climb_uses_zigzag() {
        let stock = Stock::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(50.0, 30.0, 10.0),
            "wood",
        )
        .unwrap();
        let tool = Tool::new(1, "FM10", ToolKind::FaceMill, 10.0, 30.0, 4, "carbide").unwrap();
        let params = FaceParams {
            step_over: 4.0,
            step_down: 0.5,
            depth: 0.5,
            climb: true,
            ..Default::default()
        };
        let tp = generate(&stock, &params, &tool).unwrap();
        // Climb / zigzag → only 1 plunge per pass.
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert_eq!(
            plunges, 1,
            "zigzag should plunge once per pass; got {plunges}"
        );
    }

    #[test]
    fn face_validate_rejects_zero_depth() {
        let bad = FaceParams {
            depth: 0.0,
            ..Default::default()
        };
        assert!(validate(&bad).is_err());
    }
}
