//! Profile operation — cut the boundary of the source mesh's
//! cross-section at each Z step-down level.
//!
//! ## Algorithm
//!
//! 1. For each Z level (top of stock, decreasing by `step_down`
//!    until reaching `top_z() - depth`):
//!    1. Cut the source mesh by a horizontal plane at that Z.
//!    2. Stitch the resulting line segments into a closed polygon.
//!    3. Offset the polygon outward by `tool.radius_mm()` so the
//!       cutter centre traces *outside* the part boundary.
//!    4. Reverse the polygon winding if `climb` is false
//!       (conventional cut).
//!    5. Emit: rapid up → rapid to start XY → plunge to Z → cut
//!       around the polygon → rapid up.
//!
//! ## v1 simplifications
//!
//! - **Single closed polygon per Z level.** Multi-region cross-
//!   sections (a Z plane that hits two disjoint parts of the mesh)
//!   produce undefined behaviour — the segment-stitching grabs one
//!   loop and ignores the rest.
//! - **Polygon offset** is per-vertex bisector (see
//!   [`crate::offset`] for limitations).

use nalgebra::Vector3;
use valenx_mesh::cut::{intersect_plane, LineSegment};
use valenx_mesh::Mesh;

use crate::error::CamError;
use crate::offset;
use crate::operation::ProfileParams;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Generate a profile-cut toolpath. See module docs for the algorithm
/// and v1 simplifications.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &ProfileParams,
    tool: &Tool,
) -> Result<Toolpath, CamError> {
    validate(params)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    // Move to safe Z first.
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));

    let n_passes = crate::op::compute_n_passes(params.depth, params.step_down, "profile")?;
    let mut had_polygon = false;
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
        let offsets = offset::polygon(&polygon, tool.radius_mm());
        if offsets.is_empty() {
            continue;
        }
        let mut ring = offsets.into_iter().next().unwrap();
        if !params.climb {
            ring.reverse();
        }
        // Lift the polygon Z to the current pass Z.
        for v in &mut ring {
            v.z = z;
        }
        // Wrap with rapid + plunge + cut + rapid up.
        let start = ring[0];
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
        for v in ring.iter().skip(1) {
            tp.push(Move::new(MoveKind::Cut, *v, params.feed_mm_per_min));
        }
        // Close the loop back to start.
        tp.push(Move::new(MoveKind::Cut, start, params.feed_mm_per_min));
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(start.x, start.y, safe_z),
            0.0,
        ));
        had_polygon = true;
    }
    if !had_polygon {
        return Err(CamError::BadOperation {
            name: "profile".into(),
            reason: "no cross-section polygon could be extracted from the source mesh".into(),
        });
    }
    Ok(tp)
}

fn validate(params: &ProfileParams) -> Result<(), CamError> {
    if !(params.step_down > 0.0) {
        return Err(CamError::BadOperation {
            name: "profile".into(),
            reason: format!("step_down must be > 0 (got {})", params.step_down),
        });
    }
    if !(params.depth > 0.0) {
        return Err(CamError::BadOperation {
            name: "profile".into(),
            reason: format!("depth must be > 0 (got {})", params.depth),
        });
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(CamError::BadOperation {
            name: "profile".into(),
            reason: format!("feed must be > 0 (got {})", params.feed_mm_per_min),
        });
    }
    Ok(())
}

/// Walk a list of line segments and stitch them into a closed
/// polygon by greedily matching endpoints.
///
/// Returns `None` if no closed loop can be assembled. v1 grabs the
/// first loop it finds and ignores other components — multi-region
/// cross-sections are out of scope.
pub fn stitch_segments(segments: &[LineSegment]) -> Option<Vec<Vector3<f64>>> {
    if segments.is_empty() {
        return None;
    }
    let mut remaining: Vec<LineSegment> = segments.to_vec();
    let start = remaining.swap_remove(0);
    let mut polygon = vec![start.a, start.b];
    let mut current = start.b;
    let start_pt = start.a;
    const EPS_SQ: f64 = 1e-8;
    loop {
        // Find a segment whose endpoint matches current.
        let mut found: Option<usize> = None;
        let mut reverse = false;
        for (i, s) in remaining.iter().enumerate() {
            if (s.a - current).norm_squared() < EPS_SQ {
                found = Some(i);
                reverse = false;
                break;
            }
            if (s.b - current).norm_squared() < EPS_SQ {
                found = Some(i);
                reverse = true;
                break;
            }
        }
        let idx = match found {
            Some(i) => i,
            None => break, // Open chain — bail.
        };
        let next = remaining.swap_remove(idx);
        let next_pt = if reverse { next.a } else { next.b };
        // Check if we've closed the loop.
        if (next_pt - start_pt).norm_squared() < EPS_SQ {
            // Don't push duplicate of start.
            return Some(polygon);
        }
        polygon.push(next_pt);
        current = next_pt;
        if polygon.len() > segments.len() + 2 {
            // Sanity bail.
            break;
        }
    }
    if polygon.len() >= 3 {
        Some(polygon)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operation::ProfileParams;
    use crate::stock::Stock;
    use crate::tool::{Tool, ToolKind};
    use valenx_mesh::element::{ElementBlock, ElementType};
    use valenx_mesh::Mesh;

    /// Build a unit-cube triangle mesh centred at origin, +Z up.
    fn cube_mesh(size: f64) -> Mesh {
        let s = size * 0.5;
        let nodes = vec![
            // bottom face (z = -s)
            Vector3::new(-s, -s, -s),
            Vector3::new(s, -s, -s),
            Vector3::new(s, s, -s),
            Vector3::new(-s, s, -s),
            // top face (z = +s)
            Vector3::new(-s, -s, s),
            Vector3::new(s, -s, s),
            Vector3::new(s, s, s),
            Vector3::new(-s, s, s),
        ];
        // 12 triangles, 2 per face.
        let conn: Vec<u32> = vec![
            // bottom (CW from below)
            0, 2, 1, 0, 3, 2, // top
            4, 5, 6, 4, 6, 7, // front (y=-s)
            0, 1, 5, 0, 5, 4, // back (y=+s)
            2, 3, 7, 2, 7, 6, // left (x=-s)
            0, 4, 7, 0, 7, 3, // right (x=+s)
            1, 2, 6, 1, 6, 5,
        ];
        let mut mesh = Mesh::new("cube");
        mesh.nodes = nodes;
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = conn;
        mesh.element_blocks.push(block);
        mesh
    }

    #[test]
    fn profile_unit_cube_produces_toolpath() {
        let mesh = cube_mesh(10.0); // 10×10×10 cube centred at origin
                                    // Stock placed so the cube sits on top (cube top z = +5, stock
                                    // top z = +5).
        let stock = Stock::new(
            Vector3::new(-6.0, -6.0, -5.0),
            Vector3::new(12.0, 12.0, 10.0),
            "wood",
        )
        .unwrap();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let params = ProfileParams {
            tool_id: 1,
            step_down: 2.0,
            depth: 6.0,
            ..Default::default()
        };
        let tp = generate(&stock, &mesh, &params, &tool).expect("profile should generate");
        assert!(!tp.is_empty(), "toolpath should have moves");
        // 3 passes (z=3, 1, -1) ⇒ at least 3 plunges.
        let plunges = tp
            .moves
            .iter()
            .filter(|m| m.kind == MoveKind::Plunge)
            .count();
        assert!(plunges >= 3, "expected ≥3 plunges, got {plunges}");
        let cuts = tp.moves.iter().filter(|m| m.kind == MoveKind::Cut).count();
        assert!(cuts >= 12, "expected ≥12 cuts, got {cuts}");
    }

    #[test]
    fn stitch_simple_square() {
        let segments = vec![
            LineSegment {
                a: Vector3::new(0.0, 0.0, 0.0),
                b: Vector3::new(1.0, 0.0, 0.0),
            },
            LineSegment {
                a: Vector3::new(1.0, 0.0, 0.0),
                b: Vector3::new(1.0, 1.0, 0.0),
            },
            LineSegment {
                a: Vector3::new(1.0, 1.0, 0.0),
                b: Vector3::new(0.0, 1.0, 0.0),
            },
            LineSegment {
                a: Vector3::new(0.0, 1.0, 0.0),
                b: Vector3::new(0.0, 0.0, 0.0),
            },
        ];
        let poly = stitch_segments(&segments).expect("should close");
        assert_eq!(poly.len(), 4);
    }

    #[test]
    fn validate_rejects_zero_step_down() {
        let bad = ProfileParams {
            step_down: 0.0,
            ..Default::default()
        };
        let err = validate(&bad).unwrap_err();
        assert_eq!(err.code(), "cam.bad_operation");
    }

    #[test]
    fn climb_reverses_polygon() {
        let mesh = cube_mesh(10.0);
        let stock = Stock::new(
            Vector3::new(-6.0, -6.0, -5.0),
            Vector3::new(12.0, 12.0, 10.0),
            "wood",
        )
        .unwrap();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let p_climb = ProfileParams {
            step_down: 2.0,
            depth: 2.0,
            climb: true,
            ..Default::default()
        };
        let p_conv = ProfileParams {
            step_down: 2.0,
            depth: 2.0,
            climb: false,
            ..Default::default()
        };
        let tp_climb = generate(&stock, &mesh, &p_climb, &tool).unwrap();
        let tp_conv = generate(&stock, &mesh, &p_conv, &tool).unwrap();
        // Find the first cut moves in each — they should differ (reversed direction).
        let first_cut_climb = tp_climb
            .moves
            .iter()
            .find(|m| m.kind == MoveKind::Cut)
            .unwrap();
        let first_cut_conv = tp_conv
            .moves
            .iter()
            .find(|m| m.kind == MoveKind::Cut)
            .unwrap();
        // Cut positions differ because the polygon was reversed.
        assert_ne!(first_cut_climb.position, first_cut_conv.position);
    }
}
