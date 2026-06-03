//! Toolpath simulation helpers — group moves into polylines for
//! rendering, estimate cycle time, estimate material removed.
//!
//! The actual viewport painting lives in the host app
//! (`crates/valenx-app/src/cam_overlay.rs`); this module produces
//! the data structures the renderer consumes.

use nalgebra::Vector3;

use crate::{
    stock::Stock,
    tool::Tool,
    toolpath::{MoveKind, Toolpath},
    voxel::{Voxel, VoxelError},
};

/// Rapid traversal feed assumed for cycle-time estimation when no
/// explicit feed is set (mm/min). Typical hobby spindle rapids at
/// 5000 mm/min; production machines go much faster.
pub const RAPID_FEED_MM_PER_MIN: f64 = 5000.0;

/// Group consecutive moves of the same kind into polylines.
///
/// Returns a list of `(kind, points)`. Each polyline starts from
/// the *end* of the previous move (or the very first position for
/// the first polyline), so the painter can draw `points.windows(2)`
/// as line segments.
pub fn toolpath_polylines(toolpath: &Toolpath) -> Vec<(MoveKind, Vec<Vector3<f64>>)> {
    let mut out: Vec<(MoveKind, Vec<Vector3<f64>>)> = Vec::new();
    let mut prev_pos: Option<Vector3<f64>> = None;
    let mut current_kind: Option<MoveKind> = None;
    let mut current_pts: Vec<Vector3<f64>> = Vec::new();
    for m in &toolpath.moves {
        match current_kind {
            Some(k) if k == m.kind => {
                current_pts.push(m.position);
            }
            _ => {
                // Flush previous polyline.
                if let Some(k) = current_kind {
                    out.push((k, std::mem::take(&mut current_pts)));
                }
                // Start a new polyline. Seed with the previous
                // end-of-segment so the new polyline starts where
                // the previous one finished.
                if let Some(p) = prev_pos {
                    current_pts.push(p);
                }
                current_pts.push(m.position);
                current_kind = Some(m.kind);
            }
        }
        prev_pos = Some(m.position);
    }
    if let Some(k) = current_kind {
        if !current_pts.is_empty() {
            out.push((k, current_pts));
        }
    }
    out
}

/// Estimated cycle time in minutes.
///
/// For each move segment, time = distance / feed. Rapid moves use
/// [`RAPID_FEED_MM_PER_MIN`]; other moves use their per-move feed.
/// Zero-feed cut/plunge moves are skipped (would be `inf` time).
pub fn estimated_time(toolpath: &Toolpath) -> f64 {
    let mut total = 0.0;
    for w in toolpath.moves.windows(2) {
        let dist = (w[1].position - w[0].position).norm();
        let feed = match w[1].kind {
            MoveKind::Rapid => RAPID_FEED_MM_PER_MIN,
            _ => {
                if w[1].feed > 0.0 {
                    w[1].feed
                } else {
                    continue;
                }
            }
        };
        total += dist / feed;
    }
    total
}

/// Default voxel resolution for material-removal simulation. 64³
/// keeps memory at ~32 KiB per stock — interactive enough for
/// preview while still useful for visualisation.
pub const DEFAULT_VOXEL_RES: u32 = 64;

/// Initialise a voxel grid from a stock block at the given
/// resolution.
///
/// # Errors
///
/// Forwards [`VoxelError`] from [`Voxel::from_aabb`] when the cubed
/// resolution overflows or exceeds the cell-count cap.
pub fn voxel_from_stock(stock: &Stock, resolution: u32) -> Result<Voxel, VoxelError> {
    let (min, max) = stock.aabb();
    Voxel::from_aabb(min, max, (resolution, resolution, resolution))
}

/// Produce `n_frames` intermediate voxel states for a swept toolpath.
///
/// Frame 0 = full stock, frame N-1 = final state. Cut and plunge
/// segments remove material; rapid segments are ignored (the tool is
/// not in contact).
///
/// Returns a `Vec` of length `n_frames` containing extracted boundary
/// meshes (one per frame).
///
/// # Errors
///
/// Forwards [`VoxelError`] from voxel grid construction when the
/// requested resolution is too large for the per-axis cap.
pub fn animate(
    toolpath: &Toolpath,
    stock: &Stock,
    tool: &Tool,
    n_frames: u32,
    resolution: u32,
) -> Result<Vec<valenx_mesh::Mesh>, VoxelError> {
    let mut voxel = voxel_from_stock(stock, resolution)?;
    let mut frames = Vec::with_capacity(n_frames as usize);
    // Frame 0 = full stock.
    frames.push(voxel.to_mesh());
    if n_frames <= 1 || toolpath.moves.is_empty() {
        return Ok(frames);
    }
    let n_moves = toolpath.moves.len();
    let n_intervals = n_frames - 1;
    let moves_per_interval = ((n_moves as f64) / (n_intervals as f64)).ceil() as usize;
    let r = tool.radius_mm();
    let mut i = 1_usize;
    let mut placed = 1_u32;
    while i < n_moves && placed < n_frames {
        let end = (i + moves_per_interval).min(n_moves);
        // Sweep this slice of moves into the voxel grid.
        for k in i..end {
            let prev = toolpath.moves[k - 1].position;
            let cur = toolpath.moves[k].position;
            match toolpath.moves[k].kind {
                MoveKind::Cut | MoveKind::Plunge => {
                    voxel.cut_segment(prev, cur, r);
                }
                MoveKind::Arc { centre_xy, dir } => {
                    // Approximate the swept arc volume by carving along
                    // the polyline expansion of the arc.
                    for (a, b) in arc_polyline_segments(prev, cur, centre_xy, dir) {
                        voxel.cut_segment(a, b, r);
                    }
                }
                MoveKind::Rapid => {}
            }
        }
        frames.push(voxel.to_mesh());
        i = end;
        placed += 1;
    }
    // Pad with the final state if we ran short.
    while frames.len() < n_frames as usize {
        frames.push(voxel.to_mesh());
    }
    Ok(frames)
}

/// Single-shot final-state voxel mesh — faster than [`animate`] when
/// the host only needs the end result.
///
/// # Errors
///
/// Forwards [`VoxelError`] from voxel grid construction when the
/// requested resolution is too large for the per-axis cap.
pub fn final_state(
    toolpath: &Toolpath,
    stock: &Stock,
    tool: &Tool,
    resolution: u32,
) -> Result<valenx_mesh::Mesh, VoxelError> {
    let mut voxel = voxel_from_stock(stock, resolution)?;
    let r = tool.radius_mm();
    for w in toolpath.moves.windows(2) {
        match w[1].kind {
            MoveKind::Cut | MoveKind::Plunge => {
                voxel.cut_segment(w[0].position, w[1].position, r);
            }
            MoveKind::Arc { centre_xy, dir } => {
                for (a, b) in arc_polyline_segments(w[0].position, w[1].position, centre_xy, dir) {
                    voxel.cut_segment(a, b, r);
                }
            }
            MoveKind::Rapid => {}
        }
    }
    Ok(voxel.to_mesh())
}

/// Tessellate an XY-plane arc move into a polyline of `~5°` segments
/// for material-removal + length estimates. Z is linearly
/// interpolated from start to end.
pub fn arc_polyline_segments(
    start: Vector3<f64>,
    end: Vector3<f64>,
    centre_xy: nalgebra::Vector2<f64>,
    dir: crate::arcfit::ArcDir,
) -> Vec<(Vector3<f64>, Vector3<f64>)> {
    let cx = centre_xy.x;
    let cy = centre_xy.y;
    let r0_x = start.x - cx;
    let r0_y = start.y - cy;
    let r1_x = end.x - cx;
    let r1_y = end.y - cy;
    let r = (r0_x * r0_x + r0_y * r0_y).sqrt().max(1e-9);
    let theta0 = r0_y.atan2(r0_x);
    let theta1 = r1_y.atan2(r1_x);
    // Compute signed arc span.
    let mut delta = theta1 - theta0;
    let two_pi = std::f64::consts::TAU;
    match dir {
        crate::arcfit::ArcDir::Counterclockwise => {
            while delta < 0.0 {
                delta += two_pi;
            }
        }
        crate::arcfit::ArcDir::Clockwise => {
            while delta > 0.0 {
                delta -= two_pi;
            }
        }
    }
    let n_seg = (delta.abs() / (5.0_f64).to_radians()).ceil().max(1.0) as usize;
    let mut segs = Vec::with_capacity(n_seg);
    let mut prev = start;
    for k in 1..=n_seg {
        let t = (k as f64) / (n_seg as f64);
        let theta = theta0 + delta * t;
        let x = cx + r * theta.cos();
        let y = cy + r * theta.sin();
        let z = start.z + (end.z - start.z) * t;
        let p = Vector3::new(x, y, z);
        segs.push((prev, p));
        prev = p;
    }
    segs
}

/// Estimated material removed in mm³.
///
/// v1 approximation: each `Cut` move's removed volume ≈
/// `tool.diameter × cut_distance × cut_depth_per_pass`. We don't
/// know the per-pass depth here, so we use the full tool length as
/// an upper bound — callers should treat this as an *upper bound*
/// estimate only.
///
/// Plunge moves add `π × (tool.diameter / 2)² × distance`
/// (cylindrical hole volume). Rapids contribute nothing.
pub fn removed_volume_mm3(toolpath: &Toolpath, tool: &Tool) -> f64 {
    let mut total = 0.0;
    let r = tool.radius_mm();
    let area_cylinder = std::f64::consts::PI * r * r;
    for w in toolpath.moves.windows(2) {
        let dist = (w[1].position - w[0].position).norm();
        match w[1].kind {
            MoveKind::Cut => {
                // Approximate as a swept rectangle of width = tool diameter
                // and an unknown axial-engagement depth. Using r as a
                // conservative engagement depth so the estimate scales
                // with tool size.
                total += dist * tool.diameter_mm * r;
            }
            MoveKind::Arc { centre_xy, dir } => {
                // Use the arc length (tessellated) instead of the
                // straight-line chord distance.
                let segs =
                    arc_polyline_segments(w[0].position, w[1].position, centre_xy, dir);
                let arc_len: f64 = segs.iter().map(|(a, b)| (b - a).norm()).sum();
                total += arc_len * tool.diameter_mm * r;
            }
            MoveKind::Plunge => {
                total += dist * area_cylinder;
            }
            MoveKind::Rapid => {}
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::{Tool, ToolKind};
    use crate::toolpath::Move;

    fn p(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    #[test]
    fn polylines_group_by_kind() {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Rapid, p(10.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Plunge, p(10.0, 0.0, 0.0), 200.0));
        tp.push(Move::new(MoveKind::Cut, p(20.0, 0.0, 0.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(20.0, 10.0, 0.0), 500.0));
        let polys = toolpath_polylines(&tp);
        // 3 polylines: rapid (2 pts), plunge (2 pts), cut (3 pts incl. start)
        assert_eq!(polys.len(), 3);
        assert_eq!(polys[0].0, MoveKind::Rapid);
        assert_eq!(polys[0].1.len(), 2);
        assert_eq!(polys[1].0, MoveKind::Plunge);
        assert_eq!(polys[2].0, MoveKind::Cut);
        assert_eq!(polys[2].1.len(), 3);
    }

    #[test]
    fn estimated_time_positive() {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Plunge, p(0.0, 0.0, 0.0), 200.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 0.0), 500.0));
        let t = estimated_time(&tp);
        // 5mm at 200mm/min = 0.025min, 10mm at 500mm/min = 0.02min => 0.045min
        assert!((t - 0.045).abs() < 1e-9, "got {t}");
    }

    #[test]
    fn estimated_time_skips_zero_feed_cut() {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 5.0), 0.0));
        let t = estimated_time(&tp);
        assert!(t.is_finite(), "should not be infinite");
        // Time = 0 since the only candidate move has feed 0.
        assert!(t < 1e-12);
    }

    #[test]
    fn removed_volume_positive_for_cut_path() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Plunge, p(0.0, 0.0, 0.0), 200.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 0.0), 500.0));
        let v = removed_volume_mm3(&tp, &tool);
        assert!(v > 0.0, "expected positive removed volume, got {v}");
    }

    #[test]
    fn empty_toolpath_yields_zeros() {
        let tp = Toolpath::new();
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        assert!(toolpath_polylines(&tp).is_empty());
        assert!(estimated_time(&tp) < 1e-12);
        assert!(removed_volume_mm3(&tp, &tool) < 1e-12);
    }

    #[test]
    fn animate_produces_n_frames() {
        let stock =
            crate::stock::Stock::new(Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), "alu")
                .unwrap();
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 5.0, 11.0), 0.0));
        tp.push(Move::new(MoveKind::Plunge, p(0.0, 5.0, 5.0), 200.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 5.0, 5.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 5.0, 0.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(0.0, 5.0, 0.0), 500.0));
        let frames = animate(&tp, &stock, &tool, 5, 16).unwrap();
        assert_eq!(frames.len(), 5);
        // Frame 0 = full stock; final frame has fewer tris.
        let n0 = frames[0].element_blocks[0].connectivity.len();
        let nn = frames[4].element_blocks[0].connectivity.len();
        assert!(
            nn != n0,
            "expected animation frames to differ from start ({n0} vs {nn})"
        );
    }

    #[test]
    fn final_state_reduces_voxel_count() {
        let stock =
            crate::stock::Stock::new(Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), "alu")
                .unwrap();
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "").unwrap();
        let v0 = voxel_from_stock(&stock, 16).unwrap();
        let initial = v0.solid_count();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 5.0, 11.0), 0.0));
        tp.push(Move::new(MoveKind::Plunge, p(0.0, 5.0, 5.0), 200.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 5.0, 5.0), 500.0));
        let mesh = final_state(&tp, &stock, &tool, 16).unwrap();
        // Mesh should have triangles (boundary of remaining material).
        assert!(!mesh.element_blocks.is_empty());
        // Final voxel state has fewer solids than initial.
        let mut v_final = voxel_from_stock(&stock, 16).unwrap();
        for w in tp.moves.windows(2) {
            if w[1].kind != MoveKind::Rapid {
                v_final.cut_segment(w[0].position, w[1].position, tool.radius_mm());
            }
        }
        assert!(v_final.solid_count() < initial);
    }

    #[test]
    fn square_pocket_toolpath_time_estimate_reasonable() {
        // Synthesise a small pocket: plunge + 4 cuts around a square.
        let tool = Tool::new(1, "EM2", ToolKind::EndMill, 2.0, 20.0, 2, "carbide").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Plunge, p(0.0, 0.0, -1.0), 200.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, -1.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 10.0, -1.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(0.0, 10.0, -1.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(0.0, 0.0, -1.0), 500.0));
        let t = estimated_time(&tp);
        assert!(
            t > 0.0 && t < 1.0,
            "1-min upper bound for this pocket; got {t}"
        );
        let v = removed_volume_mm3(&tp, &tool);
        assert!(v > 0.0);
    }
}
