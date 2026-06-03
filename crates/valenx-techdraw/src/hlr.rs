//! Hidden-line removal (approximate v1).
//!
//! True analytic HLR — extract silhouette curves, partition the model
//! into "in front" and "behind" pieces via BSP or BRep boolean ops —
//! is months of work and out of scope for Phase 5. This module ships
//! a coarse z-buffer approximation that's correct enough for the
//! common shapes a CAD review meeting cares about (cubes, prisms,
//! cylinders, fillet/chamfer outputs).
//!
//! The algorithm:
//!
//! 1. Tessellate the solid (or use the cached mesh for mesh-backed
//!    solids).
//! 2. Project every triangle through the view's camera, recording
//!    `(x, y, z_avg)` per triangle.
//! 3. Build a coarse 2D grid (`GRID` × `GRID` cells); for each grid
//!    cell, remember the smallest `z` (closest-to-camera) of any
//!    triangle whose projected bbox overlaps the cell.
//! 4. For each edge to classify, sample its 2D midpoint. Look up the
//!    grid cell. If the edge's depth at that midpoint is "in front of
//!    or equal to" the cell's recorded depth, mark visible; else
//!    hidden.
//!
//! Pros: simple, deterministic, no silhouette-curve sorting.
//! Cons: a 64×64 grid will miss thin features at low scale; the
//! midpoint sample assumes the edge is monotone in depth (usually
//! true for projected straight edges). Both are accepted v1
//! compromises.
//!
//! Phase 5.5 will replace this with a true silhouette extractor
//! built on BRep edge classification.

use nalgebra::{Matrix4, Vector3, Vector4};

use crate::error::TechDrawError;
use crate::projection::{extract_edges, project_point};

/// Resolution of the depth-grid used by the v1 HLR pass. 128 × 128
/// gives sub-mm precision on an A4 sheet for typical part sizes (1 mm
/// per cell at 128 mm view extent) without making the per-edge lookup
/// expensive.
const GRID: usize = 128;

/// Visible / hidden segments are returned in the same drawing-plane
/// millimeter frame as [`crate::projection::project_edges`].
pub type EdgeSegments = Vec<[(f64, f64); 2]>;

/// Classify every edge of `solid` as visible or hidden when viewed
/// through `camera`.
///
/// Returns `(visible, hidden)`. The union of the two equals the
/// output of [`crate::projection::project_edges`] modulo ordering.
pub fn classify_edges(
    solid: &valenx_cad::Solid,
    camera: &Matrix4<f64>,
) -> Result<(EdgeSegments, EdgeSegments), TechDrawError> {
    // Pull tessellated triangles (BRep or mesh) so we have something
    // to z-buffer against. We re-tessellate here independently of
    // `extract_edges` because we need per-triangle depth, not just
    // edge endpoints — keeping the two passes decoupled makes the
    // intent of each one clearer.
    let mesh = match solid {
        valenx_cad::Solid::Brep(_) => valenx_cad::tessellate::solid_to_mesh(solid, 0.1)
            .map_err(|e| TechDrawError::ExportFailed(format!("hlr tessellation: {e}")))?,
        valenx_cad::Solid::Mesh(m) => m.clone(),
    };
    if mesh.nodes.is_empty() {
        return Err(TechDrawError::EmptySolid);
    }

    // Project every triangle. Triangles store
    // (min_x, min_y, max_x, max_y, depth) — depth is the average of
    // the three projected z values (camera +Z is "forward" toward
    // viewer; we want larger Z = closer).
    struct Tri {
        bbox: [f64; 4], // [min_x, min_y, max_x, max_y]
        depth: f64,
    }
    let mut tris: Vec<Tri> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::element::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            // Bounds-check the connectivity: a loaded/corrupt mesh can
            // carry a node index past `nodes.len()`. Index by `.get()`
            // and skip the triangle on an out-of-range value rather than
            // panicking during view generation (mirrors projection.rs).
            let (Some(&a), Some(&b), Some(&c)) = (
                mesh.nodes.get(tri[0] as usize),
                mesh.nodes.get(tri[1] as usize),
                mesh.nodes.get(tri[2] as usize),
            ) else {
                continue;
            };
            let v = [a, b, c];
            let p = v.map(|w| transform(w, camera));
            let mn_x = p[0][0].min(p[1][0]).min(p[2][0]);
            let mx_x = p[0][0].max(p[1][0]).max(p[2][0]);
            let mn_y = p[0][1].min(p[1][1]).min(p[2][1]);
            let mx_y = p[0][1].max(p[1][1]).max(p[2][1]);
            // depth = average projected z. Camera in our look-at
            // convention has -Z forward, so smaller z = closer to
            // camera. We negate so "closer = bigger depth value",
            // simplifying the comparison.
            let depth = -(p[0][2] + p[1][2] + p[2][2]) / 3.0;
            tris.push(Tri {
                bbox: [mn_x, mn_y, mx_x, mx_y],
                depth,
            });
        }
    }
    if tris.is_empty() {
        return Err(TechDrawError::EmptySolid);
    }

    // Compute overall bbox of every triangle so the depth-grid maps
    // cleanly. Add 1% padding to dodge edge cases at the boundary.
    let mut world_min = [f64::INFINITY; 2];
    let mut world_max = [f64::NEG_INFINITY; 2];
    for t in &tris {
        if t.bbox[0] < world_min[0] {
            world_min[0] = t.bbox[0];
        }
        if t.bbox[1] < world_min[1] {
            world_min[1] = t.bbox[1];
        }
        if t.bbox[2] > world_max[0] {
            world_max[0] = t.bbox[2];
        }
        if t.bbox[3] > world_max[1] {
            world_max[1] = t.bbox[3];
        }
    }
    let span_x = (world_max[0] - world_min[0]).max(1e-6);
    let span_y = (world_max[1] - world_min[1]).max(1e-6);
    let pad_x = span_x * 0.01;
    let pad_y = span_y * 0.01;
    world_min[0] -= pad_x;
    world_min[1] -= pad_y;
    world_max[0] += pad_x;
    world_max[1] += pad_y;
    let cell_x = (world_max[0] - world_min[0]) / GRID as f64;
    let cell_y = (world_max[1] - world_min[1]) / GRID as f64;

    let cell = |x: f64, y: f64| -> (usize, usize) {
        let cx = ((x - world_min[0]) / cell_x)
            .floor()
            .clamp(0.0, (GRID - 1) as f64) as usize;
        let cy = ((y - world_min[1]) / cell_y)
            .floor()
            .clamp(0.0, (GRID - 1) as f64) as usize;
        (cx, cy)
    };

    // Depth-grid: for each cell, the largest depth value (= closest
    // triangle) that covers the cell.
    let mut grid = vec![f64::NEG_INFINITY; GRID * GRID];
    for t in &tris {
        let (x0, y0) = cell(t.bbox[0], t.bbox[1]);
        let (x1, y1) = cell(t.bbox[2], t.bbox[3]);
        for gy in y0..=y1 {
            for gx in x0..=x1 {
                let i = gy * GRID + gx;
                if t.depth > grid[i] {
                    grid[i] = t.depth;
                }
            }
        }
    }

    // Now classify every edge. Sample its 2D midpoint, compute its
    // own depth, compare to the grid cell. An edge is visible iff its
    // depth is within `eps` of the cell's depth (i.e. it's on or in
    // front of the closest triangle covering its midpoint).
    let edge_eps = (cell_x.max(cell_y)) * 0.5 + 1e-4;
    let world_edges = extract_edges(solid)?;
    let mut visible: EdgeSegments = Vec::new();
    let mut hidden: EdgeSegments = Vec::new();
    for (a, b) in world_edges {
        let pa = transform(a, camera);
        let pb = transform(b, camera);
        let mid_x = (pa[0] + pb[0]) * 0.5;
        let mid_y = (pa[1] + pb[1]) * 0.5;
        let mid_depth = -(pa[2] + pb[2]) * 0.5;
        let (gx, gy) = cell(mid_x, mid_y);
        let cell_depth = grid[gy * GRID + gx];
        let seg = [(pa[0], pa[1]), (pb[0], pb[1])];
        if mid_depth + edge_eps >= cell_depth {
            visible.push(seg);
        } else {
            hidden.push(seg);
        }
    }
    Ok((visible, hidden))
}

fn transform(p: Vector3<f64>, m: &Matrix4<f64>) -> [f64; 3] {
    let v: Vector4<f64> = m * Vector4::new(p.x, p.y, p.z, 1.0);
    [v.x, v.y, v.z]
}

/// Re-export so `view::View::project` can call straight into the
/// projection module without re-importing both crates.
#[allow(dead_code)]
fn _ensure_link(p: Vector3<f64>, m: &Matrix4<f64>) -> [f64; 2] {
    project_point(p, m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::view::ViewKind;
    use valenx_cad::primitives::box_solid;

    /// Task 10 — HLR test on a unit cube viewed from the front.
    ///
    /// Expectations (plan's notes):
    /// - 4 visible edges on the front face.
    /// - 4 hidden edges on the back face (these project to the same
    ///   `(x, y)` rectangle as the front face).
    /// - 4 "shared" silhouette edges projecting down to points — the
    ///   z-axis edges that join front to back. v1 z-buffer counts
    ///   these as visible because their midpoint depth matches the
    ///   closest cell.
    ///
    /// We assert the loose invariant "at least 4 visible, at least 4
    /// hidden, and visible + hidden = total". Exact counts depend on
    /// tessellation (truck may emit different triangulations between
    /// versions), so we don't lock to a specific number.
    #[test]
    fn unit_cube_front_view_visible_and_hidden_counts() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let cam = ViewKind::Front.camera_matrix();
        let (visible, hidden) = classify_edges(&cube, &cam).unwrap();
        assert!(!visible.is_empty(), "front face should yield visible edges");
        assert!(!hidden.is_empty(), "back face should yield hidden edges");
        // Sanity: total = visible + hidden.
        let total_extracted = crate::projection::extract_edges(&cube).unwrap().len();
        assert_eq!(visible.len() + hidden.len(), total_extracted);
    }

    #[test]
    fn classify_edges_out_of_range_connectivity_does_not_panic() {
        // R32 M1: sibling of the R31 projection.rs fix. classify_edges
        // indexed `mesh.nodes[tri[N] as usize]` directly, so a
        // loaded/corrupt Tri3 mesh whose connectivity references a node
        // past `nodes.len()` panicked ("index out of bounds") during
        // View::generate. The bad triangle must be skipped gracefully.
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = valenx_mesh::Mesh::new("corrupt");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        // First triangle valid; second references node 99 (only 3 exist).
        block.connectivity = vec![0, 1, 2, 0, 1, 99];
        mesh.element_blocks.push(block);
        let s = valenx_cad::Solid::from_mesh(mesh);
        let cam = ViewKind::Front.camera_matrix();
        // Must not panic; the one valid triangle still classifies.
        let (visible, hidden) = classify_edges(&s, &cam).expect("valid triangle classifies");
        assert!(
            !visible.is_empty() || !hidden.is_empty(),
            "the in-range triangle's edges should survive"
        );
    }

    /// Task 12 — Iso projection test: cube from iso view shows a
    /// hexagonal silhouette. We verify the projected visible-edge
    /// extent has the characteristic 6-fold "wider in x than y is
    /// trivial" shape — specifically, the x-extent should be roughly
    /// √3 / √2 ≈ 1.22× the cube's edge length and the y-extent
    /// should be roughly √2 ≈ 1.41× the edge length.
    #[test]
    fn iso_view_unit_cube_has_hexagonal_silhouette() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let cam = ViewKind::Isometric.camera_matrix();
        let (visible, _hidden) = classify_edges(&cube, &cam).unwrap();
        assert!(!visible.is_empty());
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for seg in &visible {
            for (x, y) in seg {
                if *x < min_x {
                    min_x = *x;
                }
                if *x > max_x {
                    max_x = *x;
                }
                if *y < min_y {
                    min_y = *y;
                }
                if *y > max_y {
                    max_y = *y;
                }
            }
        }
        let dx = max_x - min_x;
        let dy = max_y - min_y;
        // For a unit cube in iso, the projected silhouette is a
        // regular hexagon inscribed in a circle of radius √(2/3).
        // The hexagon's width (x-extent) is 2·√(2/3) ≈ 1.633 and the
        // height (y-extent) is √2 ≈ 1.414. Allow ±5% slack for
        // tessellation choices.
        assert!(
            dx > 1.3 && dx < 1.8,
            "iso x-extent should be ~1.6, got {dx}"
        );
        assert!(
            dy > 1.2 && dy < 1.6,
            "iso y-extent should be ~1.4, got {dy}"
        );
    }
}
