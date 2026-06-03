//! 2D Delaunay triangulation via the Bowyer-Watson incremental
//! algorithm.

use crate::error::CgalError;

/// 2D Delaunay triangulation.
///
/// Returns a list of triangles as `(i, j, k)` index triples into the
/// input `points` array (CCW). Requires ≥ 3 points and rejects fully
/// collinear input.
pub fn triangulate_2d(points: &[[f64; 2]]) -> Result<Vec<(usize, usize, usize)>, CgalError> {
    if points.len() < 3 {
        return Err(CgalError::NotEnoughPoints {
            needed: 3,
            given: points.len(),
        });
    }
    // Reject collinear input early.
    let mut collinear = true;
    let p0 = points[0];
    let p1 = points[1];
    for p in points.iter().skip(2) {
        let cross = (p1[0] - p0[0]) * (p[1] - p0[1]) - (p1[1] - p0[1]) * (p[0] - p0[0]);
        if cross.abs() > 1e-12 {
            collinear = false;
            break;
        }
    }
    if collinear {
        return Err(CgalError::Degenerate);
    }
    // Compute super-triangle covering all points (vertex indices are
    // `points.len() + 0/1/2`).
    let (xmin, xmax, ymin, ymax) = aabb(points);
    let dx = xmax - xmin;
    let dy = ymax - ymin;
    let mid_x = (xmin + xmax) / 2.0;
    let mid_y = (ymin + ymax) / 2.0;
    let m = dx.max(dy) * 20.0;
    let st0 = [mid_x - m, mid_y - m];
    let st1 = [mid_x + m, mid_y - m];
    let st2 = [mid_x, mid_y + m];
    let mut work: Vec<[f64; 2]> = points.to_vec();
    work.push(st0);
    work.push(st1);
    work.push(st2);
    let i0 = points.len();
    let i1 = points.len() + 1;
    let i2 = points.len() + 2;
    let mut tris: Vec<(usize, usize, usize)> = vec![(i0, i1, i2)];

    for (i, p) in points.iter().enumerate() {
        // Find every triangle whose circumcircle contains p.
        let bad: Vec<usize> = tris
            .iter()
            .enumerate()
            .filter(|(_, t)| in_circumcircle(work[t.0], work[t.1], work[t.2], *p))
            .map(|(idx, _)| idx)
            .collect();
        // Polygon hole = edges of bad tris that aren't shared with
        // another bad tri.
        let mut edge_count: std::collections::HashMap<(usize, usize), i32> =
            std::collections::HashMap::new();
        for &bi in &bad {
            let t = tris[bi];
            for &(a, b) in &[(t.0, t.1), (t.1, t.2), (t.2, t.0)] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        // Remove bad tris (reverse so indices stay valid).
        for &bi in bad.iter().rev() {
            tris.swap_remove(bi);
        }
        // Re-stitch with p.
        for ((a, b), c) in edge_count.iter() {
            if *c == 1 {
                tris.push((*a, *b, i));
            }
        }
    }
    // Drop triangles touching the super-triangle vertices.
    let final_tris = tris
        .into_iter()
        .filter(|t| t.0 < points.len() && t.1 < points.len() && t.2 < points.len())
        .collect();
    Ok(final_tris)
}

fn aabb(pts: &[[f64; 2]]) -> (f64, f64, f64, f64) {
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for p in pts {
        if p[0] < xmin {
            xmin = p[0];
        }
        if p[0] > xmax {
            xmax = p[0];
        }
        if p[1] < ymin {
            ymin = p[1];
        }
        if p[1] > ymax {
            ymax = p[1];
        }
    }
    (xmin, xmax, ymin, ymax)
}

fn in_circumcircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], p: [f64; 2]) -> bool {
    // The in-circle determinant
    //   |a-p|² ((b-p)×(c-p)) + |b-p|² ((c-p)×(a-p)) + |c-p|² ((a-p)×(b-p))
    // is positive when `p` is inside the circumcircle **only for a
    // CCW-wound triangle**; for a CW triangle the sign flips. The
    // Bowyer-Watson re-stitch builds triangles from sorted edge keys,
    // so their winding is arbitrary — multiply the determinant by the
    // triangle's signed orientation to get a winding-independent test.
    let ax = a[0] - p[0];
    let ay = a[1] - p[1];
    let bx = b[0] - p[0];
    let by = b[1] - p[1];
    let cx = c[0] - p[0];
    let cy = c[1] - p[1];
    let det = (ax * ax + ay * ay) * (bx * cy - cx * by)
        - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay);
    // Signed area ×2 of triangle (a, b, c): >0 CCW, <0 CW.
    let orient = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    det * orient > 0.0
}
