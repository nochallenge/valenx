//! Alpha shapes — concave hull family.
//!
//! v1 approximation: run [`crate::delaunay::triangulate_2d`], then
//! keep only triangles whose circumradius < `alpha`. The boundary
//! edges of the kept set form the alpha-shape boundary. Returned as
//! an ordered polyline (caller can close if needed).

use crate::delaunay::triangulate_2d;
use crate::error::CgalError;

/// Compute the 2D alpha-shape boundary of `points`.
pub fn alpha_shape_2d(points: &[[f64; 2]], alpha: f64) -> Result<Vec<[f64; 2]>, CgalError> {
    if alpha <= 0.0 {
        return Err(CgalError::BadParameter {
            name: "alpha",
            reason: format!("must be > 0, got {alpha}"),
        });
    }
    let tris = triangulate_2d(points)?;
    let mut edge_count: std::collections::HashMap<(usize, usize), i32> =
        std::collections::HashMap::new();
    for t in &tris {
        let r = circumradius(points[t.0], points[t.1], points[t.2]);
        if r > alpha {
            continue;
        }
        for &(a, b) in &[(t.0, t.1), (t.1, t.2), (t.2, t.0)] {
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    // Boundary = edges with count == 1.
    let mut boundary_edges: Vec<(usize, usize)> = edge_count
        .into_iter()
        .filter_map(|(k, c)| if c == 1 { Some(k) } else { None })
        .collect();
    if boundary_edges.is_empty() {
        return Err(CgalError::Numerical {
            algo: "alpha_shape_2d",
            reason: "no triangles below alpha threshold".into(),
        });
    }
    // Walk the edges into an ordered polyline. Best-effort — if the
    // boundary isn't a single closed loop the output is partial.
    let mut out = Vec::new();
    let first = boundary_edges.swap_remove(0);
    out.push(points[first.0]);
    let mut current = first.1;
    let mut safety = boundary_edges.len() + 1;
    loop {
        out.push(points[current]);
        let pos = boundary_edges.iter().position(|(a, b)| *a == current || *b == current);
        match pos {
            Some(i) => {
                let (a, b) = boundary_edges.swap_remove(i);
                current = if a == current { b } else { a };
                if current == first.0 {
                    break;
                }
            }
            None => break,
        }
        safety -= 1;
        if safety == 0 {
            break;
        }
    }
    Ok(out)
}

fn circumradius(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let ax = b[0] - a[0];
    let ay = b[1] - a[1];
    let bx = c[0] - a[0];
    let by = c[1] - a[1];
    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-18 {
        return f64::INFINITY;
    }
    let ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    let uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    (ux * ux + uy * uy).sqrt()
}
