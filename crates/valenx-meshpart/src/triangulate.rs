//! Ear-clipping triangulation for a 2D simple polygon.

use crate::error::MeshPartError;

/// Triangulate a 2D simple polygon `polygon` via ear-clipping.
///
/// `polygon` is the closed loop of `[u, v]` coordinates (the last
/// point is *not* an explicit duplicate of the first). Returns one
/// `(i, j, k)` triple per triangle indexing back into `polygon`.
///
/// Errors:
/// - [`MeshPartError::Empty`] when fewer than 3 points are supplied.
/// - [`MeshPartError::BadPolygon`] when the polygon is non-simple
///   (self-intersecting) or degenerate so the ear-clipping loop
///   bottoms out before producing all triangles.
pub fn triangulate_polygon(polygon: &[[f64; 2]]) -> Result<Vec<(usize, usize, usize)>, MeshPartError> {
    if polygon.len() < 3 {
        return Err(MeshPartError::Empty("polygon (<3 points)"));
    }
    let n = polygon.len();
    // Orientation: positive area → CCW.
    let signed = signed_area(polygon);
    // Indices of remaining vertices (CCW order).
    let mut indices: Vec<usize> = if signed >= 0.0 {
        (0..n).collect()
    } else {
        (0..n).rev().collect()
    };

    let mut out: Vec<(usize, usize, usize)> = Vec::with_capacity(n.saturating_sub(2));
    let mut guard = n * n; // Hard upper-bound on iterations.

    while indices.len() > 3 {
        if guard == 0 {
            return Err(MeshPartError::BadPolygon(
                "ear-clipping failed to terminate (non-simple polygon?)".into(),
            ));
        }
        guard -= 1;

        let m = indices.len();
        let mut ear = None;
        for i in 0..m {
            let prev = indices[(i + m - 1) % m];
            let cur = indices[i];
            let nxt = indices[(i + 1) % m];
            let a = polygon[prev];
            let b = polygon[cur];
            let c = polygon[nxt];
            if !is_convex(a, b, c) {
                continue;
            }
            // No other vertex should sit inside the (a, b, c) triangle.
            let mut clean = true;
            for &j in &indices {
                if j == prev || j == cur || j == nxt {
                    continue;
                }
                if point_in_triangle(polygon[j], a, b, c) {
                    clean = false;
                    break;
                }
            }
            if clean {
                ear = Some((prev, cur, nxt, i));
                break;
            }
        }
        let Some((p, c, nx, i)) = ear else {
            return Err(MeshPartError::BadPolygon("no ear found".into()));
        };
        out.push((p, c, nx));
        indices.remove(i);
    }
    if indices.len() == 3 {
        out.push((indices[0], indices[1], indices[2]));
    }
    Ok(out)
}

fn signed_area(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = p[i];
        let b = p[(i + 1) % n];
        s += (b[0] - a[0]) * (b[1] + a[1]);
    }
    -0.5 * s // negative because we accumulated the shoelace dual
}

fn is_convex(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    cross > 0.0
}

fn point_in_triangle(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let s1 = sign(p, a, b);
    let s2 = sign(p, b, c);
    let s3 = sign(p, c, a);
    let neg = (s1 < 0.0) || (s2 < 0.0) || (s3 < 0.0);
    let pos = (s1 > 0.0) || (s2 > 0.0) || (s3 > 0.0);
    !(neg && pos)
}

fn sign(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    (p[0] - b[0]) * (a[1] - b[1]) - (a[0] - b[0]) * (p[1] - b[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_errors() {
        assert!(matches!(
            triangulate_polygon(&[]),
            Err(MeshPartError::Empty(_))
        ));
    }

    #[test]
    fn triangle_passes_through() {
        let tri = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let t = triangulate_polygon(&tri).unwrap();
        assert_eq!(t.len(), 1);
    }

    #[test]
    fn square_gives_two_triangles() {
        let sq = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let t = triangulate_polygon(&sq).unwrap();
        assert_eq!(t.len(), 2);
    }

    #[test]
    fn cw_input_is_flipped_to_ccw() {
        let sq = vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]; // CW
        let t = triangulate_polygon(&sq).unwrap();
        assert_eq!(t.len(), 2);
    }
}
