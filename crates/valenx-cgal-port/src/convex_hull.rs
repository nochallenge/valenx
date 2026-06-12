//! Convex hull algorithms — 2D Graham scan + 3D incremental hull.

use nalgebra::Vector3;

use crate::error::CgalError;

/// 2D convex hull via Graham scan.
///
/// Returns the hull vertices in CCW order. Duplicate-tolerant. Requires
/// ≥ 3 input points.
pub fn hull_2d(points: &[[f64; 2]]) -> Result<Vec<[f64; 2]>, CgalError> {
    if points.len() < 3 {
        return Err(CgalError::NotEnoughPoints {
            needed: 3,
            given: points.len(),
        });
    }
    let mut pts: Vec<[f64; 2]> = points.to_vec();
    // Lowest Y (then lowest X) is anchor — Graham invariant.
    pts.sort_by(|a, b| {
        a[1].partial_cmp(&b[1])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a[0].partial_cmp(&b[0]).unwrap_or(std::cmp::Ordering::Equal))
    });
    let anchor = pts[0];
    // Sort remainder by polar angle from anchor.
    pts[1..].sort_by(|a, b| {
        let ang_a = (a[1] - anchor[1]).atan2(a[0] - anchor[0]);
        let ang_b = (b[1] - anchor[1]).atan2(b[0] - anchor[0]);
        ang_a
            .partial_cmp(&ang_b)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut stack: Vec<[f64; 2]> = Vec::with_capacity(pts.len());
    for p in pts {
        while stack.len() >= 2 && cross_2d(stack[stack.len() - 2], stack[stack.len() - 1], p) <= 0.0
        {
            stack.pop();
        }
        stack.push(p);
    }
    if stack.len() < 3 {
        return Err(CgalError::Degenerate);
    }
    Ok(stack)
}

fn cross_2d(o: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

/// 3D convex hull via the incremental algorithm.
///
/// Returns the hull triangles as `(i, j, k)` index triples into the
/// input point list. The first 4 non-coplanar points seed a tetrahedron;
/// every subsequent point is added by removing visible faces and
/// re-stitching the horizon.
///
/// v1 returns the face list — vertex indices may have duplicates; the
/// caller can dedupe.
pub fn hull_3d(points: &[Vector3<f64>]) -> Result<Vec<(usize, usize, usize)>, CgalError> {
    if points.len() < 4 {
        return Err(CgalError::NotEnoughPoints {
            needed: 4,
            given: points.len(),
        });
    }
    // Bootstrap a tetrahedron from the first 4 non-coplanar points.
    let mut tet = [0usize; 4];
    tet[0] = 0;
    tet[1] = 1;
    tet[2] = 2;
    let mut found_3 = false;
    for (i, p) in points.iter().enumerate().take(points.len()).skip(3) {
        let a = points[tet[0]];
        let b = points[tet[1]];
        let c = points[tet[2]];
        let n = (b - a).cross(&(c - a));
        if (p - a).dot(&n).abs() > 1e-9 {
            tet[3] = i;
            found_3 = true;
            break;
        }
    }
    if !found_3 {
        return Err(CgalError::Degenerate);
    }
    // Initial 4 faces of the tetrahedron, oriented outward.
    let mut faces = tetrahedron_faces(points, tet);

    // Add the remaining points one at a time.
    for (idx, p) in points.iter().enumerate() {
        if tet.contains(&idx) {
            continue;
        }
        // Find faces visible from p.
        let visible: Vec<usize> = faces
            .iter()
            .enumerate()
            .filter(|(_, f)| {
                let a = points[f.0];
                let b = points[f.1];
                let c = points[f.2];
                let n = (b - a).cross(&(c - a));
                (p - a).dot(&n) > 1e-9
            })
            .map(|(i, _)| i)
            .collect();
        if visible.is_empty() {
            // Interior point; skip.
            continue;
        }
        // Horizon = edges on the boundary of the visible region.
        let mut edge_count: std::collections::HashMap<(usize, usize), i32> =
            std::collections::HashMap::new();
        for &vi in &visible {
            let f = faces[vi];
            for &(a, b) in &[(f.0, f.1), (f.1, f.2), (f.2, f.0)] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        // Remove visible faces (in reverse so indices stay valid).
        for &vi in visible.iter().rev() {
            faces.swap_remove(vi);
        }
        // Stitch in new faces over the horizon edges (count == 1).
        for ((a, b), c) in edge_count.iter() {
            if *c == 1 {
                // Need to find an orientation that points away from
                // the centroid of all existing hull vertices. Picking
                // ((a, b, idx)) by default; if it points inward, swap.
                let ca = points[*a];
                let cb = points[*b];
                let cp = points[idx];
                let n = (cb - ca).cross(&(cp - ca));
                // Centre estimate = avg of tet points.
                let centre =
                    (points[tet[0]] + points[tet[1]] + points[tet[2]] + points[tet[3]]) / 4.0;
                if (ca - centre).dot(&n) > 0.0 {
                    faces.push((*a, *b, idx));
                } else {
                    faces.push((*b, *a, idx));
                }
            }
        }
    }
    Ok(faces)
}

fn tetrahedron_faces(points: &[Vector3<f64>], tet: [usize; 4]) -> Vec<(usize, usize, usize)> {
    let centre = (points[tet[0]] + points[tet[1]] + points[tet[2]] + points[tet[3]]) / 4.0;
    // 4 candidate faces; orient each so its normal points outward
    // (away from `centre`).
    let raw: [(usize, usize, usize); 4] = [
        (tet[0], tet[1], tet[2]),
        (tet[0], tet[1], tet[3]),
        (tet[0], tet[2], tet[3]),
        (tet[1], tet[2], tet[3]),
    ];
    raw.iter()
        .map(|(a, b, c)| {
            let pa = points[*a];
            let pb = points[*b];
            let pc = points[*c];
            let n = (pb - pa).cross(&(pc - pa));
            if (pa - centre).dot(&n) > 0.0 {
                (*a, *b, *c)
            } else {
                (*a, *c, *b)
            }
        })
        .collect()
}
