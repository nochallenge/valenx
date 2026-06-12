//! Geometric closed-form calculations over 2D drafting entities.
//!
//! A pure-math layer computing derived scalars (arc length, …) from raw
//! vertex coordinates via classical formulas.

/// Total length of a polyline — the sum of the Euclidean distances between consecutive
/// vertices, `Σ √((xᵢ₊₁−xᵢ)² + (yᵢ₊₁−yᵢ)²)`. The closing segment of a closed polyline is NOT
/// included (closure is a topology flag, not part of the open path). Returns `0.0` for fewer
/// than two vertices. Segments with a non-finite endpoint are skipped so a stray `NaN` cannot
/// poison the sum.
pub fn polyline_length(vertices: &[[f64; 2]]) -> f64 {
    vertices
        .windows(2)
        .filter_map(|w| {
            let [x0, y0] = w[0];
            let [x1, y1] = w[1];
            if x0.is_finite() && y0.is_finite() && x1.is_finite() && y1.is_finite() {
                let dx = x1 - x0;
                let dy = y1 - y0;
                Some((dx * dx + dy * dy).sqrt())
            } else {
                None
            }
        })
        .sum()
}

/// Enclosed area of a polygon via the shoelace formula, `½·|Σ (xᵢ·yᵢ₊₁ − xᵢ₊₁·yᵢ)|` with
/// cyclic closure (the last vertex wraps to the first). Always non-negative (orientation-
/// independent). Returns `0.0` for fewer than three vertices or if any coordinate is non-finite.
pub fn polygon_area(vertices: &[[f64; 2]]) -> f64 {
    if vertices.len() < 3
        || !vertices
            .iter()
            .all(|v| v[0].is_finite() && v[1].is_finite())
    {
        return 0.0;
    }
    let n = vertices.len();
    let sum: f64 = vertices
        .iter()
        .enumerate()
        .map(|(i, &[xi, yi])| {
            let [xj, yj] = vertices[(i + 1) % n];
            xi * yj - xj * yi
        })
        .sum();
    sum.abs() / 2.0
}

/// Diagonal of the axis-aligned bounding box over the vertices,
/// `√((xmax−xmin)² + (ymax−ymin)²)`. Translation-invariant; distinct from polyline_length
/// (perimeter) and polygon_area. Returns `0.0` for fewer than two vertices or if any coordinate
/// is non-finite.
pub fn bounding_box_diagonal_2d(vertices: &[[f64; 2]]) -> f64 {
    if vertices.len() < 2
        || !vertices
            .iter()
            .all(|v| v[0].is_finite() && v[1].is_finite())
    {
        return 0.0;
    }
    let (mut xmin, mut ymin) = (vertices[0][0], vertices[0][1]);
    let (mut xmax, mut ymax) = (xmin, ymin);
    for &[x, y] in vertices {
        xmin = xmin.min(x);
        ymin = ymin.min(y);
        xmax = xmax.max(x);
        ymax = ymax.max(y);
    }
    let dx = xmax - xmin;
    let dy = ymax - ymin;
    (dx * dx + dy * dy).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polyline_length_worked() {
        // Two points (0,0)→(3,4): the 3-4-5 hypotenuse = 5.0.
        assert!((polyline_length(&[[0.0, 0.0], [3.0, 4.0]]) - 5.0).abs() < 1e-9);
        // Open square path (0,0)→(1,0)→(1,1)→(0,1): 1+1+1 = 3.0 (no closing segment).
        assert!(
            (polyline_length(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]) - 3.0).abs() < 1e-9
        );
        // Right-triangle legs (0,0)→(3,0)→(3,4): 3 + 4 = 7.0.
        assert!((polyline_length(&[[0.0, 0.0], [3.0, 0.0], [3.0, 4.0]]) - 7.0).abs() < 1e-9);
        // Degenerate: 0 or 1 vertex → 0.0.
        assert_eq!(polyline_length(&[]), 0.0);
        assert_eq!(polyline_length(&[[5.0, 5.0]]), 0.0);
        // A non-finite endpoint skips only the touching segments (the first leg, 1.0, survives).
        assert!(
            (polyline_length(&[[0.0, 0.0], [1.0, 0.0], [f64::NAN, 1.0], [1.0, 1.0]]) - 1.0).abs()
                < 1e-9
        );
    }

    #[test]
    fn polygon_area_worked() {
        // Unit square → 1.0; 2×3 rectangle → 6.0; right triangle (0,0),(3,0),(0,4) → 6.0.
        assert!(
            (polygon_area(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]) - 1.0).abs() < 1e-9
        );
        assert!(
            (polygon_area(&[[0.0, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]]) - 6.0).abs() < 1e-9
        );
        assert!((polygon_area(&[[0.0, 0.0], [3.0, 0.0], [0.0, 4.0]]) - 6.0).abs() < 1e-9);
        // Orientation-independent: reversing the winding gives the same positive area.
        let ccw = polygon_area(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
        let cw = polygon_area(&[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0]]);
        assert!((ccw - cw).abs() < 1e-9 && ccw > 0.0);
        // Collinear → 0.0; <3 verts → 0.0; non-finite → 0.0.
        assert_eq!(polygon_area(&[[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]), 0.0);
        assert_eq!(polygon_area(&[[0.0, 0.0], [1.0, 0.0]]), 0.0);
        assert_eq!(
            polygon_area(&[[0.0, 0.0], [f64::NAN, 0.0], [0.0, 1.0]]),
            0.0
        );
    }

    #[test]
    fn bounding_box_diagonal_2d_worked() {
        // Extent (3,4) → diagonal 5.0.
        assert!((bounding_box_diagonal_2d(&[[0.0, 0.0], [3.0, 4.0]]) - 5.0).abs() < 1e-9);
        assert!(
            (bounding_box_diagonal_2d(&[[0.0, 0.0], [3.0, 0.0], [0.0, 4.0]]) - 5.0).abs() < 1e-9
        );
        // Translation-invariant: shifting all vertices leaves the diagonal unchanged.
        assert!(
            (bounding_box_diagonal_2d(&[[10.0, 20.0], [13.0, 24.0]])
                - bounding_box_diagonal_2d(&[[0.0, 0.0], [3.0, 4.0]]))
            .abs()
                < 1e-9
        );
        // Collinear horizontal extent 10 → 10.0.
        assert!(
            (bounding_box_diagonal_2d(&[[0.0, 0.0], [5.0, 0.0], [10.0, 0.0]]) - 10.0).abs() < 1e-9
        );
        // <2 verts → 0.0; empty → 0.0; non-finite → 0.0.
        assert_eq!(bounding_box_diagonal_2d(&[[5.0, 5.0]]), 0.0);
        assert_eq!(bounding_box_diagonal_2d(&[]), 0.0);
        assert_eq!(
            bounding_box_diagonal_2d(&[[f64::NAN, 0.0], [1.0, 1.0]]),
            0.0
        );
    }
}
