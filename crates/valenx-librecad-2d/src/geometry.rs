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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polyline_length_worked() {
        // Two points (0,0)→(3,4): the 3-4-5 hypotenuse = 5.0.
        assert!((polyline_length(&[[0.0, 0.0], [3.0, 4.0]]) - 5.0).abs() < 1e-9);
        // Open square path (0,0)→(1,0)→(1,1)→(0,1): 1+1+1 = 3.0 (no closing segment).
        assert!(
            (polyline_length(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]) - 3.0).abs()
                < 1e-9
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
}
