//! Standard scientific colormaps for visualising scalar fields.
//!
//! No image library, no LUT files — each colormap is a small
//! piecewise-linear interpolation between RGB control points. Good
//! enough for residual charts, surface contour previews, and the
//! mesh-wireframe overlay in the viewport. Production-quality
//! colour-bar rendering (with proper perceptual viridis, gamma
//! correction, log-scale support, NaN handling) lands when we add
//! the wgpu surface-rendering pass; for now the cool-to-warm ramp
//! covers most "what does this field look like?" questions.

/// Map a normalised value `t ∈ [0, 1]` to an RGB colour using a
/// 5-stop **cool-to-warm** divergent ramp:
///
/// ```text
/// t = 0.00 → dark blue   (38, 41, 86)
/// t = 0.25 → light blue  (88, 162, 215)
/// t = 0.50 → near-white  (242, 242, 240)
/// t = 0.75 → orange      (242, 168, 90)
/// t = 1.00 → dark red    (165, 30, 30)
/// ```
///
/// Values outside `[0, 1]` are clamped — useful when the caller has
/// already remapped a field to its `[min, max]` range and just wants
/// to colour everything sensibly without special-casing the endpoints.
pub fn cool_to_warm(t: f32) -> [u8; 3] {
    const STOPS: &[(f32, [f32; 3])] = &[
        (0.00, [38.0, 41.0, 86.0]),
        (0.25, [88.0, 162.0, 215.0]),
        (0.50, [242.0, 242.0, 240.0]),
        (0.75, [242.0, 168.0, 90.0]),
        (1.00, [165.0, 30.0, 30.0]),
    ];
    let t = t.clamp(0.0, 1.0);
    // Find the bracket [STOPS[i], STOPS[i+1]] containing t.
    for i in 0..STOPS.len() - 1 {
        let (t0, c0) = STOPS[i];
        let (t1, c1) = STOPS[i + 1];
        if t <= t1 {
            let span = (t1 - t0).max(1e-9);
            let f = (t - t0) / span;
            return [
                lerp(c0[0], c1[0], f) as u8,
                lerp(c0[1], c1[1], f) as u8,
                lerp(c0[2], c1[2], f) as u8,
            ];
        }
    }
    // t == 1.0 (or > 1.0 after the bracket loop misses) — last stop.
    let (_, c) = STOPS[STOPS.len() - 1];
    [c[0] as u8, c[1] as u8, c[2] as u8]
}

/// Map a value `v ∈ [min, max]` to an RGB colour via [`cool_to_warm`],
/// normalising first. Degenerate ranges (`min == max`) collapse the
/// whole field to the midpoint colour rather than dividing by zero.
pub fn cool_to_warm_in_range(v: f64, min: f64, max: f64) -> [u8; 3] {
    let span = max - min;
    let t = if span.abs() < 1e-30 {
        0.5
    } else {
        ((v - min) / span) as f32
    };
    cool_to_warm(t)
}

fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cool_to_warm_endpoints_match_stops() {
        // t=0 → dark blue, t=1 → dark red. Verify exact stop colours
        // come back so the colour bar's labels match what the user sees.
        assert_eq!(cool_to_warm(0.0), [38, 41, 86]);
        assert_eq!(cool_to_warm(1.0), [165, 30, 30]);
    }

    #[test]
    fn cool_to_warm_midpoint_is_near_white() {
        // The 0.5 stop is the divergent map's neutral point. Should
        // be visually close to white so users instantly recognise
        // "this is the average / no signal" region.
        let mid = cool_to_warm(0.5);
        assert!(mid[0] > 200 && mid[1] > 200 && mid[2] > 200, "got {mid:?}");
    }

    #[test]
    fn cool_to_warm_clamps_out_of_range() {
        // Values below 0 collapse to the first stop colour; values
        // above 1 collapse to the last. No panics, no garbage RGB.
        assert_eq!(cool_to_warm(-5.0), cool_to_warm(0.0));
        assert_eq!(cool_to_warm(7.0), cool_to_warm(1.0));
        assert_eq!(cool_to_warm(f32::NEG_INFINITY), cool_to_warm(0.0));
        assert_eq!(cool_to_warm(f32::INFINITY), cool_to_warm(1.0));
    }

    #[test]
    fn cool_to_warm_is_monotonic_in_brightness_at_endpoints() {
        // Going from t=0 (dark blue) to t=0.5 (near white), at least
        // one channel must increase strictly. Same for t=0.5 → t=1
        // (going to dark red, so the red channel should drop and
        // blue/green should drop too — divergent map crests at the
        // mid). We check the "cool half" here.
        let a = cool_to_warm(0.0);
        let b = cool_to_warm(0.5);
        // At least one channel increased when moving toward white.
        assert!(a[0] < b[0] || a[1] < b[1] || a[2] < b[2]);
    }

    #[test]
    fn cool_to_warm_in_range_handles_degenerate_range() {
        // min == max → divisor is zero. Don't panic; return the
        // mid colour so the user sees "uniform field" cleanly.
        let constant_field = cool_to_warm_in_range(7.0, 7.0, 7.0);
        let mid = cool_to_warm(0.5);
        assert_eq!(constant_field, mid);
    }

    #[test]
    fn cool_to_warm_in_range_normalises_correctly() {
        // A field ranging 100..200, sampled at 100, should give
        // the cold endpoint; at 200 the warm; at 150 the middle.
        assert_eq!(
            cool_to_warm_in_range(100.0, 100.0, 200.0),
            cool_to_warm(0.0)
        );
        assert_eq!(
            cool_to_warm_in_range(200.0, 100.0, 200.0),
            cool_to_warm(1.0)
        );
        assert_eq!(
            cool_to_warm_in_range(150.0, 100.0, 200.0),
            cool_to_warm(0.5)
        );
    }
}
