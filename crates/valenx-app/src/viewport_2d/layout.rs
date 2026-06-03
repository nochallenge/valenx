//! Pure coordinate-mapping math for the 2D DNA viewport.
//!
//! No egui types — every function is testable without a render context.

use std::f32::consts::{FRAC_PI_2, TAU};

// ─────────────────────────────────────────────────────────────────────────────
// Linear track
// ─────────────────────────────────────────────────────────────────────────────

/// Pixel x-coordinate of base-pair index `bp` within a horizontal track.
///
/// `pan` is the bp index at the **left edge** of the visible window.
/// `bases_per_pixel` is the zoom level (lower = more zoomed in).
/// `track_left_px` is the pixel x of the track's left edge.
pub fn bp_to_x(bp: f32, pan: f32, bases_per_pixel: f32, track_left_px: f32) -> f32 {
    track_left_px + (bp - pan) / bases_per_pixel
}

/// Inverse: pixel x-coordinate → fractional base-pair index.
pub fn x_to_bp(x: f32, pan: f32, bases_per_pixel: f32, track_left_px: f32) -> f32 {
    (x - track_left_px) * bases_per_pixel + pan
}

/// Initial zoom so the entire sequence fits in `track_px` pixels.
///
/// Returns `1.0` for empty sequences or zero-width tracks.
pub fn fit_zoom(seq_len: usize, track_px: f32) -> f32 {
    if track_px <= 0.0 || seq_len == 0 {
        return 1.0;
    }
    (seq_len as f32 / track_px).max(0.01)
}

/// Clamp `bases_per_pixel` between `0.25` (4 px/base, very zoomed in)
/// and `2 × seq_len` (two screen-widths of sequence).
pub fn clamp_zoom(bases_per_pixel: f32, seq_len: usize) -> f32 {
    let min = 0.25;
    let max = (seq_len as f32 * 2.0).max(1.0);
    bases_per_pixel.clamp(min, max)
}

/// Clamp the pan value so the sequence doesn't scroll off both sides.
/// Allows half a screen of padding on each edge.
pub fn clamp_pan(pan: f32, seq_len: usize, bases_per_pixel: f32, track_px: f32) -> f32 {
    let half_screen = track_px * bases_per_pixel * 0.5;
    let min_pan = -half_screen;
    let max_pan = seq_len as f32 - half_screen.max(1.0);
    pan.clamp(min_pan, max_pan.max(0.0))
}

// ─────────────────────────────────────────────────────────────────────────────
// Circular map
// ─────────────────────────────────────────────────────────────────────────────

/// Angle (radians) on the plasmid ring for base-pair index `bp`.
///
/// Position 0 is at the **top** (12 o'clock). Angles increase clockwise
/// because egui's y-axis points down, so sin(angle) increases downward.
pub fn bp_to_angle(bp: usize, seq_len: usize) -> f32 {
    if seq_len == 0 {
        return -FRAC_PI_2;
    }
    TAU * (bp as f32 / seq_len as f32) - FRAC_PI_2
}

/// `(dx, dy)` pixel offset from the circle centre for `bp` on a ring of
/// `radius` pixels. Positive y points down (egui convention).
pub fn bp_to_ring_offset(bp: usize, seq_len: usize, radius: f32) -> (f32, f32) {
    let a = bp_to_angle(bp, seq_len);
    (radius * a.cos(), radius * a.sin())
}

/// Arc angle (radians) swept by a feature spanning `feature_len` bases on
/// a sequence of `seq_len` bases.
pub fn feature_arc_span(feature_len: usize, seq_len: usize) -> f32 {
    if seq_len == 0 {
        return 0.0;
    }
    TAU * (feature_len as f32 / seq_len as f32)
}

// ─────────────────────────────────────────────────────────────────────────────
// Feature colours — SnapGene / Benchling / ApE-inspired palette
// ─────────────────────────────────────────────────────────────────────────────

/// sRGB colour for a DNA feature type.
///
/// Returns `(r, g, b)` in `0..=255`. Palette mirrors SnapGene / Benchling /
/// ApE conventions so molecular biologists feel at home.
pub fn feature_rgb(feature_type: &str) -> (u8, u8, u8) {
    match feature_type.to_lowercase().as_str() {
        "cds" | "gene" => (52, 120, 246),              // blue — coding sequences
        "promoter" => (46, 186, 64),                   // green
        "terminator" => (220, 70, 70),                 // red
        "rep_origin" | "origin" => (220, 140, 20),     // amber
        "primer_bind" => (240, 200, 0),                // yellow
        "misc_feature" | "misc_rna" => (160, 120, 220), // lavender
        "regulatory" => (0, 180, 180),                 // teal
        "ltr" => (200, 100, 50),                       // burnt orange
        "rrna" | "trna" | "ncrna" => (100, 200, 200),  // cyan
        "intron" => (140, 140, 140),                   // grey
        "exon" => (80, 160, 80),                       // forest green
        _ => (130, 130, 130),                          // default grey
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ruler tick spacing
// ─────────────────────────────────────────────────────────────────────────────

/// Tick interval in base pairs so ticks are `~target_gap_px` pixels apart.
///
/// Rounds to the nearest "nice" interval (1, 2, 5, 10, 20, 50, 100, …).
pub fn nice_tick_spacing_bp(bases_per_pixel: f32, target_gap_px: f32) -> usize {
    let raw = (target_gap_px * bases_per_pixel).max(1.0);
    let magnitude = 10.0_f32.powf(raw.log10().floor());
    let normalised = raw / magnitude;
    let nice: f32 = if normalised < 1.5 {
        1.0
    } else if normalised < 3.5 {
        2.0
    } else if normalised < 7.5 {
        5.0
    } else {
        10.0
    };
    ((nice * magnitude) as usize).max(1)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bp_to_x_at_origin() {
        let x = bp_to_x(0.0, 0.0, 1.0, 100.0);
        assert!((x - 100.0).abs() < 1e-5, "got {x}");
    }

    #[test]
    fn bp_to_x_offset_from_left() {
        let x = bp_to_x(100.0, 0.0, 1.0, 0.0);
        assert!((x - 100.0).abs() < 1e-5, "got {x}");
    }

    #[test]
    fn bp_to_x_with_zoom() {
        // 2 bp/px → bp 100 should be at left + 50
        let x = bp_to_x(100.0, 0.0, 2.0, 0.0);
        assert!((x - 50.0).abs() < 1e-5, "got {x}");
    }

    #[test]
    fn x_to_bp_round_trips() {
        let pan = 50.0;
        let bpp = 2.5;
        let track_left = 20.0;
        let bp = 300.0;
        let x = bp_to_x(bp, pan, bpp, track_left);
        let back = x_to_bp(x, pan, bpp, track_left);
        assert!((back - bp).abs() < 1e-4, "round-trip: got {back}, expected {bp}");
    }

    #[test]
    fn fit_zoom_fills_track() {
        let z = fit_zoom(1000, 500.0);
        assert!((z - 2.0).abs() < 1e-5, "got {z}");
    }

    #[test]
    fn fit_zoom_zero_len_returns_one() {
        assert_eq!(fit_zoom(0, 800.0), 1.0);
    }

    #[test]
    fn fit_zoom_zero_track_returns_one() {
        assert_eq!(fit_zoom(500, 0.0), 1.0);
    }

    #[test]
    fn clamp_zoom_lower_bound() {
        assert_eq!(clamp_zoom(0.0, 1000), 0.25);
    }

    #[test]
    fn clamp_zoom_upper_bound() {
        assert_eq!(clamp_zoom(99999.0, 500), 1000.0);
    }

    #[test]
    fn bp_to_angle_zero_is_top() {
        let a = bp_to_angle(0, 1000);
        assert!((a - (-FRAC_PI_2)).abs() < 1e-5, "got {a}");
    }

    #[test]
    fn bp_to_angle_halfway_is_bottom() {
        let a = bp_to_angle(500, 1000);
        assert!((a - FRAC_PI_2).abs() < 1e-5, "got {a}");
    }

    #[test]
    fn bp_to_angle_full_wraps_to_start() {
        let a_full = bp_to_angle(1000, 1000);
        let a_zero = bp_to_angle(0, 1000);
        // Should differ by exactly TAU
        assert!((a_full - a_zero - TAU).abs() < 1e-4, "got diff {}", a_full - a_zero);
    }

    #[test]
    fn feature_arc_full_sequence() {
        let s = feature_arc_span(1000, 1000);
        assert!((s - TAU).abs() < 1e-5, "got {s}");
    }

    #[test]
    fn feature_arc_half_sequence() {
        let s = feature_arc_span(500, 1000);
        assert!((s - std::f32::consts::PI).abs() < 1e-5, "got {s}");
    }

    #[test]
    fn feature_arc_zero_seq_len() {
        assert_eq!(feature_arc_span(100, 0), 0.0);
    }

    #[test]
    fn feature_rgb_cds_is_blue_dominant() {
        let (r, _, b) = feature_rgb("CDS");
        assert!(b > r, "CDS should be blue-dominant: r={r}, b={b}");
    }

    #[test]
    fn feature_rgb_promoter_is_green_dominant() {
        let (r, g, b) = feature_rgb("promoter");
        assert!(g > r && g > b, "promoter should be green-dominant: {r},{g},{b}");
    }

    #[test]
    fn feature_rgb_unknown_is_grey() {
        let (r, g, b) = feature_rgb("unknown_xyz_abc_qrs");
        assert!(r == g && g == b, "unknown should be grey: {r},{g},{b}");
    }

    #[test]
    fn feature_rgb_case_insensitive() {
        assert_eq!(feature_rgb("CDS"), feature_rgb("cds"));
        assert_eq!(feature_rgb("PROMOTER"), feature_rgb("promoter"));
    }

    #[test]
    fn nice_tick_spacing_reasonable_at_1_bpp() {
        let sp = nice_tick_spacing_bp(1.0, 100.0);
        // 100 bp raw → rounds to 100
        assert!((50..=200).contains(&sp), "got {sp}");
    }

    #[test]
    fn nice_tick_spacing_zoomed_in() {
        // 0.1 bpp (10 px/base), target 80 px → ~8 bp → rounds to 5 or 10
        let sp = nice_tick_spacing_bp(0.1, 80.0);
        assert!((5..=20).contains(&sp), "got {sp}");
    }

    #[test]
    fn nice_tick_spacing_at_least_one() {
        assert!(nice_tick_spacing_bp(0.001, 1.0) >= 1);
    }
}
