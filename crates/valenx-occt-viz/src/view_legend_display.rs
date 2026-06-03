//! Phase 198 — `Prs3d_LegendItem` — display the field-mapping legend
//! (colour → value, units).
//!
//! ## What OCCT does
//!
//! When a scalar-field overlay is active (typically a finite-element
//! stress / temperature result colour-mapped onto the geometry),
//! `Prs3d_LegendItem` paints a vertical colour bar in the corner of
//! the viewport with tick labels showing the value at each colour
//! stop and the units (e.g. "MPa" or "°C"). The colour ramp is
//! configurable — OCCT ships rainbow + greyscale + cool-to-warm;
//! Valenx uses cool-to-warm by default (`valenx_app::viewport`'s
//! field overlay).
//!
//! ## v1 status
//!
//! **Honest v1.** Computes a [`LegendStops`] table — N equally-spaced
//! (value, RGBA) pairs from the cool-to-warm ramp matching the
//! existing `valenx_fields` palette. The caller paints the bar
//! using egui's `Painter` directly (rectangles for each stop +
//! text labels). The cool-to-warm formula is the diverging
//! Kindlmann/Moreland 2007 ramp: hue 240° at min (blue),
//! desaturated at midpoint, hue 0° at max (red).

use crate::error::OcctVizError;
use crate::prs3d_drawer_face_color::{prs3d_drawer_face_color, FaceColorRgba};

/// One stop in the legend: a value + the RGBA the renderer assigned.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct LegendStop {
    /// Scalar value this colour represents (in the field's units).
    pub value: f32,
    /// Colour the ramp produces at `value`.
    pub color: FaceColorRgba,
}

/// Full legend payload: N stops spanning `[min_value, max_value]`,
/// plus the units string.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendStops {
    /// Ordered stops from low value → high value.
    pub stops: Vec<LegendStop>,
    /// Units label rendered next to the bar (e.g. "MPa").
    pub units: String,
}

/// Build a legend with `n_stops` equally-spaced colour stops between
/// `min_value` and `max_value`.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `n_stops < 2`, `n_stops > 16`,
///   either value is non-finite, or `min_value > max_value`.
pub fn view_legend_display(
    min_value: f32,
    max_value: f32,
    n_stops: usize,
    units: &str,
) -> Result<LegendStops, OcctVizError> {
    if !min_value.is_finite() || !max_value.is_finite() {
        return Err(OcctVizError::bad_input("value_range", "must be finite"));
    }
    if min_value > max_value {
        return Err(OcctVizError::bad_input(
            "value_range",
            format!("min={min_value} > max={max_value}"),
        ));
    }
    if !(2..=16).contains(&n_stops) {
        return Err(OcctVizError::bad_input(
            "n_stops",
            format!("must be in [2, 16] (got {n_stops})"),
        ));
    }

    let mut stops = Vec::with_capacity(n_stops);
    for i in 0..n_stops {
        let t = i as f32 / (n_stops - 1) as f32;
        let value = min_value + (max_value - min_value) * t;
        let (r, g, b) = cool_to_warm(t);
        let color = prs3d_drawer_face_color(r, g, b, 1.0)?;
        stops.push(LegendStop { value, color });
    }
    Ok(LegendStops {
        stops,
        units: units.to_owned(),
    })
}

/// Diverging cool-to-warm ramp (Moreland 2007): blue → white → red.
/// `t ∈ [0, 1]`; clamps out-of-range inputs.
fn cool_to_warm(t: f32) -> (f32, f32, f32) {
    let t = t.clamp(0.0, 1.0);
    // Endpoints: blue (0.230, 0.299, 0.754), red (0.706, 0.016, 0.150)
    // Midpoint: near-white (0.865, 0.865, 0.865).
    if t < 0.5 {
        let s = t * 2.0; // 0..1 in lower half
        let r = 0.230 + (0.865 - 0.230) * s;
        let g = 0.299 + (0.865 - 0.299) * s;
        let b = 0.754 + (0.865 - 0.754) * s;
        (r, g, b)
    } else {
        let s = (t - 0.5) * 2.0;
        let r = 0.865 + (0.706 - 0.865) * s;
        let g = 0.865 + (0.016 - 0.865) * s;
        let b = 0.865 + (0.150 - 0.865) * s;
        (r, g, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_inverted_range() {
        let err = view_legend_display(10.0, 0.0, 5, "MPa").unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_one_stop() {
        let err = view_legend_display(0.0, 10.0, 1, "MPa").unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_too_many_stops() {
        let err = view_legend_display(0.0, 10.0, 100, "MPa").unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn stops_span_range() {
        let l = view_legend_display(0.0, 100.0, 5, "MPa").unwrap();
        assert_eq!(l.stops.len(), 5);
        assert!((l.stops[0].value - 0.0).abs() < 1e-4);
        assert!((l.stops[4].value - 100.0).abs() < 1e-4);
        assert!((l.stops[2].value - 50.0).abs() < 1e-4);
    }

    #[test]
    fn min_is_blue_max_is_red() {
        let l = view_legend_display(0.0, 1.0, 2, "x").unwrap();
        // First stop = blue-ish (b > r).
        assert!(l.stops[0].color.b > l.stops[0].color.r);
        // Last stop = red-ish (r > b).
        assert!(l.stops[1].color.r > l.stops[1].color.b);
    }

    #[test]
    fn units_round_trip() {
        let l = view_legend_display(0.0, 1.0, 3, "°C").unwrap();
        assert_eq!(l.units, "°C");
    }
}
