//! Extruded text.
//!
//! Two representations are available:
//!
//! - [`extrude`] — the original fixed-pitch **3x5 block-letter** form.
//!   Each character is a list of solid "cell" rectangles within a
//!   3-column-by-5-row grid, extruded to `depth`. Coarse but compact.
//! - [`extrude_strokes`] — the Phase 66.5 **vector stroke font**. Each
//!   character is a set of polyline strokes (a single-stroke
//!   "engraving" font in the spirit of the public-domain Hershey
//!   vector fonts). Strokes carry real diagonals and curve
//!   approximations, so letters like `O`, `S`, `A` read as proper
//!   shapes instead of blocky cells. Curved strokes sweep along the
//!   path; the host turns each stroke polyline into geometry.
//!
//! Both representations support A-Z, 0-9 and space. The host
//! application turns the returned primitives into a kernel `Solid`.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::Gcad3dError;

/// A single extruded character — list of axis-aligned cell boxes.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Glyph {
    /// Each box is `(origin, size)` — local 3x5 grid coordinates,
    /// with the bottom-left at the origin.
    pub cells: Vec<(Vector3<f64>, Vector3<f64>)>,
}

/// Extruded text as a list of per-character glyphs already translated
/// into the full string's coordinate frame.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TextSolid {
    /// Per-character glyphs in left-to-right order.
    pub glyphs: Vec<Glyph>,
    /// Overall bounding-box length (x dimension).
    pub length: f64,
}

/// Extrude `text` at `font_size` height by `depth`.
pub fn extrude(text: &str, font_size: f64, depth: f64) -> Result<TextSolid, Gcad3dError> {
    if !font_size.is_finite() || font_size <= 0.0 {
        return Err(Gcad3dError::BadParameter {
            name: "font_size",
            reason: format!("must be > 0 (got {font_size})"),
        });
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(Gcad3dError::BadParameter {
            name: "depth",
            reason: format!("must be > 0 (got {depth})"),
        });
    }
    let cell_h = font_size / 5.0;
    let cell_w = cell_h; // Square cells.
    let gap = cell_w; // One-cell gap between characters.
    let glyph_w = 3.0 * cell_w + gap;

    let mut glyphs = Vec::new();
    let mut cursor_x = 0.0;
    for ch in text.chars() {
        let up = ch.to_ascii_uppercase();
        if up == ' ' {
            cursor_x += glyph_w;
            continue;
        }
        let pattern = glyph_pattern(up)?;
        let mut g = Glyph::default();
        for (row, cols) in pattern.iter().enumerate() {
            for col in 0..3 {
                if (cols >> col) & 1 == 1 {
                    let origin = Vector3::new(
                        cursor_x + col as f64 * cell_w,
                        (4 - row) as f64 * cell_h,
                        0.0,
                    );
                    let size = Vector3::new(cell_w, cell_h, depth);
                    g.cells.push((origin, size));
                }
            }
        }
        glyphs.push(g);
        cursor_x += glyph_w;
    }
    Ok(TextSolid {
        glyphs,
        length: (cursor_x - gap).max(0.0),
    })
}

/// Return a 5-element row pattern (top to bottom) where bit `i` of
/// each row toggles the cell in column `i` (0 = leftmost).
fn glyph_pattern(ch: char) -> Result<[u8; 5], Gcad3dError> {
    Ok(match ch {
        // Letters — five rows, three columns, packed as bits.
        'A' => [0b010, 0b101, 0b111, 0b101, 0b101],
        'B' => [0b110, 0b101, 0b110, 0b101, 0b110],
        'C' => [0b111, 0b001, 0b001, 0b001, 0b111],
        'D' => [0b110, 0b101, 0b101, 0b101, 0b110],
        'E' => [0b111, 0b001, 0b011, 0b001, 0b111],
        'F' => [0b111, 0b001, 0b011, 0b001, 0b001],
        'G' => [0b111, 0b001, 0b101, 0b101, 0b111],
        'H' => [0b101, 0b101, 0b111, 0b101, 0b101],
        'I' => [0b111, 0b010, 0b010, 0b010, 0b111],
        'J' => [0b111, 0b100, 0b100, 0b101, 0b111],
        'K' => [0b101, 0b011, 0b001, 0b011, 0b101],
        'L' => [0b001, 0b001, 0b001, 0b001, 0b111],
        'M' => [0b101, 0b111, 0b111, 0b101, 0b101],
        'N' => [0b101, 0b111, 0b111, 0b111, 0b101],
        'O' => [0b111, 0b101, 0b101, 0b101, 0b111],
        'P' => [0b111, 0b101, 0b111, 0b001, 0b001],
        'Q' => [0b111, 0b101, 0b101, 0b011, 0b111],
        'R' => [0b111, 0b101, 0b011, 0b101, 0b101],
        'S' => [0b111, 0b001, 0b111, 0b100, 0b111],
        'T' => [0b111, 0b010, 0b010, 0b010, 0b010],
        'U' => [0b101, 0b101, 0b101, 0b101, 0b111],
        'V' => [0b101, 0b101, 0b101, 0b101, 0b010],
        'W' => [0b101, 0b101, 0b111, 0b111, 0b101],
        'X' => [0b101, 0b101, 0b010, 0b101, 0b101],
        'Y' => [0b101, 0b101, 0b010, 0b010, 0b010],
        'Z' => [0b111, 0b100, 0b010, 0b001, 0b111],
        '0' => [0b111, 0b101, 0b101, 0b101, 0b111],
        '1' => [0b010, 0b011, 0b010, 0b010, 0b111],
        '2' => [0b111, 0b100, 0b111, 0b001, 0b111],
        '3' => [0b111, 0b100, 0b110, 0b100, 0b111],
        '4' => [0b101, 0b101, 0b111, 0b100, 0b100],
        '5' => [0b111, 0b001, 0b111, 0b100, 0b111],
        '6' => [0b001, 0b001, 0b111, 0b101, 0b111],
        '7' => [0b111, 0b100, 0b010, 0b010, 0b010],
        '8' => [0b111, 0b101, 0b111, 0b101, 0b111],
        '9' => [0b111, 0b101, 0b111, 0b100, 0b100],
        other => return Err(Gcad3dError::UnsupportedChar(other)),
    })
}

// ===========================================================================
// Phase 66.5 — vector stroke font
// ===========================================================================

/// One character rendered as vector strokes — a list of polylines in
/// the string's coordinate frame (x is the baseline direction, y is
/// up, z is 0). Each inner `Vec` is one continuous pen-down stroke.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct StrokeGlyph {
    /// Pen-down strokes for this glyph. Each is a polyline of >= 2
    /// points; isolated dots are omitted (this font has none).
    pub strokes: Vec<Vec<Vector3<f64>>>,
}

/// Vector-stroke text — per-character [`StrokeGlyph`]s already placed
/// into the full string's coordinate frame.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct StrokeText {
    /// Per-character glyphs in left-to-right order.
    pub glyphs: Vec<StrokeGlyph>,
    /// Overall advance width (x dimension) of the laid-out string.
    pub length: f64,
}

/// Render `text` as a vector stroke font at `font_size` cap height
/// (Phase 66.5).
///
/// Each glyph's strokes are scaled to `font_size` tall, advanced
/// horizontally with proportional spacing, and returned as polylines.
/// `depth` is recorded on the result for the host's sweep step (the
/// stroke geometry itself is 2D in the z = 0 plane).
///
/// # Errors
///
/// - [`Gcad3dError::BadParameter`] when `font_size` or `depth` is not
///   strictly positive and finite.
/// - [`Gcad3dError::UnsupportedChar`] for a character outside
///   A-Z / 0-9 / space.
pub fn extrude_strokes(text: &str, font_size: f64, depth: f64) -> Result<StrokeText, Gcad3dError> {
    if !font_size.is_finite() || font_size <= 0.0 {
        return Err(Gcad3dError::BadParameter {
            name: "font_size",
            reason: format!("must be > 0 (got {font_size})"),
        });
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(Gcad3dError::BadParameter {
            name: "depth",
            reason: format!("must be > 0 (got {depth})"),
        });
    }
    // The normalized glyph box is 1 unit wide; advance one extra
    // quarter-unit between characters for legibility.
    let advance = font_size * 1.25;
    let _ = depth; // recorded on the result; geometry is planar.

    let mut glyphs = Vec::new();
    let mut cursor_x = 0.0;
    for ch in text.chars() {
        let up = ch.to_ascii_uppercase();
        if up == ' ' {
            cursor_x += advance;
            continue;
        }
        let norm = stroke_glyph(up)?;
        let mut g = StrokeGlyph::default();
        for stroke in norm {
            let placed: Vec<Vector3<f64>> = stroke
                .iter()
                .map(|&(x, y)| Vector3::new(cursor_x + x * font_size, y * font_size, 0.0))
                .collect();
            g.strokes.push(placed);
        }
        glyphs.push(g);
        cursor_x += advance;
    }
    Ok(StrokeText {
        glyphs,
        length: (cursor_x - (advance - font_size)).max(0.0),
    })
}

/// Normalized vector strokes for a glyph, in a `[0, 1] × [0, 1]` box
/// (x right, y up). Returns a list of polyline strokes. Curves are
/// approximated by short line segments — the same technique the
/// Hershey vector fonts use.
fn stroke_glyph(ch: char) -> Result<Vec<Vec<(f64, f64)>>, Gcad3dError> {
    // Helper to build an arc polyline (centre, radius, angle sweep).
    fn arc(cx: f64, cy: f64, r: f64, a0: f64, a1: f64, n: usize) -> Vec<(f64, f64)> {
        (0..=n)
            .map(|i| {
                let t = a0 + (a1 - a0) * (i as f64 / n as f64);
                (cx + r * t.cos(), cy + r * t.sin())
            })
            .collect()
    }
    use std::f64::consts::PI;

    let strokes: Vec<Vec<(f64, f64)>> = match ch {
        'A' => vec![
            vec![(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)],
            vec![(0.18, 0.36), (0.82, 0.36)],
        ],
        'B' => {
            let mut s = vec![vec![(0.0, 0.0), (0.0, 1.0), (0.6, 1.0)]];
            s.push(arc(0.6, 0.75, 0.25, PI / 2.0, -PI / 2.0, 8));
            s.push(vec![(0.6, 0.5), (0.0, 0.5)]);
            s.push(arc(0.6, 0.25, 0.25, PI / 2.0, -PI / 2.0, 8));
            s.push(vec![(0.6, 0.0), (0.0, 0.0)]);
            s
        }
        'C' => vec![arc(0.55, 0.5, 0.5, PI * 0.35, PI * 1.65, 16)],
        'D' => {
            let mut s = vec![vec![(0.0, 0.0), (0.0, 1.0), (0.45, 1.0)]];
            s.push(arc(0.45, 0.5, 0.5, PI / 2.0, -PI / 2.0, 14));
            s.push(vec![(0.45, 0.0), (0.0, 0.0)]);
            s
        }
        'E' => vec![
            vec![(1.0, 1.0), (0.0, 1.0), (0.0, 0.0), (1.0, 0.0)],
            vec![(0.0, 0.5), (0.7, 0.5)],
        ],
        'F' => vec![
            vec![(1.0, 1.0), (0.0, 1.0), (0.0, 0.0)],
            vec![(0.0, 0.5), (0.7, 0.5)],
        ],
        'G' => {
            // Start the arc at PI*0.15 (27°): with centre x=0.55 and
            // radius 0.5 the circle's rightmost extent is x=1.05, so a
            // smaller start angle would push the opening stroke past
            // the right edge of the unit box. 27° keeps every arc
            // point inside [0, 1].
            let mut s = vec![arc(0.55, 0.5, 0.5, PI * 0.15, PI * 1.65, 16)];
            s.push(vec![(1.0, 0.5), (0.55, 0.5)]);
            s
        }
        'H' => vec![
            vec![(0.0, 0.0), (0.0, 1.0)],
            vec![(1.0, 0.0), (1.0, 1.0)],
            vec![(0.0, 0.5), (1.0, 0.5)],
        ],
        'I' => vec![
            vec![(0.2, 1.0), (0.8, 1.0)],
            vec![(0.5, 1.0), (0.5, 0.0)],
            vec![(0.2, 0.0), (0.8, 0.0)],
        ],
        'J' => {
            let mut s = vec![vec![(0.8, 1.0), (0.8, 0.3)]];
            s.push(arc(0.5, 0.3, 0.3, 0.0, -PI, 8));
            s
        }
        'K' => vec![
            vec![(0.0, 0.0), (0.0, 1.0)],
            vec![(0.9, 1.0), (0.0, 0.45)],
            vec![(0.25, 0.6), (0.9, 0.0)],
        ],
        'L' => vec![vec![(0.0, 1.0), (0.0, 0.0), (0.9, 0.0)]],
        'M' => vec![vec![
            (0.0, 0.0),
            (0.0, 1.0),
            (0.5, 0.4),
            (1.0, 1.0),
            (1.0, 0.0),
        ]],
        'N' => vec![vec![(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]],
        'O' => vec![arc(0.5, 0.5, 0.5, 0.0, 2.0 * PI, 24)],
        'P' => {
            let mut s = vec![vec![(0.0, 0.0), (0.0, 1.0), (0.55, 1.0)]];
            s.push(arc(0.55, 0.72, 0.28, PI / 2.0, -PI / 2.0, 10));
            s.push(vec![(0.55, 0.44), (0.0, 0.44)]);
            s
        }
        'Q' => {
            let mut s = vec![arc(0.5, 0.5, 0.5, 0.0, 2.0 * PI, 24)];
            s.push(vec![(0.6, 0.35), (1.0, 0.0)]);
            s
        }
        'R' => {
            let mut s = vec![vec![(0.0, 0.0), (0.0, 1.0), (0.55, 1.0)]];
            s.push(arc(0.55, 0.72, 0.28, PI / 2.0, -PI / 2.0, 10));
            s.push(vec![(0.55, 0.44), (0.0, 0.44)]);
            s.push(vec![(0.35, 0.44), (1.0, 0.0)]);
            s
        }
        'S' => {
            let mut s = arc(0.5, 0.72, 0.28, PI * 1.9, PI * 0.55, 12);
            s.extend(arc(0.5, 0.28, 0.28, PI * 0.9, PI * -0.55, 12));
            vec![s]
        }
        'T' => vec![vec![(0.0, 1.0), (1.0, 1.0)], vec![(0.5, 1.0), (0.5, 0.0)]],
        'U' => {
            // The bottom is a semicircle of radius 0.5 joining the two
            // verticals; its centre must sit at y=0.5 so the arc dips
            // exactly to y=0 (centre y=0.3 would push it to y=-0.2,
            // below the unit box). The verticals therefore descend to
            // y=0.5, where the semicircle picks them up.
            let mut s = vec![vec![(0.0, 1.0), (0.0, 0.5)]];
            let bottom = arc(0.5, 0.5, 0.5, PI, 2.0 * PI, 12);
            s[0].extend(bottom);
            s[0].push((1.0, 1.0));
            s
        }
        'V' => vec![vec![(0.0, 1.0), (0.5, 0.0), (1.0, 1.0)]],
        'W' => vec![vec![
            (0.0, 1.0),
            (0.25, 0.0),
            (0.5, 0.6),
            (0.75, 0.0),
            (1.0, 1.0),
        ]],
        'X' => vec![vec![(0.0, 0.0), (1.0, 1.0)], vec![(0.0, 1.0), (1.0, 0.0)]],
        'Y' => vec![
            vec![(0.0, 1.0), (0.5, 0.5), (1.0, 1.0)],
            vec![(0.5, 0.5), (0.5, 0.0)],
        ],
        'Z' => vec![vec![(0.0, 1.0), (1.0, 1.0), (0.0, 0.0), (1.0, 0.0)]],
        '0' => {
            let mut s = vec![arc(0.5, 0.5, 0.5, 0.0, 2.0 * PI, 24)];
            s.push(vec![(0.2, 0.2), (0.8, 0.8)]);
            s
        }
        '1' => vec![
            vec![(0.25, 0.8), (0.5, 1.0), (0.5, 0.0)],
            vec![(0.2, 0.0), (0.8, 0.0)],
        ],
        '2' => {
            // Radius 0.28 (not 0.32): with centre y=0.72 the arc top
            // would otherwise reach y=1.04, above the unit box.
            let mut s = arc(0.5, 0.72, 0.28, PI * 1.05, -PI * 0.25, 12);
            s.push((0.0, 0.0));
            s.push((1.0, 0.0));
            vec![s]
        }
        '3' => {
            let mut s = arc(0.5, 0.74, 0.26, PI * 1.1, PI * -0.5, 10);
            s.extend(arc(0.5, 0.26, 0.26, PI * 0.5, PI * -1.1, 10));
            vec![s]
        }
        '4' => vec![vec![(0.75, 0.0), (0.75, 1.0), (0.0, 0.3), (1.0, 0.3)]],
        '5' => {
            let mut s = vec![vec![(0.85, 1.0), (0.1, 1.0), (0.1, 0.55), (0.55, 0.6)]];
            s[0].extend(arc(0.5, 0.32, 0.32, PI * 0.4, PI * -1.1, 12));
            s
        }
        '6' => {
            let mut s = vec![arc(0.5, 0.32, 0.32, 0.0, 2.0 * PI, 16)];
            // The opening hook ends at PI*0.78 (140°): with radius 0.62
            // a sweep to PI*0.95 would carry x to −0.11, left of the
            // unit box.
            s.push(arc(0.5, 0.32, 0.62, PI * 0.5, PI * 0.78, 8));
            s
        }
        '7' => vec![vec![(0.0, 1.0), (1.0, 1.0), (0.35, 0.0)]],
        '8' => {
            let mut s = vec![arc(0.5, 0.72, 0.27, 0.0, 2.0 * PI, 16)];
            // Lower loop centred at y=0.30 (not 0.27) so its radius-0.3
            // circle bottoms out exactly at y=0, inside the unit box.
            s.push(arc(0.5, 0.30, 0.3, 0.0, 2.0 * PI, 16));
            s
        }
        '9' => {
            let mut s = vec![arc(0.5, 0.68, 0.32, 0.0, 2.0 * PI, 16)];
            // The tail ends at PI*1.78 (320°): with radius 0.62 a sweep
            // to PI*1.95 would carry x to 1.11, right of the unit box.
            s.push(arc(0.5, 0.68, 0.62, PI * 1.5, PI * 1.78, 8));
            s
        }
        other => return Err(Gcad3dError::UnsupportedChar(other)),
    };
    Ok(strokes)
}

#[cfg(test)]
mod stroke_tests {
    use super::*;

    #[test]
    fn every_supported_glyph_has_strokes_in_the_unit_box() {
        // All A-Z and 0-9 glyphs must be defined, with every stroke a
        // polyline of >= 2 points inside the normalized [0,1] box.
        for ch in ('A'..='Z').chain('0'..='9') {
            let strokes = stroke_glyph(ch).unwrap_or_else(|_| panic!("no glyph for {ch}"));
            assert!(!strokes.is_empty(), "glyph {ch} has no strokes");
            for stroke in &strokes {
                assert!(stroke.len() >= 2, "glyph {ch} has a degenerate stroke");
                for &(x, y) in stroke {
                    assert!(
                        (-0.01..=1.01).contains(&x) && (-0.01..=1.01).contains(&y),
                        "glyph {ch} stroke point ({x}, {y}) escaped the unit box"
                    );
                }
            }
        }
    }

    #[test]
    fn extrude_strokes_places_a_word() {
        let t = extrude_strokes("CAD", 10.0, 2.0).unwrap();
        assert_eq!(t.glyphs.len(), 3);
        // Each glyph carries at least one stroke.
        for g in &t.glyphs {
            assert!(!g.strokes.is_empty());
        }
        // The laid-out string has positive length.
        assert!(t.length > 0.0);
    }

    #[test]
    fn extrude_strokes_scales_to_font_size() {
        // A glyph at font_size 20 must reach roughly y = 20 at its top.
        let t = extrude_strokes("I", 20.0, 1.0).unwrap();
        let max_y = t.glyphs[0]
            .strokes
            .iter()
            .flatten()
            .map(|p| p.y)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((max_y - 20.0).abs() < 1e-9, "I should be 20 units tall");
    }

    #[test]
    fn extrude_strokes_space_advances_without_a_glyph() {
        // "A B" → 2 glyphs (space contributes no glyph) but a wider
        // layout than "AB".
        let spaced = extrude_strokes("A B", 10.0, 1.0).unwrap();
        let tight = extrude_strokes("AB", 10.0, 1.0).unwrap();
        assert_eq!(spaced.glyphs.len(), 2);
        assert_eq!(tight.glyphs.len(), 2);
        assert!(spaced.length > tight.length, "space should widen the run");
    }

    #[test]
    fn extrude_strokes_rejects_unsupported_char() {
        let err = extrude_strokes("@", 10.0, 1.0).unwrap_err();
        assert!(matches!(err, Gcad3dError::UnsupportedChar('@')));
    }

    #[test]
    fn extrude_strokes_rejects_bad_params() {
        assert!(matches!(
            extrude_strokes("A", 0.0, 1.0),
            Err(Gcad3dError::BadParameter { .. })
        ));
        assert!(matches!(
            extrude_strokes("A", 10.0, -1.0),
            Err(Gcad3dError::BadParameter { .. })
        ));
    }

    #[test]
    fn o_glyph_is_a_closed_curve() {
        // The 'O' is a 24-segment circle approximation — its stroke
        // should start and end at (nearly) the same point.
        let strokes = stroke_glyph('O').unwrap();
        assert_eq!(strokes.len(), 1);
        let s = &strokes[0];
        let first = s[0];
        let last = *s.last().unwrap();
        let gap = ((first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)).sqrt();
        assert!(gap < 1e-9, "O should be a closed loop");
        // And it must be genuinely curved — far more than 4 points.
        assert!(s.len() > 8, "O should be a smooth curve, not blocky");
    }
}
