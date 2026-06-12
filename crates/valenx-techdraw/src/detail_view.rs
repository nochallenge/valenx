//! Detail views — circular detail callouts with magnified close-ups.
//!
//! A detail view is the standard drawing convention for calling out a
//! small region of a parent view and showing it again, elsewhere on
//! the sheet, at a higher scale. The convention has three parts:
//!
//! 1. A **detail bubble** (circle) drawn on the parent view marking
//!    the region to expand, with a label letter (`A`, `B`, `C`…).
//! 2. The **detail view itself** — a magnified view of the geometry
//!    inside the bubble, placed at a separate location on the sheet,
//!    with its own label (e.g. `Detail A — 4:1`).
//! 3. The **magnification factor** so the reader knows the scale
//!    ratio.
//!
//! # Geometric model
//!
//! A [`DetailView`] carries:
//! - `parent_view_idx` — index of the parent [`crate::view::View`].
//! - `center` + `radius` — the detail bubble in **parent-view local
//!   mm** (same frame as the parent's `visible_edges` /
//!   `hidden_edges`).
//! - `position` — where on the sheet to place the magnified output.
//! - `magnification` — multiplier on the parent view's scale.
//! - `label` — single letter (A/B/C…) used for cross-reference.
//!
//! [`DetailView::clip_and_magnify`] takes the parent view's edges and
//! produces the *magnified, recentered* edge list ready to drop into a
//! standalone [`crate::view::View`] at the detail's `position`. Edges
//! entirely outside the bubble are dropped; edges crossing the
//! boundary are clipped at the circle.
//!
//! [`DetailView::bubble_segments`] returns a polygonal approximation
//! of the detail bubble (drawn on the parent), plus a small
//! label-leader tick.

use serde::{Deserialize, Serialize};

/// One detail view tied to a parent.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DetailView {
    /// Stable id within the drawing (sequential by [`crate::Drawing::add_detail_view`]).
    pub id: usize,
    /// Index of the parent view in the drawing's `views` vector.
    pub parent_view_idx: usize,
    /// Center of the detail bubble in **parent-view local mm**.
    pub center: [f64; 2],
    /// Radius of the detail bubble in **parent-view local mm**.
    pub radius: f64,
    /// Where on the sheet (mm) the magnified output is drawn — places
    /// the detail's own local origin here, same convention as
    /// [`crate::view::View::position`].
    pub position: [f64; 2],
    /// Magnification factor relative to the *parent's* scale. A value
    /// of `2.0` means edges are drawn twice as large as in the parent.
    pub magnification: f64,
    /// Label letter (A, B, C, …) used to identify this detail in
    /// cross-references. The auto-numbering helper
    /// [`crate::Drawing::add_detail_view`] picks the next letter.
    pub label: String,
}

impl DetailView {
    /// Construct a detail view. The `label` is whatever caller wants;
    /// [`crate::Drawing::add_detail_view`] auto-picks letters in
    /// sequence (`"A"`, `"B"`, …).
    pub fn new(
        parent_view_idx: usize,
        center: [f64; 2],
        radius: f64,
        position: [f64; 2],
        magnification: f64,
        label: impl Into<String>,
    ) -> Self {
        Self {
            id: 0,
            parent_view_idx,
            center,
            radius,
            position,
            magnification,
            label: label.into(),
        }
    }

    /// Clip the parent's edge list against this detail's bubble and
    /// magnify the surviving fragments to be placed at `position`.
    ///
    /// Output coordinates are in **detail-view local mm**, with the
    /// detail bubble's center mapped to the origin (and edges scaled
    /// by `magnification`). The renderer then translates by `position`
    /// to place the result on the sheet.
    pub fn clip_and_magnify(&self, edges: &[[(f64, f64); 2]]) -> Vec<[(f64, f64); 2]> {
        let mut out: Vec<[(f64, f64); 2]> = Vec::new();
        for seg in edges {
            for piece in clip_segment_to_circle(*seg, self.center, self.radius) {
                let a = (
                    (piece[0].0 - self.center[0]) * self.magnification,
                    (piece[0].1 - self.center[1]) * self.magnification,
                );
                let b = (
                    (piece[1].0 - self.center[0]) * self.magnification,
                    (piece[1].1 - self.center[1]) * self.magnification,
                );
                out.push([a, b]);
            }
        }
        out
    }

    /// Polygonal approximation of the detail bubble drawn on the
    /// parent, plus a leader-tick from the bubble's perimeter to a
    /// label position 1.5 mm beyond.
    ///
    /// Coordinates are in **parent-view local mm** so the caller can
    /// pass them through the same per-view transform as the parent's
    /// edges.
    pub fn bubble_segments(&self) -> Vec<[(f64, f64); 2]> {
        let mut out: Vec<[(f64, f64); 2]> = Vec::new();
        let n = 32; // 32-gon: visually round at sheet scale
        let two_pi = std::f64::consts::TAU;
        for i in 0..n {
            let a = (i as f64) * two_pi / n as f64;
            let b = ((i + 1) as f64) * two_pi / n as f64;
            let p0 = (
                self.center[0] + self.radius * a.cos(),
                self.center[1] + self.radius * a.sin(),
            );
            let p1 = (
                self.center[0] + self.radius * b.cos(),
                self.center[1] + self.radius * b.sin(),
            );
            out.push([p0, p1]);
        }
        // Leader tick at 45° outside the circle pointing up-right.
        let tick_a = (
            self.center[0] + self.radius * std::f64::consts::FRAC_1_SQRT_2,
            self.center[1] + self.radius * std::f64::consts::FRAC_1_SQRT_2,
        );
        let tick_b = (tick_a.0 + 1.5, tick_a.1 + 1.5);
        out.push([tick_a, tick_b]);
        out
    }

    /// Suggested label text for the magnified output ("Detail A —
    /// 4:1"). Embeds the magnification as a `N:1` or `1:N` ratio
    /// depending on whether it's a zoom-in (>1) or zoom-out (<1).
    pub fn detail_caption(&self) -> String {
        let scale = if self.magnification >= 1.0 {
            format!("{:.0}:1", self.magnification)
        } else {
            format!("1:{:.0}", 1.0 / self.magnification.max(1e-9))
        };
        format!("Detail {} — {scale}", self.label)
    }
}

/// Clip a 2D segment against a circle. Returns the portion(s) inside
/// the circle; either 0 (entirely outside), 1 (one fragment), or up to
/// 2 fragments (if the segment enters + exits twice — geometrically
/// possible only if it grazes a chord).
fn clip_segment_to_circle(
    seg: [(f64, f64); 2],
    center: [f64; 2],
    radius: f64,
) -> Vec<[(f64, f64); 2]> {
    let r2 = radius * radius;
    let inside = |p: (f64, f64)| -> bool {
        let dx = p.0 - center[0];
        let dy = p.1 - center[1];
        dx * dx + dy * dy <= r2
    };
    let a_in = inside(seg[0]);
    let b_in = inside(seg[1]);
    if a_in && b_in {
        return vec![seg];
    }
    // Compute intersection of segment with circle.
    let ax = seg[0].0 - center[0];
    let ay = seg[0].1 - center[1];
    let bx = seg[1].0 - center[0];
    let by = seg[1].1 - center[1];
    let dx = bx - ax;
    let dy = by - ay;
    // |a + t * d|^2 = r^2 → quadratic in t.
    let a_q = dx * dx + dy * dy;
    let b_q = 2.0 * (ax * dx + ay * dy);
    let c_q = ax * ax + ay * ay - r2;
    let disc = b_q * b_q - 4.0 * a_q * c_q;
    if disc < 0.0 || a_q.abs() < 1e-12 {
        return Vec::new();
    }
    let sq = disc.sqrt();
    let t1 = (-b_q - sq) / (2.0 * a_q);
    let t2 = (-b_q + sq) / (2.0 * a_q);
    let lerp = |t: f64| -> (f64, f64) {
        (
            seg[0].0 + t * (seg[1].0 - seg[0].0),
            seg[0].1 + t * (seg[1].1 - seg[0].1),
        )
    };
    let in_range = |t: f64| (0.0..=1.0).contains(&t);
    match (a_in, b_in) {
        (true, false) => {
            // A inside, B outside: trim B back to t2 (the exit point ≥ 0).
            let t = if in_range(t2) { t2 } else { t1 };
            if !in_range(t) {
                return Vec::new();
            }
            vec![[seg[0], lerp(t)]]
        }
        (false, true) => {
            // A outside, B inside: trim A forward to t1 (entry point).
            let t = if in_range(t1) { t1 } else { t2 };
            if !in_range(t) {
                return Vec::new();
            }
            vec![[lerp(t), seg[1]]]
        }
        (false, false) => {
            // Both outside — either misses the circle or enters and
            // exits as a chord. Both intersections must be in [0, 1].
            if in_range(t1) && in_range(t2) {
                vec![[lerp(t1), lerp(t2)]]
            } else {
                Vec::new()
            }
        }
        (true, true) => vec![seg], // handled at top
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Segment fully inside the bubble survives untouched.
    #[test]
    fn segment_inside_bubble_kept() {
        let dv = DetailView::new(0, [10.0, 10.0], 5.0, [100.0, 100.0], 2.0, "A");
        let edges = vec![[(8.0, 8.0), (12.0, 12.0)]];
        let out = dv.clip_and_magnify(&edges);
        assert_eq!(out.len(), 1);
        // Should be magnified by 2x and recentered around bubble origin.
        // Original A=(8,8), B=(12,12), center=(10,10)
        // After centering: A=(-2,-2), B=(2,2)
        // After mag x2: A=(-4,-4), B=(4,4)
        assert!((out[0][0].0 - -4.0).abs() < 1e-6);
        assert!((out[0][0].1 - -4.0).abs() < 1e-6);
        assert!((out[0][1].0 - 4.0).abs() < 1e-6);
        assert!((out[0][1].1 - 4.0).abs() < 1e-6);
    }

    /// Segment fully outside the bubble dropped.
    #[test]
    fn segment_outside_bubble_dropped() {
        let dv = DetailView::new(0, [10.0, 10.0], 5.0, [100.0, 100.0], 2.0, "A");
        let edges = vec![[(50.0, 50.0), (60.0, 60.0)]];
        let out = dv.clip_and_magnify(&edges);
        assert!(out.is_empty());
    }

    /// Segment partially inside the bubble is clipped at the boundary.
    #[test]
    fn segment_crossing_bubble_clipped_at_circle() {
        let dv = DetailView::new(0, [0.0, 0.0], 5.0, [100.0, 100.0], 1.0, "A");
        // Horizontal line at y=0 from x=-10 to x=10. Bubble center (0,0)
        // radius 5. Should clip to [-5, 5].
        let edges = vec![[(-10.0, 0.0), (10.0, 0.0)]];
        let out = dv.clip_and_magnify(&edges);
        assert_eq!(out.len(), 1);
        // At magnification=1, x ranges from -5 to 5.
        assert!((out[0][0].0 + 5.0).abs() < 1e-6, "got {}", out[0][0].0);
        assert!((out[0][1].0 - 5.0).abs() < 1e-6, "got {}", out[0][1].0);
    }

    /// Chord case: both endpoints outside, but the segment crosses
    /// through the circle.
    #[test]
    fn segment_passing_through_bubble_chord_clipped() {
        let dv = DetailView::new(0, [0.0, 0.0], 5.0, [100.0, 100.0], 1.0, "A");
        let edges = vec![[(-20.0, 0.0), (20.0, 0.0)]];
        let out = dv.clip_and_magnify(&edges);
        assert_eq!(out.len(), 1);
        // Chord from -5 to 5.
        assert!((out[0][0].0 + 5.0).abs() < 1e-6);
        assert!((out[0][1].0 - 5.0).abs() < 1e-6);
    }

    /// Magnification scales the output.
    #[test]
    fn magnification_scales_output_edges() {
        let dv = DetailView::new(0, [0.0, 0.0], 5.0, [100.0, 100.0], 4.0, "A");
        let edges = vec![[(0.0, 0.0), (1.0, 0.0)]];
        let out = dv.clip_and_magnify(&edges);
        // 1mm edge at 4x = 4mm output.
        assert_eq!(out.len(), 1);
        assert!((out[0][1].0 - 4.0).abs() < 1e-6);
    }

    /// Bubble segments form a closed polygon + leader tick.
    #[test]
    fn bubble_segments_emit_closed_polygon() {
        let dv = DetailView::new(0, [10.0, 10.0], 5.0, [100.0, 100.0], 2.0, "A");
        let segs = dv.bubble_segments();
        // 32 polygon edges + 1 leader tick.
        assert_eq!(segs.len(), 33);
        // The polygon's first and last polygon segments should share an
        // endpoint forming a closed loop.
        assert!((segs[0][0].0 - segs[31][1].0).abs() < 1e-6);
        assert!((segs[0][0].1 - segs[31][1].1).abs() < 1e-6);
    }

    /// Caption format follows "Detail A — N:1" convention.
    #[test]
    fn caption_format_for_magnification_above_one() {
        let dv = DetailView::new(0, [0.0, 0.0], 5.0, [100.0, 100.0], 4.0, "A");
        let cap = dv.detail_caption();
        assert!(cap.contains("Detail A"));
        assert!(cap.contains("4:1"));
    }

    /// Caption format for magnification below one ("1:N").
    #[test]
    fn caption_format_for_magnification_below_one() {
        let dv = DetailView::new(0, [0.0, 0.0], 5.0, [100.0, 100.0], 0.5, "B");
        let cap = dv.detail_caption();
        assert!(cap.contains("Detail B"));
        assert!(cap.contains("1:2"));
    }
}
