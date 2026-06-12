//! Broken views — long views with one or more break regions.
//!
//! A broken view is the standard drawing convention for objects too
//! long to fit on the sheet at the desired scale: a strip of the view
//! is "removed" (the omitted region is collapsed) and replaced with a
//! zigzag break symbol. The two surviving fragments are then shifted
//! together so the final drawing is shorter.
//!
//! Linear breaks only (horizontal or vertical strips); radial breaks
//! are explicitly out of scope for the v1 commercial-depth pass.
//!
//! # Geometric model
//!
//! A [`BreakRegion`] picks an axis ([`BreakAxis`]) and an interval
//! `[lo, hi]` in **view-local mm** (the same frame as
//! [`crate::view::View::visible_edges`]). [`apply_breaks`] is the
//! workhorse: given a list of input edges (in view-local mm) and a
//! list of break regions, it returns:
//!
//! 1. The surviving edge fragments after each edge is clipped against
//!    every break region, with the post-break fragments shifted by
//!    `-Σ (hi - lo)` so the break collapses cleanly.
//! 2. The zigzag break symbol segments at each break's center.
//!
//! Edges that *cross* a break are split at the boundaries; the
//! fragment inside the break is dropped; the fragment past the break
//! is shifted by the break's collapsed width on its axis. Edges
//! entirely outside every break are kept (just shifted by the total
//! collapsed width if they're past the break).
//!
//! Dimensions whose extent crosses a break carry the omitted span
//! plus an asterisk per the standard convention — exposed through
//! [`break_aware_dimension_label`].

use serde::{Deserialize, Serialize};

/// Which axis a [`BreakRegion`] cuts perpendicular to.
///
/// `Horizontal` removes a horizontal strip (the break lines run
/// horizontally; surviving fragments are above + below the strip).
/// `Vertical` removes a vertical strip (break lines run vertically;
/// surviving fragments are left + right of the strip).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BreakAxis {
    /// Break perpendicular to the x-axis (removes a *vertical* strip
    /// of the view, leaving the left and right portions intact).
    Vertical,
    /// Break perpendicular to the y-axis (removes a *horizontal* strip,
    /// leaving the top and bottom portions intact).
    Horizontal,
}

/// Style of break symbol drawn in the gap. v1 supports zigzag only.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BreakStyle {
    /// Standard zigzag (mechanical-drawing convention).
    Zigzag,
}

/// One break region in a view.
///
/// `axis` chooses the cut direction; `lo` and `hi` are the
/// removed-strip's range in view-local mm. `lo < hi` is enforced by
/// [`Self::new`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BreakRegion {
    /// Which axis to cut perpendicular to.
    pub axis: BreakAxis,
    /// Lower coordinate of the removed strip (mm).
    pub lo: f64,
    /// Upper coordinate of the removed strip (mm).
    pub hi: f64,
    /// Symbol drawn in the gap.
    pub style: BreakStyle,
}

impl BreakRegion {
    /// Construct a region. Swaps `lo` / `hi` if given out of order.
    pub fn new(axis: BreakAxis, lo: f64, hi: f64) -> Self {
        let (lo, hi) = if lo <= hi { (lo, hi) } else { (hi, lo) };
        Self {
            axis,
            lo,
            hi,
            style: BreakStyle::Zigzag,
        }
    }

    /// Width of this break's removed strip on its cut axis.
    pub fn span(&self) -> f64 {
        self.hi - self.lo
    }
}

/// Result of [`apply_breaks`]: the shifted-and-clipped edges plus the
/// break symbol segments to overlay.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct BrokenEdges {
    /// Edge fragments after every break has been applied. Same frame
    /// (view-local mm) as the input edges, but shifted by the
    /// cumulative collapsed widths so fragments on the "high" side of
    /// each break sit flush against the "low" side.
    pub edges: Vec<[(f64, f64); 2]>,
    /// Break symbol segments (zigzag) at every break, drawn at the
    /// break's center on the cut axis and spanning the *full
    /// orthogonal extent* of the view (caller passes the extent —
    /// derived from the input edges' bounding box).
    pub break_symbols: Vec<[(f64, f64); 2]>,
}

/// Clip + shift the input edges through the given break regions.
///
/// Regions don't need to be sorted; they're sorted internally by `lo`
/// per-axis. Overlapping regions on the same axis are merged
/// (the merge takes the union). Edges fully inside any region are
/// dropped; edges that cross a boundary are split.
///
/// The returned [`BrokenEdges::break_symbols`] are zigzag polylines
/// drawn at each break's *collapsed* center, sized by the orthogonal
/// bbox of the input.
pub fn apply_breaks(edges: &[[(f64, f64); 2]], regions: &[BreakRegion]) -> BrokenEdges {
    if regions.is_empty() {
        return BrokenEdges {
            edges: edges.to_vec(),
            break_symbols: Vec::new(),
        };
    }

    // Compute orthogonal extents of the input edges per axis — used
    // for drawing the zigzag symbol later.
    let (xmin, xmax, ymin, ymax) = bbox(edges);

    // Process horizontal breaks first, then vertical (or vice versa —
    // order is deterministic). Each call to `apply_axis` collapses one
    // axis's worth of breaks and returns the new edges.
    let h_regions: Vec<&BreakRegion> = regions
        .iter()
        .filter(|r| r.axis == BreakAxis::Horizontal)
        .collect();
    let v_regions: Vec<&BreakRegion> = regions
        .iter()
        .filter(|r| r.axis == BreakAxis::Vertical)
        .collect();

    let mut current: Vec<[(f64, f64); 2]> = edges.to_vec();
    // Sort + merge per axis.
    let h_merged = merge_regions(&h_regions);
    let v_merged = merge_regions(&v_regions);
    current = apply_axis(&current, &h_merged, BreakAxis::Horizontal);
    current = apply_axis(&current, &v_merged, BreakAxis::Vertical);

    // Build zigzag break symbols. Each symbol is at the (collapsed)
    // mid-point of its axis range, with the orthogonal span set to the
    // *original* extent so the zigzag visually spans the view's full
    // width (or height).
    let mut break_symbols: Vec<[(f64, f64); 2]> = Vec::new();
    // Horizontal breaks: the cut is at y=lo (after collapse). Zigzag
    // runs horizontally across [xmin, xmax].
    let mut h_shift_acc = 0.0;
    for r in &h_merged {
        let y_at = r.lo - h_shift_acc;
        push_zigzag(
            &mut break_symbols,
            (xmin, y_at),
            (xmax, y_at),
            6,   // 6 teeth across the width — looks good at typical scales
            1.5, // amplitude (mm)
            true,
        );
        h_shift_acc += r.span();
    }
    let mut v_shift_acc = 0.0;
    for r in &v_merged {
        let x_at = r.lo - v_shift_acc;
        push_zigzag(
            &mut break_symbols,
            (x_at, ymin),
            (x_at, ymax),
            6,
            1.5,
            false,
        );
        v_shift_acc += r.span();
    }

    BrokenEdges {
        edges: current,
        break_symbols,
    }
}

/// Apply every break on a single axis in sorted order.
///
/// `regions` are already merged + sorted by `lo`. For each region,
/// every edge is processed: parts inside the region drop, parts past
/// the region get shifted left/down by the cumulative collapsed span.
fn apply_axis(
    edges: &[[(f64, f64); 2]],
    regions: &[BreakRegion],
    axis: BreakAxis,
) -> Vec<[(f64, f64); 2]> {
    if regions.is_empty() {
        return edges.to_vec();
    }
    let mut out: Vec<[(f64, f64); 2]> = edges.to_vec();
    let mut shift_total = 0.0;
    for r in regions {
        let mut next: Vec<[(f64, f64); 2]> = Vec::with_capacity(out.len());
        for seg in &out {
            // Apply the running shift to the segment first (the segment
            // was already shifted by previous breaks — we add this
            // break's shift only to the parts on its "high" side).
            for piece in clip_segment_through_break(*seg, r, axis) {
                next.push(piece);
            }
        }
        shift_total += r.span();
        let _ = shift_total; // kept for clarity / future debug use
        out = next;
    }
    out
}

/// Split a single segment by one break and shift the high-side
/// fragment back. Returns 0, 1, or 2 surviving fragments.
fn clip_segment_through_break(
    seg: [(f64, f64); 2],
    r: &BreakRegion,
    axis: BreakAxis,
) -> Vec<[(f64, f64); 2]> {
    // Pick the axis coordinate accessor + mutator.
    let coord = |p: (f64, f64)| -> f64 {
        match axis {
            BreakAxis::Vertical => p.0,
            BreakAxis::Horizontal => p.1,
        }
    };
    let span = r.span();
    let shift = |p: (f64, f64)| -> (f64, f64) {
        match axis {
            BreakAxis::Vertical => (p.0 - span, p.1),
            BreakAxis::Horizontal => (p.0, p.1 - span),
        }
    };
    let ca = coord(seg[0]);
    let cb = coord(seg[1]);
    let (lo, hi) = (r.lo, r.hi);
    // Three regions: below lo, inside [lo, hi], above hi.
    // Categorise each endpoint.
    let cat = |c: f64| -> i32 {
        if c < lo {
            -1
        } else if c > hi {
            1
        } else {
            0
        }
    };
    let ka = cat(ca);
    let kb = cat(cb);
    // Both below: keep as-is, no shift.
    if ka == -1 && kb == -1 {
        return vec![seg];
    }
    // Both inside: drop.
    if ka == 0 && kb == 0 {
        return vec![];
    }
    // Both above: shift both endpoints down.
    if ka == 1 && kb == 1 {
        return vec![[shift(seg[0]), shift(seg[1])]];
    }
    // Mixed cases: parametric interpolation at the break boundaries.
    // Solve t for axis(P(t)) = boundary along the segment direction.
    // P(t) = a + t * (b - a), 0 ≤ t ≤ 1.
    let t_at = |boundary: f64| -> Option<f64> {
        let denom = cb - ca;
        if denom.abs() < 1e-12 {
            return None;
        }
        Some(((boundary - ca) / denom).clamp(0.0, 1.0))
    };
    let interp = |t: f64| -> (f64, f64) {
        (
            seg[0].0 + t * (seg[1].0 - seg[0].0),
            seg[0].1 + t * (seg[1].1 - seg[0].1),
        )
    };
    let mut out: Vec<[(f64, f64); 2]> = Vec::new();
    // Below-side fragment (one endpoint below lo).
    if ka == -1 || kb == -1 {
        let t = t_at(lo).unwrap_or(0.5);
        let below_end = if ka == -1 { seg[0] } else { seg[1] };
        let p_lo = interp(t);
        out.push([below_end, p_lo]);
    }
    // Above-side fragment (one endpoint above hi), shifted down.
    if ka == 1 || kb == 1 {
        let t = t_at(hi).unwrap_or(0.5);
        let above_end = if ka == 1 { seg[0] } else { seg[1] };
        let p_hi = interp(t);
        out.push([shift(p_hi), shift(above_end)]);
    }
    out
}

/// Merge overlapping break regions on the same axis. Output is sorted
/// by `lo`. Idempotent for already-disjoint inputs.
fn merge_regions(regions: &[&BreakRegion]) -> Vec<BreakRegion> {
    if regions.is_empty() {
        return Vec::new();
    }
    let mut v: Vec<BreakRegion> = regions.iter().map(|r| **r).collect();
    v.sort_by(|a, b| a.lo.partial_cmp(&b.lo).unwrap_or(std::cmp::Ordering::Equal));
    let mut merged: Vec<BreakRegion> = Vec::with_capacity(v.len());
    for r in v {
        if let Some(last) = merged.last_mut() {
            if r.lo <= last.hi {
                last.hi = last.hi.max(r.hi);
                continue;
            }
        }
        merged.push(r);
    }
    merged
}

fn bbox(edges: &[[(f64, f64); 2]]) -> (f64, f64, f64, f64) {
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for seg in edges {
        for p in seg {
            if p.0 < xmin {
                xmin = p.0;
            }
            if p.0 > xmax {
                xmax = p.0;
            }
            if p.1 < ymin {
                ymin = p.1;
            }
            if p.1 > ymax {
                ymax = p.1;
            }
        }
    }
    if !xmin.is_finite() {
        (0.0, 0.0, 0.0, 0.0)
    } else {
        (xmin, xmax, ymin, ymax)
    }
}

/// Push a zigzag polyline between two points. `n_teeth` controls how
/// many V-segments, `amp` is the perpendicular amplitude in mm.
/// `horizontal=true` means the line runs horizontally and the
/// amplitude is in y.
fn push_zigzag(
    out: &mut Vec<[(f64, f64); 2]>,
    a: (f64, f64),
    b: (f64, f64),
    n_teeth: usize,
    amp: f64,
    horizontal: bool,
) {
    if n_teeth == 0 {
        out.push([a, b]);
        return;
    }
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let n = (n_teeth * 2).max(2);
    let mut prev = a;
    for i in 1..=n {
        let t = i as f64 / n as f64;
        let mut p = (a.0 + dx * t, a.1 + dy * t);
        if i < n {
            // Alternate +amp / -amp perpendicular to the line.
            let sign = if i % 2 == 1 { 1.0 } else { -1.0 };
            if horizontal {
                p.1 += amp * sign;
            } else {
                p.0 += amp * sign;
            }
        }
        out.push([prev, p]);
        prev = p;
    }
}

/// Render the standard "X*" label for a dimension whose extent crosses
/// at least one break. The asterisk warns the reader that the printed
/// value is *the omitted span included* (since the dimension still
/// measures the true distance pre-break).
///
/// For dimensions that don't cross any break, returns the input value
/// formatted without the asterisk.
pub fn break_aware_dimension_label(
    from: [f64; 2],
    to: [f64; 2],
    true_value: f64,
    regions: &[BreakRegion],
) -> String {
    let crosses = regions.iter().any(|r| dim_crosses_region(from, to, r));
    if crosses {
        format!("{true_value:.2}*")
    } else {
        format!("{true_value:.2}")
    }
}

/// Does a linear dimension between `from` and `to` overlap a break
/// region on its measurement axis?
fn dim_crosses_region(from: [f64; 2], to: [f64; 2], r: &BreakRegion) -> bool {
    let (a, b) = match r.axis {
        BreakAxis::Vertical => (from[0], to[0]),
        BreakAxis::Horizontal => (from[1], to[1]),
    };
    let (lo_dim, hi_dim) = if a <= b { (a, b) } else { (b, a) };
    // Overlap between [lo_dim, hi_dim] and [r.lo, r.hi].
    lo_dim < r.hi && hi_dim > r.lo
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Edge entirely below the break survives untouched.
    #[test]
    fn edge_below_break_untouched() {
        let r = BreakRegion::new(BreakAxis::Vertical, 10.0, 20.0);
        let seg = [(0.0, 0.0), (5.0, 0.0)];
        let out = clip_segment_through_break(seg, &r, BreakAxis::Vertical);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], seg);
    }

    /// Edge fully inside the break dropped.
    #[test]
    fn edge_inside_break_dropped() {
        let r = BreakRegion::new(BreakAxis::Vertical, 10.0, 20.0);
        let seg = [(12.0, 0.0), (18.0, 0.0)];
        let out = clip_segment_through_break(seg, &r, BreakAxis::Vertical);
        assert!(out.is_empty());
    }

    /// Edge fully above the break shifted down by span.
    #[test]
    fn edge_above_break_shifted_down_by_span() {
        let r = BreakRegion::new(BreakAxis::Vertical, 10.0, 20.0);
        let seg = [(30.0, 5.0), (40.0, 5.0)];
        let out = clip_segment_through_break(seg, &r, BreakAxis::Vertical);
        assert_eq!(out.len(), 1);
        // Shift by 10 mm (the break span).
        assert!((out[0][0].0 - 20.0).abs() < 1e-9);
        assert!((out[0][1].0 - 30.0).abs() < 1e-9);
    }

    /// Edge crossing the break splits into two fragments — below part
    /// keeps its left end + clipped to break lo, above part is shifted.
    #[test]
    fn edge_crossing_break_splits_into_two() {
        let r = BreakRegion::new(BreakAxis::Vertical, 10.0, 20.0);
        let seg = [(0.0, 0.0), (30.0, 0.0)];
        let out = clip_segment_through_break(seg, &r, BreakAxis::Vertical);
        assert_eq!(out.len(), 2, "should split into below + above");
        // First fragment: 0 → 10.
        assert!((out[0][0].0 - 0.0).abs() < 1e-9);
        assert!((out[0][1].0 - 10.0).abs() < 1e-9);
        // Second fragment: was 20 → 30, shifted to 10 → 20.
        assert!((out[1][0].0 - 10.0).abs() < 1e-9);
        assert!((out[1][1].0 - 20.0).abs() < 1e-9);
    }

    /// apply_breaks integration: a 100 mm-long edge with one 50 mm
    /// vertical break becomes a 50 mm composite edge.
    #[test]
    fn apply_breaks_shrinks_long_edge() {
        let edges = vec![[(0.0, 0.0), (100.0, 0.0)]];
        let regions = vec![BreakRegion::new(BreakAxis::Vertical, 30.0, 80.0)];
        let out = apply_breaks(&edges, &regions);
        // Two fragments, total length = 30 + (100 - 80) = 50.
        assert_eq!(out.edges.len(), 2);
        let span: f64 = out.edges.iter().map(|s| (s[1].0 - s[0].0).abs()).sum();
        assert!((span - 50.0).abs() < 1e-6, "post-break length got {span}");
        // Should emit zigzag break symbols.
        assert!(!out.break_symbols.is_empty());
    }

    /// Two vertical breaks compose: total collapsed width = sum of
    /// spans.
    #[test]
    fn apply_two_vertical_breaks_composes() {
        let edges = vec![[(0.0, 0.0), (100.0, 0.0)]];
        let regions = vec![
            BreakRegion::new(BreakAxis::Vertical, 20.0, 30.0),
            BreakRegion::new(BreakAxis::Vertical, 60.0, 70.0),
        ];
        let out = apply_breaks(&edges, &regions);
        // Three fragments, total length = 20 + 30 + 30 = 80.
        let span: f64 = out.edges.iter().map(|s| (s[1].0 - s[0].0).abs()).sum();
        assert!((span - 80.0).abs() < 1e-6, "total post-break got {span}");
        assert_eq!(out.break_symbols.len() / 12, 2, "two zigzag breaks");
    }

    /// merge_regions collapses overlapping ranges.
    #[test]
    fn merge_regions_unions_overlapping_ranges() {
        let r1 = BreakRegion::new(BreakAxis::Vertical, 10.0, 30.0);
        let r2 = BreakRegion::new(BreakAxis::Vertical, 20.0, 40.0);
        let m = merge_regions(&[&r1, &r2]);
        assert_eq!(m.len(), 1);
        assert!((m[0].lo - 10.0).abs() < 1e-9);
        assert!((m[0].hi - 40.0).abs() < 1e-9);
    }

    /// merge_regions sorts by `lo`.
    #[test]
    fn merge_regions_sorts_output() {
        let r1 = BreakRegion::new(BreakAxis::Vertical, 50.0, 60.0);
        let r2 = BreakRegion::new(BreakAxis::Vertical, 10.0, 20.0);
        let m = merge_regions(&[&r1, &r2]);
        assert_eq!(m.len(), 2);
        assert!(m[0].lo < m[1].lo);
    }

    /// Horizontal break shrinks vertical edges.
    #[test]
    fn horizontal_break_shrinks_vertical_edge() {
        let edges = vec![[(5.0, 0.0), (5.0, 100.0)]];
        let regions = vec![BreakRegion::new(BreakAxis::Horizontal, 30.0, 80.0)];
        let out = apply_breaks(&edges, &regions);
        let span: f64 = out.edges.iter().map(|s| (s[1].1 - s[0].1).abs()).sum();
        assert!((span - 50.0).abs() < 1e-6);
    }

    /// new() swaps lo/hi if given out of order.
    #[test]
    fn break_region_normalises_lo_hi_order() {
        let r = BreakRegion::new(BreakAxis::Vertical, 20.0, 10.0);
        assert_eq!(r.lo, 10.0);
        assert_eq!(r.hi, 20.0);
    }

    /// break_aware_dimension_label adds an asterisk when the dim
    /// crosses a break.
    #[test]
    fn dimension_crossing_break_gets_asterisk() {
        let r = BreakRegion::new(BreakAxis::Vertical, 30.0, 80.0);
        let lbl = break_aware_dimension_label([10.0, 0.0], [100.0, 0.0], 100.0, &[r]);
        assert!(lbl.ends_with('*'), "expected asterisk, got `{lbl}`");
    }

    /// Dimension fully below the break — no asterisk.
    #[test]
    fn dimension_not_crossing_break_no_asterisk() {
        let r = BreakRegion::new(BreakAxis::Vertical, 30.0, 80.0);
        let lbl = break_aware_dimension_label([0.0, 0.0], [20.0, 0.0], 20.0, &[r]);
        assert!(!lbl.contains('*'), "no asterisk expected, got `{lbl}`");
    }

    /// apply_breaks with empty regions is a passthrough.
    #[test]
    fn empty_regions_passthrough() {
        let edges = vec![[(0.0, 0.0), (1.0, 1.0)]];
        let out = apply_breaks(&edges, &[]);
        assert_eq!(out.edges, edges);
        assert!(out.break_symbols.is_empty());
    }

    /// Integration: take a real `View` from a 100 mm-long box, apply a
    /// 50 mm break in the middle, verify the resulting view edges fit
    /// in a ~50 mm footprint.
    #[test]
    fn integration_break_long_box_to_half_length() {
        use crate::view::{View, ViewKind};
        use valenx_cad::primitives::box_solid;
        let long_box = box_solid(100.0, 10.0, 10.0).unwrap();
        let mut v = View::new(ViewKind::Front, 1.0, [0.0, 0.0]);
        v.generate(&long_box).unwrap();
        let (min_b, max_b) = v.bbox().unwrap();
        let orig_w = max_b[0] - min_b[0];
        assert!(
            (orig_w - 100.0).abs() < 1.0,
            "front-view width pre-break {orig_w}"
        );
        // Apply one vertical break of width 50 in the middle.
        let regions = vec![BreakRegion::new(BreakAxis::Vertical, 25.0, 75.0)];
        let out = apply_breaks(&v.visible_edges, &regions);
        let (xmin, xmax, _, _) = bbox(&out.edges);
        let post_w = xmax - xmin;
        assert!(
            post_w < orig_w * 0.6,
            "post-break width should be roughly half, got {post_w}"
        );
        assert!(
            !out.break_symbols.is_empty(),
            "should emit zigzag break symbol"
        );
    }

    /// Dimension crossing the break + true_value: label format is
    /// "100.00*" — the standard convention.
    #[test]
    fn integration_dimension_label_format_with_break() {
        let regions = vec![BreakRegion::new(BreakAxis::Vertical, 25.0, 75.0)];
        let label = break_aware_dimension_label([0.0, -5.0], [100.0, -5.0], 100.0, &regions);
        assert_eq!(label, "100.00*");
    }
}
