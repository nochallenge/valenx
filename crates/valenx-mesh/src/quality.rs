//! Element-level quality metrics.
//!
//! Two metrics ship in the first pass — **volume / area** (signed,
//! so inverted elements stand out) and **aspect ratio**. The broader
//! CFD-flavoured quantities (orthogonality, non-orthogonality angle,
//! skewness) land when the OpenFOAM prepare pipeline needs them for
//! mesh-quality-gated warnings. All computations are pure math and
//! fully test-covered; they don't depend on any mesher being
//! installed.
//!
//! Numerical conventions:
//!
//! - 3D volumes use signed tetrahedral or hexahedral volume; negative
//!   values indicate inverted connectivity.
//! - 2D areas use the cross-product formula (Tri3) or the
//!   shoelace/bilinear-patch area (Quad4). 2D areas are returned
//!   positive — unlike 3D, 2D orientation is a per-mesh convention.
//! - Aspect ratio is `edge_max / edge_min` for simplicial elements,
//!   and `longest_diagonal / shortest_diagonal` for Hex8.

use nalgebra::Vector3;

use crate::adjacency::build_face_adjacency;
use crate::{ElementBlock, ElementType, Mesh};

/// Summary statistics the UI surfaces as "mesh quality".
#[derive(Clone, Copy, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct QualityReport {
    pub element_count: u64,
    pub min_size: Option<f64>,
    pub max_size: Option<f64>,
    pub mean_size: Option<f64>,
    pub max_aspect_ratio: Option<f64>,
    /// Worst (largest) per-element equiangle skewness across the
    /// whole mesh. Range `[0, 1]`: `0` = all elements perfectly
    /// regular, approaches `1` as faces degenerate. See
    /// [`equiangle_skewness`].
    pub max_skewness: Option<f64>,
    /// Worst (smallest) cell-face orthogonality across all interior
    /// faces. Range `[0, 1]`: `1` = perfectly orthogonal mesh
    /// (cell-to-cell vector parallel to face normal everywhere),
    /// approaches `0` as cells skew away. `None` when the mesh has
    /// no interior faces. See [`min_orthogonality`].
    pub min_orthogonality: Option<f64>,
    /// Number of elements with negative signed volume (inverted).
    pub inverted_count: u64,
}

impl QualityReport {
    /// `true` when no elements are inverted, every quality field that
    /// was sampled is finite, and the minimum element size (if any)
    /// is positive — the "all-green" sanity check downstream solvers
    /// rely on.
    pub fn is_healthy(&self) -> bool {
        self.inverted_count == 0
            && self.min_size.map(|v| v > 0.0).unwrap_or(true)
            && self.max_aspect_ratio.map(|v| v.is_finite()).unwrap_or(true)
            && self.max_skewness.map(|v| v.is_finite()).unwrap_or(true)
            && self
                .min_orthogonality
                .map(|v| v.is_finite())
                .unwrap_or(true)
    }
}

impl std::fmt::Display for QualityReport {
    /// Multi-line text summary for CLI / log / debug output. Each
    /// populated field gets one line; `None` values are omitted so
    /// the output stays compact on partial reports (e.g. a mesh with
    /// no interior faces skips the orthogonality line). Format is
    /// stable enough for grep / awk pipelines but not a serialization
    /// format — use serde for that.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "elements: {}", self.element_count)?;
        if let Some(min) = self.min_size {
            writeln!(f, "min size: {min:.4e}")?;
        }
        if let Some(max) = self.max_size {
            writeln!(f, "max size: {max:.4e}")?;
        }
        if let Some(mean) = self.mean_size {
            writeln!(f, "mean size: {mean:.4e}")?;
        }
        if let Some(ar) = self.max_aspect_ratio {
            writeln!(f, "max aspect: {ar:.3}")?;
        }
        if let Some(sk) = self.max_skewness {
            writeln!(f, "max skew: {sk:.3}")?;
        }
        if let Some(orth) = self.min_orthogonality {
            writeln!(f, "min orthogonality: {orth:.3}")?;
        }
        // Always emit the inverted-count line — even at zero it's
        // useful for a quick "I checked, none inverted" confirmation.
        if self.element_count > 0 {
            writeln!(f, "inverted: {}", self.inverted_count)?;
        }
        Ok(())
    }
}

/// Run the quality pass on every element block in `mesh` and roll
/// the results into one `QualityReport`. Unsupported element types
/// contribute their element count but no size / aspect data.
pub fn report(mesh: &Mesh) -> QualityReport {
    let mut r = QualityReport::default();
    let mut size_sum = 0.0f64;
    let mut size_count: u64 = 0;
    for block in &mesh.element_blocks {
        analyse_block(block, &mesh.nodes, &mut r, &mut size_sum, &mut size_count);
    }
    if size_count > 0 {
        r.mean_size = Some(size_sum / size_count as f64);
    }
    r.min_orthogonality = min_orthogonality(mesh);
    r
}

fn analyse_block(
    block: &ElementBlock,
    nodes: &[Vector3<f64>],
    r: &mut QualityReport,
    size_sum: &mut f64,
    size_count: &mut u64,
) {
    let npe = block.element_type.nodes_per_element();
    if npe == 0 {
        return;
    }
    let element_count = block.connectivity.len() / npe;
    r.element_count = r.element_count.saturating_add(element_count as u64);

    for i in 0..element_count {
        let start = i * npe;
        let idxs = &block.connectivity[start..start + npe];
        // Bail this element if any index is out of range.
        let pts = idxs
            .iter()
            .map(|&j| nodes.get(j as usize).copied())
            .collect::<Option<Vec<_>>>();
        let Some(pts) = pts else {
            continue;
        };
        if let Some(size) = signed_size(block.element_type, &pts) {
            if size < 0.0 {
                r.inverted_count = r.inverted_count.saturating_add(1);
            }
            let abs = size.abs();
            r.min_size = Some(r.min_size.map(|m| m.min(size)).unwrap_or(size));
            r.max_size = Some(r.max_size.map(|m| m.max(abs)).unwrap_or(abs));
            *size_sum += abs;
            *size_count += 1;
        }
        if let Some(ar) = aspect_ratio(block.element_type, &pts) {
            if ar.is_finite() {
                r.max_aspect_ratio = Some(r.max_aspect_ratio.map(|m| m.max(ar)).unwrap_or(ar));
            }
        }
        if let Some(sk) = equiangle_skewness(block.element_type, &pts) {
            if sk.is_finite() {
                r.max_skewness = Some(r.max_skewness.map(|m| m.max(sk)).unwrap_or(sk));
            }
        }
    }
}

/// Default aspect-ratio histogram buckets — the upper edge of each
/// bucket. `1.0` lands "perfect" (equilateral) elements; `1e6`
/// catches the far tail of degenerate meshes that the default
/// reporter would clamp to "max" and lose visibility on.
///
/// Bucket semantics: histogram value at index `i` counts elements
/// with `aspect_ratio <= buckets[i]` AND (i == 0 OR `aspect_ratio
/// > buckets[i-1]`). The last bucket also catches +∞.
pub const DEFAULT_AR_BUCKETS: &[f64] = &[1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 1e6];

/// One element-quality histogram over [`aspect_ratio`] values, with
/// per-bucket counts + the over-cap "outliers" count for the long
/// tail. Use [`aspect_ratio_histogram`] to compute one.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct AspectRatioHistogram {
    /// Upper edge of each bucket. Always sorted ascending.
    pub buckets: Vec<f64>,
    /// `counts[i]` = elements with `aspect_ratio <= buckets[i]`
    /// AND `aspect_ratio > buckets[i-1]` (for i > 0). Length =
    /// `buckets.len()`.
    pub counts: Vec<u64>,
    /// Elements with `aspect_ratio > buckets.last()`. Always 0
    /// when the last bucket is +∞ but non-zero on user-bounded
    /// histograms.
    pub overflow: u64,
    /// Elements where aspect_ratio couldn't be computed (Hex8 +
    /// other unsupported types or degenerate elements with
    /// non-finite ratios).
    pub uncategorised: u64,
}

impl AspectRatioHistogram {
    /// Total elements counted across every bucket + overflow.
    pub fn total(&self) -> u64 {
        self.counts.iter().sum::<u64>() + self.overflow + self.uncategorised
    }

    /// Fraction of categorised elements falling in or below
    /// bucket `i`. Range `[0.0, 1.0]`. `None` when no elements
    /// were categorised.
    pub fn cumulative_fraction(&self, i: usize) -> Option<f64> {
        let denom = self.counts.iter().sum::<u64>() + self.overflow;
        if denom == 0 {
            return None;
        }
        let through_i: u64 = self.counts.iter().take(i + 1).sum();
        Some(through_i as f64 / denom as f64)
    }
}

/// Compute an aspect-ratio histogram for every element in the
/// mesh. `buckets` lists the upper edges in ascending order — pass
/// [`DEFAULT_AR_BUCKETS`] for sensible defaults, or supply custom
/// bin boundaries.
///
/// Bucket-range convention: bucket `i` catches elements with
/// `aspect_ratio <= buckets[i]` AND (i == 0 OR > buckets[i-1]).
/// Anything beyond the last bucket counts in `overflow`. Elements
/// where the ratio can't be computed (unsupported element types,
/// inf / NaN ratios from degenerate connectivity) count in
/// `uncategorised`.
pub fn aspect_ratio_histogram(mesh: &Mesh, buckets: &[f64]) -> AspectRatioHistogram {
    let mut hist = AspectRatioHistogram {
        buckets: buckets.to_vec(),
        counts: vec![0; buckets.len()],
        overflow: 0,
        uncategorised: 0,
    };
    for block in &mesh.element_blocks {
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let element_count = block.connectivity.len() / npe;
        for i in 0..element_count {
            let start = i * npe;
            let idxs = &block.connectivity[start..start + npe];
            let pts: Option<Vec<_>> = idxs
                .iter()
                .map(|&j| mesh.nodes.get(j as usize).copied())
                .collect();
            let Some(pts) = pts else {
                hist.uncategorised += 1;
                continue;
            };
            match aspect_ratio(block.element_type, &pts) {
                Some(ar) if ar.is_finite() => {
                    let mut placed = false;
                    for (bi, &edge) in buckets.iter().enumerate() {
                        if ar <= edge {
                            hist.counts[bi] += 1;
                            placed = true;
                            break;
                        }
                    }
                    if !placed {
                        hist.overflow += 1;
                    }
                }
                _ => {
                    hist.uncategorised += 1;
                }
            }
        }
    }
    hist
}

/// Signed area for 2D elements, signed volume for 3D elements,
/// length for 1D. Returns `None` for element types the metric isn't
/// defined on yet.
/// Default equiangle-skewness histogram buckets. These match the
/// quality-band convention CFD prep tools use: `<=0.25` excellent,
/// `<=0.5` good, `<=0.75` acceptable, `<=0.9` poor, `<=1.0` very
/// poor / sliver. The metric itself is bounded in `[0, 1]` so the
/// last bucket is the cap (no separate overflow concept).
pub const DEFAULT_SKEW_BUCKETS: &[f64] = &[0.25, 0.5, 0.75, 0.9, 1.0];

/// Per-bucket counts for an [`equiangle_skewness`] histogram. Mirrors
/// [`AspectRatioHistogram`] but without an `overflow` field — skewness
/// is bounded `[0, 1]`, so the last bucket is the inclusive cap.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SkewnessHistogram {
    /// Inclusive upper edge of each bucket. Always sorted ascending.
    pub buckets: Vec<f64>,
    /// `counts[i]` = elements with `skewness <= buckets[i]` AND
    /// (i == 0 OR `skewness > buckets[i-1]`).
    pub counts: Vec<u64>,
    /// Elements where skewness couldn't be computed — Line2,
    /// degenerate elements with zero-length edges, or out-of-range
    /// node indices.
    pub uncategorised: u64,
}

impl SkewnessHistogram {
    /// Total elements counted across every bucket + uncategorised.
    pub fn total(&self) -> u64 {
        self.counts.iter().sum::<u64>() + self.uncategorised
    }

    /// Fraction of categorised elements falling in or below
    /// bucket `i`. `None` when no elements were categorised.
    pub fn cumulative_fraction(&self, i: usize) -> Option<f64> {
        let denom: u64 = self.counts.iter().sum();
        if denom == 0 {
            return None;
        }
        let through_i: u64 = self.counts.iter().take(i + 1).sum();
        Some(through_i as f64 / denom as f64)
    }
}

/// Compute an equiangle-skewness histogram for every element in the
/// mesh. `buckets` lists the inclusive upper edges in ascending
/// order — pass [`DEFAULT_SKEW_BUCKETS`] for sensible defaults.
///
/// Bucket-range convention: bucket `i` catches elements with
/// `skewness <= buckets[i]` AND (i == 0 OR > buckets[i-1]). Elements
/// where the skew can't be computed (Line2, unsupported types,
/// degenerate connectivity) count in `uncategorised`.
///
/// Note: a finite skewness `> buckets.last()` is still placed into
/// the last bucket. The metric is bounded `[0, 1]`; the only way to
/// land outside is rounding noise, which we treat as "very poor"
/// rather than dropping the element.
pub fn skewness_histogram(mesh: &Mesh, buckets: &[f64]) -> SkewnessHistogram {
    let mut hist = SkewnessHistogram {
        buckets: buckets.to_vec(),
        counts: vec![0; buckets.len()],
        uncategorised: 0,
    };
    if buckets.is_empty() {
        // No bins — every element is "uncategorised" by definition.
        for block in &mesh.element_blocks {
            hist.uncategorised += block.count() as u64;
        }
        return hist;
    }
    for block in &mesh.element_blocks {
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let element_count = block.connectivity.len() / npe;
        for i in 0..element_count {
            let start = i * npe;
            let idxs = &block.connectivity[start..start + npe];
            let pts: Option<Vec<_>> = idxs
                .iter()
                .map(|&j| mesh.nodes.get(j as usize).copied())
                .collect();
            let Some(pts) = pts else {
                hist.uncategorised += 1;
                continue;
            };
            match equiangle_skewness(block.element_type, &pts) {
                Some(sk) if sk.is_finite() => {
                    let mut placed = false;
                    for (bi, &edge) in buckets.iter().enumerate() {
                        if sk <= edge {
                            hist.counts[bi] += 1;
                            placed = true;
                            break;
                        }
                    }
                    if !placed {
                        // skew > last bucket (rounding noise above
                        // 1.0) — still "very poor", land in last.
                        *hist.counts.last_mut().unwrap() += 1;
                    }
                }
                _ => {
                    hist.uncategorised += 1;
                }
            }
        }
    }
    hist
}

/// Signed size (length / area / volume) of one element of `element_type`
/// given its `pts` in source order. Returns `None` when the element type
/// is unsupported or `pts` is too short to compute the size.
pub fn signed_size(element_type: ElementType, pts: &[Vector3<f64>]) -> Option<f64> {
    match element_type {
        ElementType::Line2 => {
            let [a, b] = [pts[0], pts[1]];
            Some((b - a).norm())
        }
        ElementType::Tri3 => {
            let [a, b, c] = [pts[0], pts[1], pts[2]];
            let cross = (b - a).cross(&(c - a));
            Some(0.5 * cross.norm())
        }
        ElementType::Quad4 => {
            let [a, b, c, d] = [pts[0], pts[1], pts[2], pts[3]];
            let a1 = (b - a).cross(&(c - a)).norm() * 0.5;
            let a2 = (c - a).cross(&(d - a)).norm() * 0.5;
            Some(a1 + a2)
        }
        ElementType::Tet4 => {
            let [a, b, c, d] = [pts[0], pts[1], pts[2], pts[3]];
            let v = (b - a).dot(&((c - a).cross(&(d - a)))) / 6.0;
            Some(v)
        }
        ElementType::Hex8 => {
            // 6-tet decomposition along the body diagonal h[0]-h[6].
            // Each tet shares that diagonal as its "spine" and wraps
            // one of the six wedge faces. For an axis-aligned unit
            // cube this sums to exactly 1.
            let h: [Vector3<f64>; 8] = [
                pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7],
            ];
            let tet = |a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>, d: Vector3<f64>| {
                (b - a).dot(&((c - a).cross(&(d - a)))) / 6.0
            };
            Some(
                tet(h[0], h[1], h[2], h[6])
                    + tet(h[0], h[2], h[3], h[6])
                    + tet(h[0], h[3], h[7], h[6])
                    + tet(h[0], h[7], h[4], h[6])
                    + tet(h[0], h[4], h[5], h[6])
                    + tet(h[0], h[5], h[1], h[6]),
            )
        }
        ElementType::Pyr5 => {
            // Split the square-based pyramid into two tets along the
            // base diagonal {0,2}. Apex is node 4.
            let p: [Vector3<f64>; 5] = [pts[0], pts[1], pts[2], pts[3], pts[4]];
            let tet = |a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>, d: Vector3<f64>| {
                (b - a).dot(&((c - a).cross(&(d - a)))) / 6.0
            };
            Some(tet(p[0], p[1], p[2], p[4]) + tet(p[0], p[2], p[3], p[4]))
        }
        ElementType::Prism6 => {
            // Split the triangular prism into three tets. Standard
            // decomposition: {0,1,2,5}, {0,1,5,4}, {0,4,5,3}.
            let p: [Vector3<f64>; 6] = [pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]];
            let tet = |a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>, d: Vector3<f64>| {
                (b - a).dot(&((c - a).cross(&(d - a)))) / 6.0
            };
            Some(
                tet(p[0], p[1], p[2], p[5])
                    + tet(p[0], p[1], p[5], p[4])
                    + tet(p[0], p[4], p[5], p[3]),
            )
        }
        ElementType::Tri6 | ElementType::Tet10 | ElementType::Hex20 => {
            // Quadratic elements: reduce to linear for a first cut.
            let linear_type = match element_type {
                ElementType::Tri6 => ElementType::Tri3,
                ElementType::Tet10 => ElementType::Tet4,
                ElementType::Hex20 => ElementType::Hex8,
                _ => return None,
            };
            let n = linear_type.nodes_per_element();
            if pts.len() < n {
                return None;
            }
            signed_size(linear_type, &pts[..n])
        }
    }
}

/// Aspect ratio: largest edge / smallest edge for simplicial
/// elements, diagonal ratio for Hex8. Quadratic elements reduce to
/// their linear corner subset — sane first cut, post-processors
/// almost always render the linearised version anyway.
pub fn aspect_ratio(element_type: ElementType, pts: &[Vector3<f64>]) -> Option<f64> {
    match element_type {
        ElementType::Line2 => Some(1.0),
        ElementType::Tri3 => {
            let edges = [
                (pts[1] - pts[0]).norm(),
                (pts[2] - pts[1]).norm(),
                (pts[0] - pts[2]).norm(),
            ];
            edge_ratio(&edges)
        }
        ElementType::Quad4 => {
            let edges = [
                (pts[1] - pts[0]).norm(),
                (pts[2] - pts[1]).norm(),
                (pts[3] - pts[2]).norm(),
                (pts[0] - pts[3]).norm(),
            ];
            edge_ratio(&edges)
        }
        ElementType::Tet4 => {
            let edges = [
                (pts[1] - pts[0]).norm(),
                (pts[2] - pts[0]).norm(),
                (pts[3] - pts[0]).norm(),
                (pts[2] - pts[1]).norm(),
                (pts[3] - pts[1]).norm(),
                (pts[3] - pts[2]).norm(),
            ];
            edge_ratio(&edges)
        }
        ElementType::Pyr5 => {
            // Standard CGNS / Gmsh ordering: nodes 0..3 form the
            // quadrilateral base, node 4 is the apex. 4 base edges
            // + 4 apex edges = 8 edges.
            let edges = [
                (pts[1] - pts[0]).norm(),
                (pts[2] - pts[1]).norm(),
                (pts[3] - pts[2]).norm(),
                (pts[0] - pts[3]).norm(),
                (pts[4] - pts[0]).norm(),
                (pts[4] - pts[1]).norm(),
                (pts[4] - pts[2]).norm(),
                (pts[4] - pts[3]).norm(),
            ];
            edge_ratio(&edges)
        }
        ElementType::Prism6 => {
            // Standard ordering: nodes 0..2 = bottom triangle,
            // nodes 3..5 = top triangle (matching corner). 3 bottom
            // + 3 top + 3 vertical edges = 9 edges.
            let edges = [
                // Bottom triangle.
                (pts[1] - pts[0]).norm(),
                (pts[2] - pts[1]).norm(),
                (pts[0] - pts[2]).norm(),
                // Top triangle.
                (pts[4] - pts[3]).norm(),
                (pts[5] - pts[4]).norm(),
                (pts[3] - pts[5]).norm(),
                // Vertical edges connecting bottom -> top.
                (pts[3] - pts[0]).norm(),
                (pts[4] - pts[1]).norm(),
                (pts[5] - pts[2]).norm(),
            ];
            edge_ratio(&edges)
        }
        ElementType::Hex8 => {
            // Four principal-body diagonals.
            let d = [
                (pts[6] - pts[0]).norm(),
                (pts[7] - pts[1]).norm(),
                (pts[4] - pts[2]).norm(),
                (pts[5] - pts[3]).norm(),
            ];
            edge_ratio(&d)
        }
        ElementType::Tri6 | ElementType::Tet10 | ElementType::Hex20 => {
            // Reduce to the linear corner subset — same approach
            // signed_size uses for these types. Mid-edge nodes don't
            // change the corner-to-corner aspect ratio meaningfully
            // for a well-formed quadratic element.
            let linear_type = match element_type {
                ElementType::Tri6 => ElementType::Tri3,
                ElementType::Tet10 => ElementType::Tet4,
                ElementType::Hex20 => ElementType::Hex8,
                _ => return None,
            };
            let n = linear_type.nodes_per_element();
            if pts.len() < n {
                return None;
            }
            aspect_ratio(linear_type, &pts[..n])
        }
    }
}

/// Equiangle skewness, normalised to `[0, 1]`. `0` is a perfectly
/// regular element (all face angles equal to the ideal — 60° for tri
/// faces, 90° for quad faces); values approach `1` as faces collapse
/// toward zero-area. The per-element value is the maximum skew across
/// every face. This is the same definition Fluent and OpenFOAM's
/// `checkMesh` use.
///
/// Returns `None` for unsupported element types or for degenerate
/// geometry (coincident face nodes, zero-length edges).
pub fn equiangle_skewness(element_type: ElementType, pts: &[Vector3<f64>]) -> Option<f64> {
    let pi3 = std::f64::consts::FRAC_PI_3;
    let pi2 = std::f64::consts::FRAC_PI_2;
    match element_type {
        ElementType::Tri3 => polygon_face_skew(&pts[..3], pi3),
        ElementType::Quad4 => polygon_face_skew(&pts[..4], pi2),
        ElementType::Tet4 => {
            // 4 triangular faces. Standard exodus / CGNS face ordering.
            let faces: [(&[usize], f64); 4] = [
                (&[0, 1, 2], pi3),
                (&[0, 1, 3], pi3),
                (&[0, 2, 3], pi3),
                (&[1, 2, 3], pi3),
            ];
            max_face_skew(pts, &faces)
        }
        ElementType::Hex8 => {
            // 6 quadrilateral faces. Standard CGNS ordering: bottom
            // = 0-1-2-3, top = 4-5-6-7, vertical edges 0-4 / 1-5 /
            // 2-6 / 3-7.
            let faces: [(&[usize], f64); 6] = [
                (&[0, 1, 2, 3], pi2), // bottom (-z)
                (&[4, 5, 6, 7], pi2), // top    (+z)
                (&[0, 1, 5, 4], pi2), // front  (-y)
                (&[1, 2, 6, 5], pi2), // right  (+x)
                (&[2, 3, 7, 6], pi2), // back   (+y)
                (&[3, 0, 4, 7], pi2), // left   (-x)
            ];
            max_face_skew(pts, &faces)
        }
        ElementType::Pyr5 => {
            // 1 quad base (0-1-2-3) + 4 triangular sides sharing
            // apex node 4.
            let faces: [(&[usize], f64); 5] = [
                (&[0, 1, 2, 3], pi2), // base
                (&[0, 1, 4], pi3),
                (&[1, 2, 4], pi3),
                (&[2, 3, 4], pi3),
                (&[3, 0, 4], pi3),
            ];
            max_face_skew(pts, &faces)
        }
        ElementType::Prism6 => {
            // 2 triangular caps + 3 quadrilateral side faces. Bottom
            // tri = 0-1-2, top tri = 3-4-5 (corner-matched).
            let faces: [(&[usize], f64); 5] = [
                (&[0, 1, 2], pi3), // bottom cap
                (&[3, 4, 5], pi3), // top cap
                (&[0, 1, 4, 3], pi2),
                (&[1, 2, 5, 4], pi2),
                (&[2, 0, 3, 5], pi2),
            ];
            max_face_skew(pts, &faces)
        }
        ElementType::Tri6 | ElementType::Tet10 | ElementType::Hex20 => {
            // Quadratic: reduce to the linear corner subset, mirroring
            // signed_size / aspect_ratio.
            let linear_type = match element_type {
                ElementType::Tri6 => ElementType::Tri3,
                ElementType::Tet10 => ElementType::Tet4,
                ElementType::Hex20 => ElementType::Hex8,
                _ => return None,
            };
            let n = linear_type.nodes_per_element();
            if pts.len() < n {
                return None;
            }
            equiangle_skewness(linear_type, &pts[..n])
        }
        ElementType::Line2 => None,
    }
}

/// Walk every face index-set, compute its polygon skew using the
/// face's `theta_eq`, and return the maximum. Returns `None` if any
/// face is degenerate or out-of-bounds — one bad face flunks the
/// whole element.
fn max_face_skew(pts: &[Vector3<f64>], faces: &[(&[usize], f64)]) -> Option<f64> {
    let mut max_s = 0.0_f64;
    for (face, theta_eq) in faces {
        let face_pts: Option<Vec<Vector3<f64>>> =
            face.iter().map(|&i| pts.get(i).copied()).collect();
        let face_pts = face_pts?;
        let s = polygon_face_skew(&face_pts, *theta_eq)?;
        if s > max_s {
            max_s = s;
        }
    }
    Some(max_s)
}

/// Equiangle skew of a single planar polygon face. `theta_eq` is the
/// ideal interior angle for the face: 60° (π/3) for triangles, 90°
/// (π/2) for quads. Returns `None` for zero-length edges.
fn polygon_face_skew(pts: &[Vector3<f64>], theta_eq: f64) -> Option<f64> {
    let n = pts.len();
    if n < 3 {
        return None;
    }
    let pi = std::f64::consts::PI;
    let mut min_a = pi;
    let mut max_a = 0.0_f64;
    for i in 0..n {
        let prev = (i + n - 1) % n;
        let next = (i + 1) % n;
        let a = pts[prev] - pts[i];
        let b = pts[next] - pts[i];
        let na = a.norm();
        let nb = b.norm();
        if !na.is_finite() || !nb.is_finite() || na == 0.0 || nb == 0.0 {
            return None;
        }
        let cos = (a.dot(&b) / (na * nb)).clamp(-1.0, 1.0);
        let theta = cos.acos();
        if theta < min_a {
            min_a = theta;
        }
        if theta > max_a {
            max_a = theta;
        }
    }
    let high = (max_a - theta_eq) / (pi - theta_eq);
    let low = (theta_eq - min_a) / theta_eq;
    Some(high.max(low).max(0.0))
}

/// Minimum cell-face **orthogonality** across every interior face of
/// the mesh, range `[0, 1]`. `1` is perfectly orthogonal — the
/// vector connecting the two cell centroids is parallel to the
/// shared face's normal, the FV best case. Values approach `0` as
/// the cell-to-cell vector tilts toward the face plane (skewed /
/// non-orthogonal mesh, the FV worst case).
///
/// Formula (per interior face):
/// ```text
/// orthogonality = |dot(d_unit, n_unit)|
/// where d = right_centroid - left_centroid
///       n = face_normal
/// ```
/// We use the absolute value so the sign of the face normal (which
/// depends on which element claims it) doesn't matter.
///
/// Returns `None` when the mesh has no interior faces (single
/// element, or only boundary faces) — there's nothing to compare.
pub fn min_orthogonality(mesh: &Mesh) -> Option<f64> {
    let adj = build_face_adjacency(mesh);
    if adj.interior_face_count() == 0 {
        return None;
    }
    let mut min = f64::INFINITY;
    for face in adj.interior_faces() {
        let Some(left_c) = global_element_centroid(mesh, face.left.global_element) else {
            continue;
        };
        let Some(right_c) = global_element_centroid(mesh, face.right.global_element) else {
            continue;
        };
        let Some(n) = face_plane_normal(mesh, &face.sorted_nodes) else {
            continue;
        };
        let d = right_c - left_c;
        let dn = d.norm();
        let nn = n.norm();
        if dn == 0.0 || nn == 0.0 {
            continue;
        }
        let cos = (d.dot(&n) / (dn * nn)).abs().min(1.0);
        if cos < min {
            min = cos;
        }
    }
    if min.is_finite() {
        Some(min)
    } else {
        None
    }
}

/// Centroid of the element with the given global index — arithmetic
/// mean of its node coordinates. Walks `mesh.element_blocks` in order
/// to find which block owns this element. Returns `None` for an
/// out-of-range index or for an element whose connectivity points
/// at an out-of-range node.
fn global_element_centroid(mesh: &Mesh, global_element: usize) -> Option<Vector3<f64>> {
    let mut cursor: usize = 0;
    for block in &mesh.element_blocks {
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let count = block.connectivity.len() / npe;
        if global_element < cursor + count {
            let local = global_element - cursor;
            let start = local * npe;
            let conn = &block.connectivity[start..start + npe];
            let mut sum = Vector3::zeros();
            for &idx in conn {
                let p = *mesh.nodes.get(idx as usize)?;
                sum += p;
            }
            return Some(sum / (npe as f64));
        }
        cursor += count;
    }
    None
}

/// Face normal (unnormalised) for a face given its sorted node index
/// list. For 3-node tri faces: cross product of two edges. For 4-node
/// quad faces: cross product of the two diagonals (Newell-equivalent
/// for planar quads). Returns `None` for degenerate faces or
/// out-of-range indices. Direction is convention-dependent — callers
/// that need orthogonality should `.abs()` the dot product.
fn face_plane_normal(mesh: &Mesh, sorted_nodes: &[u32]) -> Option<Vector3<f64>> {
    let pts: Option<Vec<Vector3<f64>>> = sorted_nodes
        .iter()
        .map(|&i| mesh.nodes.get(i as usize).copied())
        .collect();
    let pts = pts?;
    match pts.len() {
        3 => {
            let n = (pts[1] - pts[0]).cross(&(pts[2] - pts[0]));
            if n.norm() == 0.0 {
                None
            } else {
                Some(n)
            }
        }
        4 => {
            // Diagonals 0-2 and 1-3.
            let n = (pts[2] - pts[0]).cross(&(pts[3] - pts[1]));
            if n.norm() == 0.0 {
                None
            } else {
                Some(n)
            }
        }
        _ => None,
    }
}

fn edge_ratio(edges: &[f64]) -> Option<f64> {
    let mut min = f64::INFINITY;
    let mut max: f64 = 0.0;
    for &e in edges {
        if !e.is_finite() {
            return None;
        }
        if e < min {
            min = e;
        }
        if e > max {
            max = e;
        }
    }
    if min <= 0.0 {
        return Some(f64::INFINITY);
    }
    Some(max / min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tet_with_positive_volume() {
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let v = signed_size(ElementType::Tet4, &pts).unwrap();
        assert!((v - 1.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn inverted_tet_has_negative_volume() {
        // Swap v2 and v3 to reverse orientation.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let v = signed_size(ElementType::Tet4, &pts).unwrap();
        assert!(v < 0.0);
    }

    #[test]
    fn right_triangle_area_is_half() {
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let a = signed_size(ElementType::Tri3, &pts).unwrap();
        assert!((a - 0.5).abs() < 1e-12);
    }

    #[test]
    fn cube_hex_volume() {
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let v = signed_size(ElementType::Hex8, &pts).unwrap();
        assert!((v - 1.0).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn equilateral_triangle_aspect_is_one() {
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, (0.75f64).sqrt(), 0.0),
        ];
        let ar = aspect_ratio(ElementType::Tri3, &pts).unwrap();
        assert!((ar - 1.0).abs() < 1e-9, "got {ar}");
    }

    #[test]
    fn skewed_triangle_aspect_is_large() {
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let ar = aspect_ratio(ElementType::Tri3, &pts).unwrap();
        assert!(ar > 50.0, "got {ar}");
    }

    #[test]
    fn report_on_unit_box_matches_volume() {
        // Two tets that together fill the unit cube [0,1]^3 have
        // total volume 1. One tet here for simplicity.
        let mut m = Mesh::new("unit-tet");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);
        m.recompute_stats();
        let r = report(&m);
        assert_eq!(r.element_count, 1);
        assert!(r.is_healthy());
        assert!((r.max_size.unwrap() - 1.0 / 6.0).abs() < 1e-12);
        assert!(r.inverted_count == 0);
    }

    #[test]
    fn report_flags_inverted_elements() {
        let mut m = Mesh::new("inv");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);
        let r = report(&m);
        assert_eq!(r.inverted_count, 1);
        assert!(!r.is_healthy());
    }

    // -----------------------------------------------------------------
    // Aspect-ratio histogram
    // -----------------------------------------------------------------

    #[test]
    fn default_ar_buckets_are_sensibly_ordered() {
        // Sanity: ascending + last bucket is much larger than typical
        // skewed-element ratios.
        for w in DEFAULT_AR_BUCKETS.windows(2) {
            assert!(w[0] < w[1], "buckets must be ascending: {w:?}");
        }
        assert!(*DEFAULT_AR_BUCKETS.last().unwrap() >= 1e6);
    }

    #[test]
    fn aspect_ratio_histogram_groups_equilateral_into_first_bucket() {
        // One equilateral triangle — aspect ratio 1.0 lands in the
        // smallest bucket (<= 1.5).
        let mut m = Mesh::new("eq");
        let h = (3.0_f64).sqrt() / 2.0;
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let hist = aspect_ratio_histogram(&m, DEFAULT_AR_BUCKETS);
        assert_eq!(hist.counts[0], 1);
        for c in &hist.counts[1..] {
            assert_eq!(*c, 0);
        }
        assert_eq!(hist.overflow, 0);
        assert_eq!(hist.uncategorised, 0);
    }

    #[test]
    fn aspect_ratio_histogram_total_matches_element_count() {
        let mut m = Mesh::new("two");
        // Two tris: equilateral + a long skinny one.
        let h = (3.0_f64).sqrt() / 2.0;
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            // skinny: long edge of 100 + 1 + 1 -> aspect ~ 100
            Vector3::new(10.0, 10.0, 0.0),
            Vector3::new(110.0, 10.0, 0.0),
            Vector3::new(60.0, 10.5, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2, 3, 4, 5];
        m.element_blocks.push(block);
        let hist = aspect_ratio_histogram(&m, DEFAULT_AR_BUCKETS);
        assert_eq!(hist.total(), 2);
    }

    #[test]
    fn aspect_ratio_histogram_overflow_catches_very_skewed() {
        // Sliver triangle: long edge + two equal short legs that
        // collapse to ~ε. Max edge ~100, min edge ~ε -> aspect ratio
        // tens of thousands, way past any bucket.
        let mut m = Mesh::new("over");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.001, 0.0, 0.0),
            Vector3::new(100.0, 0.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        // Buckets cap at 5.0; AR ~100/0.001 = 100000 lands in overflow.
        let hist = aspect_ratio_histogram(&m, &[1.5, 2.0, 5.0]);
        assert_eq!(hist.overflow, 1);
        for c in &hist.counts {
            assert_eq!(*c, 0);
        }
    }

    #[test]
    fn aspect_ratio_histogram_handles_hex8_via_diagonals() {
        // Unit cube: all four body diagonals equal sqrt(3); ratio is
        // exactly 1.0 -> first bucket. Regression-anchors the Hex8
        // diagonal-ratio path so a future "missing implementation"
        // can't silently bucket cubes as uncategorised.
        let mut m = Mesh::new("hex");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = (0..8).collect();
        m.element_blocks.push(block);
        let hist = aspect_ratio_histogram(&m, DEFAULT_AR_BUCKETS);
        assert_eq!(hist.total(), 1);
        assert_eq!(hist.uncategorised, 0, "Hex8 must categorise cleanly");
        assert_eq!(
            hist.counts[0], 1,
            "unit cube AR=1 should land in first bucket"
        );
    }

    #[test]
    fn aspect_ratio_histogram_cumulative_fraction_zeroes_to_one() {
        // 3 equilateral tris -> all in bucket 0; cumulative
        // fraction at bucket 0 should be 1.0.
        let h = (3.0_f64).sqrt() / 2.0;
        let mut m = Mesh::new("three");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(2.5, h, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(4.5, h, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = (0..9).collect();
        m.element_blocks.push(block);
        let hist = aspect_ratio_histogram(&m, DEFAULT_AR_BUCKETS);
        let cum = hist.cumulative_fraction(0).unwrap();
        assert!((cum - 1.0).abs() < 1e-12, "got cumulative {cum}");
        // Last bucket also reads 1.0 (cumulative reaches everything).
        let last_idx = DEFAULT_AR_BUCKETS.len() - 1;
        let cum_last = hist.cumulative_fraction(last_idx).unwrap();
        assert!((cum_last - 1.0).abs() < 1e-12);
    }

    #[test]
    fn aspect_ratio_histogram_cumulative_fraction_none_for_empty_mesh() {
        let m = Mesh::new("empty");
        let hist = aspect_ratio_histogram(&m, DEFAULT_AR_BUCKETS);
        assert!(hist.cumulative_fraction(0).is_none());
    }

    // -----------------------------------------------------------------
    // Volumetric / quadratic element coverage
    // -----------------------------------------------------------------

    #[test]
    fn pyramid_signed_volume_matches_one_third_base_height() {
        // Unit-square base [0,1]^2 at z=0, apex at (0.5, 0.5, 1.0).
        // V = (1/3) * base_area * height = (1/3) * 1 * 1 = 1/3.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5, 1.0),
        ];
        let v = signed_size(ElementType::Pyr5, &pts).unwrap();
        assert!((v - 1.0 / 3.0).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn pyramid_aspect_ratio_for_regular_pyramid() {
        // Square base of edge 1, apex centred at height 1. Apex
        // edges all equal sqrt((0.5)^2 + (0.5)^2 + 1^2) = sqrt(1.5).
        // edge_max = sqrt(1.5) ~ 1.2247, edge_min = 1.0.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5, 1.0),
        ];
        let ar = aspect_ratio(ElementType::Pyr5, &pts).unwrap();
        let expected = (1.5_f64).sqrt() / 1.0;
        assert!((ar - expected).abs() < 1e-9, "got {ar}, want {expected}");
    }

    #[test]
    fn prism_signed_volume_for_right_triangular_prism() {
        // Right triangular prism: base triangle area 0.5, extruded
        // along +z by 1. Volume = 0.5.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let v = signed_size(ElementType::Prism6, &pts).unwrap();
        assert!((v - 0.5).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn prism_aspect_ratio_for_right_prism_uses_largest_edge() {
        // Same prism as above. Edges: bottom (1, sqrt(2), 1), top
        // (1, sqrt(2), 1), vertical (1, 1, 1). max = sqrt(2), min = 1.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let ar = aspect_ratio(ElementType::Prism6, &pts).unwrap();
        let expected = (2.0_f64).sqrt();
        assert!((ar - expected).abs() < 1e-9, "got {ar}");
    }

    #[test]
    fn quadratic_aspect_ratio_reduces_to_linear() {
        // Tet10: corners 0..3 form a unit-edge regular tet, mid-edge
        // nodes 4..9 sit anywhere (we don't read them). aspect_ratio
        // should match the linear Tet4 of the corners.
        let corners = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, (3.0_f64).sqrt() / 2.0, 0.0),
            Vector3::new(0.5, (3.0_f64).sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()),
        ];
        // Append 6 mid-edge nodes (positions don't matter to the
        // reduce-to-linear path).
        let mut quad_pts = corners.to_vec();
        for _ in 0..6 {
            quad_pts.push(Vector3::new(99.0, 99.0, 99.0));
        }
        let lin = aspect_ratio(ElementType::Tet4, &corners).unwrap();
        let quad = aspect_ratio(ElementType::Tet10, &quad_pts).unwrap();
        assert!((lin - quad).abs() < 1e-12, "linear {lin} vs quad {quad}");
        // Regular tet -> aspect ratio exactly 1.0.
        assert!((quad - 1.0).abs() < 1e-9, "got {quad}");
    }

    #[test]
    fn quadratic_signed_size_reduces_to_linear() {
        // Tet10 unit-tet corners + dummy mid-edge nodes. signed_size
        // should match Tet4 corners exactly.
        let corners = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut quad_pts = corners.to_vec();
        for _ in 0..6 {
            quad_pts.push(Vector3::new(99.0, 99.0, 99.0));
        }
        let lin = signed_size(ElementType::Tet4, &corners).unwrap();
        let quad = signed_size(ElementType::Tet10, &quad_pts).unwrap();
        assert!((lin - quad).abs() < 1e-12);
    }

    #[test]
    fn aspect_ratio_histogram_counts_pyr5_and_prism6() {
        // Mixed-element mesh with one Pyr5 and one Prism6 — both
        // should now categorise cleanly (no uncategorised count).
        let mut m = Mesh::new("mixed");
        m.nodes = vec![
            // Pyr5: nodes 0..4
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5, 1.0),
            // Prism6: nodes 5..10
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),
            Vector3::new(2.0, 0.0, 1.0),
            Vector3::new(3.0, 0.0, 1.0),
            Vector3::new(2.0, 1.0, 1.0),
        ];
        let mut pyr = ElementBlock::new(ElementType::Pyr5);
        pyr.connectivity = vec![0, 1, 2, 3, 4];
        m.element_blocks.push(pyr);
        let mut prism = ElementBlock::new(ElementType::Prism6);
        prism.connectivity = vec![5, 6, 7, 8, 9, 10];
        m.element_blocks.push(prism);
        let hist = aspect_ratio_histogram(&m, DEFAULT_AR_BUCKETS);
        assert_eq!(hist.total(), 2);
        assert_eq!(hist.uncategorised, 0, "all element types now supported");
    }

    #[test]
    fn report_on_pyramid_mesh_uses_volume() {
        // Quality report on a single Pyr5 element should pick up
        // the 1/3 volume + healthy state.
        let mut m = Mesh::new("pyr");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Pyr5);
        block.connectivity = vec![0, 1, 2, 3, 4];
        m.element_blocks.push(block);
        let r = report(&m);
        assert_eq!(r.element_count, 1);
        assert!(r.is_healthy());
        assert!((r.max_size.unwrap() - 1.0 / 3.0).abs() < 1e-9);
        assert_eq!(r.inverted_count, 0);
    }

    #[test]
    fn equilateral_tri_has_zero_equiangle_skewness() {
        // All three angles are 60° — the equilateral case is the
        // textbook zero-skewness baseline.
        let h = (3.0_f64).sqrt() / 2.0;
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
        ];
        let s = equiangle_skewness(ElementType::Tri3, &pts).unwrap();
        assert!(s.abs() < 1e-12, "expected ~0, got {s}");
    }

    #[test]
    fn right_isoceles_tri_skewness_is_one_quarter() {
        // Angles 90 / 45 / 45. Both branches give 0.25:
        //   high = (90-60)/(180-60) = 30/120 = 0.25
        //   low  = (60-45)/60       = 15/60  = 0.25
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let s = equiangle_skewness(ElementType::Tri3, &pts).unwrap();
        assert!((s - 0.25).abs() < 1e-12, "expected 0.25, got {s}");
    }

    #[test]
    fn degenerate_tri_skewness_is_none() {
        // Two coincident vertices — zero-length edge means we can't
        // form an angle. Should not panic, should not return Inf.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        assert!(equiangle_skewness(ElementType::Tri3, &pts).is_none());
    }

    #[test]
    fn unit_square_quad_has_zero_skewness() {
        // All four interior angles are 90°.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let s = equiangle_skewness(ElementType::Quad4, &pts).unwrap();
        assert!(s.abs() < 1e-12, "expected ~0, got {s}");
    }

    #[test]
    fn rhombus_quad_skewness_is_one_third() {
        // 60° / 120° / 60° / 120° rhombus.
        //   high = (120-90)/(180-90) = 30/90 = 1/3
        //   low  = (90-60)/90        = 30/90 = 1/3
        let h = (3.0_f64).sqrt() / 2.0;
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.5, h, 0.0),
            Vector3::new(0.5, h, 0.0),
        ];
        let s = equiangle_skewness(ElementType::Quad4, &pts).unwrap();
        assert!((s - 1.0 / 3.0).abs() < 1e-12, "expected 1/3, got {s}");
    }

    #[test]
    fn regular_tet_has_zero_skewness() {
        // Regular tetrahedron inscribed in a cube: all four faces are
        // equilateral triangles, so every face angle is exactly 60°.
        let pts = [
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(1.0, -1.0, -1.0),
            Vector3::new(-1.0, 1.0, -1.0),
            Vector3::new(-1.0, -1.0, 1.0),
        ];
        let s = equiangle_skewness(ElementType::Tet4, &pts).unwrap();
        assert!(s.abs() < 1e-12, "expected ~0, got {s}");
    }

    #[test]
    fn unit_cube_hex_has_zero_skewness() {
        // All six faces are unit squares.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let s = equiangle_skewness(ElementType::Hex8, &pts).unwrap();
        assert!(s.abs() < 1e-12, "expected ~0, got {s}");
    }

    #[test]
    fn sheared_hex_skewness_matches_parallelogram_face() {
        // Top face sheared +x by 0.5: front and back faces become
        // parallelograms with corner angles atan(2) and pi - atan(2).
        // Skew per side face = 1 - 2*atan(2)/pi (~0.2952). The 4
        // unsheared faces (top/bottom/left/right) stay at 0 skew, so
        // the per-element max equals the side-face skew.
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.0, 1.0),
            Vector3::new(1.5, 0.0, 1.0),
            Vector3::new(1.5, 1.0, 1.0),
            Vector3::new(0.5, 1.0, 1.0),
        ];
        let s = equiangle_skewness(ElementType::Hex8, &pts).unwrap();
        let expected = 1.0 - 2.0 * (2.0_f64).atan() / std::f64::consts::PI;
        assert!((s - expected).abs() < 1e-12, "expected {expected}, got {s}");
    }

    #[test]
    fn regular_pyramid_has_zero_skewness() {
        // Square base (side 1) with apex above the centre at height
        // 1/sqrt(2) — slant edge length 1, so all four side faces
        // are equilateral triangles. Base + sides all 0 skew.
        let h = 1.0 / (2.0_f64).sqrt();
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5, h),
        ];
        let s = equiangle_skewness(ElementType::Pyr5, &pts).unwrap();
        assert!(s.abs() < 1e-12, "expected ~0, got {s}");
    }

    #[test]
    fn right_equilateral_prism_has_zero_skewness() {
        // Equilateral triangular cross-section, vertical extrusion.
        // Both caps are equilateral, all 3 side faces are unit
        // rectangles — every face contributes 0 skew.
        let h = (3.0_f64).sqrt() / 2.0;
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(0.5, h, 1.0),
        ];
        let s = equiangle_skewness(ElementType::Prism6, &pts).unwrap();
        assert!(s.abs() < 1e-12, "expected ~0, got {s}");
    }

    #[test]
    fn quadratic_skewness_reduces_to_linear() {
        // Tri6 with mid-edge nodes anywhere — corner subset is the
        // right-isoceles triangle, so skew should match Tri3 result.
        let pts = [
            // Corner nodes 0..2.
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            // Mid-edge nodes 3..5 — bogus, skewness must ignore them.
            Vector3::new(0.5, 0.0, 0.0),
            Vector3::new(0.5, 0.5, 0.0),
            Vector3::new(0.0, 0.5, 0.0),
        ];
        let s = equiangle_skewness(ElementType::Tri6, &pts).unwrap();
        assert!((s - 0.25).abs() < 1e-12, "expected 0.25, got {s}");
    }

    #[test]
    fn line2_skewness_is_none() {
        // Lines have no faces — skewness undefined.
        let pts = [Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        assert!(equiangle_skewness(ElementType::Line2, &pts).is_none());
    }

    #[test]
    fn report_max_skewness_picks_worst_element() {
        // Two-tri mesh: one equilateral (0 skew) + one right-isoceles
        // (0.25 skew). The report's max_skewness should be 0.25.
        let h = (3.0_f64).sqrt() / 2.0;
        let mut m = Mesh::new("two-tris");
        m.nodes = vec![
            // Equilateral 0-1-2.
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            // Right-isoceles 3-4-5.
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2, 3, 4, 5];
        m.element_blocks.push(block);
        let r = report(&m);
        assert_eq!(r.element_count, 2);
        let s = r.max_skewness.expect("skewness should be populated");
        assert!((s - 0.25).abs() < 1e-12, "expected 0.25, got {s}");
        assert!(r.is_healthy());
    }

    #[test]
    fn skewness_histogram_empty_mesh_has_zero_counts() {
        let m = Mesh::new("empty");
        let h = skewness_histogram(&m, DEFAULT_SKEW_BUCKETS);
        assert_eq!(h.total(), 0);
        assert_eq!(h.uncategorised, 0);
        assert!(h.counts.iter().all(|&c| c == 0));
        assert_eq!(h.buckets.len(), DEFAULT_SKEW_BUCKETS.len());
    }

    #[test]
    fn skewness_histogram_places_each_element_in_one_bucket() {
        // Three Tri3s with known skewness:
        //   tri 0-1-2  : equilateral       -> skew 0    -> bucket 0 (<=0.25)
        //   tri 3-4-5  : right-isoceles    -> skew 0.25 -> bucket 0 (<=0.25)
        //   tri 6-7-8  : 30-60-90          -> skew 0.5  -> bucket 1 (<=0.5)
        let h = (3.0_f64).sqrt() / 2.0;
        let mut m = Mesh::new("mixed-skew");
        m.nodes = vec![
            // Equilateral.
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            // Right-isoceles.
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),
            // 30-60-90: legs sqrt(3) and 1, hypotenuse 2. Place at
            // (4,0)-(4+sqrt(3),0)-(4,1).
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(4.0 + (3.0_f64).sqrt(), 0.0, 0.0),
            Vector3::new(4.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7, 8];
        m.element_blocks.push(block);
        let hist = skewness_histogram(&m, DEFAULT_SKEW_BUCKETS);
        assert_eq!(hist.total(), 3);
        assert_eq!(hist.uncategorised, 0);
        // Bucket 0 (<=0.25): equilateral + right-isoceles.
        assert_eq!(hist.counts[0], 2);
        // Bucket 1 (<=0.5): 30-60-90.
        assert_eq!(hist.counts[1], 1);
        // Higher buckets empty.
        assert!(hist.counts[2..].iter().all(|&c| c == 0));
    }

    #[test]
    fn skewness_histogram_line2_counts_as_uncategorised() {
        let mut m = Mesh::new("line");
        m.nodes = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        let mut block = ElementBlock::new(ElementType::Line2);
        block.connectivity = vec![0, 1];
        m.element_blocks.push(block);
        let hist = skewness_histogram(&m, DEFAULT_SKEW_BUCKETS);
        assert_eq!(hist.uncategorised, 1);
        assert_eq!(hist.total(), 1);
        assert!(hist.counts.iter().all(|&c| c == 0));
    }

    #[test]
    fn skewness_histogram_degenerate_tri_counts_as_uncategorised() {
        // Two coincident vertices -> equiangle_skewness returns None.
        let mut m = Mesh::new("degen");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let hist = skewness_histogram(&m, DEFAULT_SKEW_BUCKETS);
        assert_eq!(hist.uncategorised, 1);
        assert_eq!(hist.total(), 1);
    }

    #[test]
    fn skewness_histogram_cumulative_fraction_sums_through_index() {
        // Two equilateral tris (bucket 0) + two 30-60-90 tris
        // (bucket 1). Cumulative fraction at i=0 = 0.5; at i=1 = 1.0.
        let h = (3.0_f64).sqrt() / 2.0;
        let mut m = Mesh::new("cumfrac");
        m.nodes = vec![
            // 2x equilateral.
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(2.5, h, 0.0),
            // 2x 30-60-90.
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(4.0 + (3.0_f64).sqrt(), 0.0, 0.0),
            Vector3::new(4.0, 1.0, 0.0),
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(6.0 + (3.0_f64).sqrt(), 0.0, 0.0),
            Vector3::new(6.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = (0..12).collect();
        m.element_blocks.push(block);
        let hist = skewness_histogram(&m, DEFAULT_SKEW_BUCKETS);
        assert_eq!(hist.cumulative_fraction(0), Some(0.5));
        assert_eq!(hist.cumulative_fraction(1), Some(1.0));
        assert_eq!(
            hist.cumulative_fraction(DEFAULT_SKEW_BUCKETS.len() - 1),
            Some(1.0)
        );
    }

    #[test]
    fn skewness_histogram_empty_buckets_uncategorises_everything() {
        // No bins -> can't place anything; tracks element count.
        let h = (3.0_f64).sqrt() / 2.0;
        let mut m = Mesh::new("nobuckets");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let hist = skewness_histogram(&m, &[]);
        assert_eq!(hist.uncategorised, 1);
        assert_eq!(hist.total(), 1);
        assert_eq!(hist.cumulative_fraction(0), None);
    }

    #[test]
    fn min_orthogonality_returns_none_when_no_interior_faces() {
        // Single tet -> no interior face -> nothing to compare.
        let mut m = Mesh::new("single-tet");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);
        assert!(min_orthogonality(&m).is_none());
    }

    #[test]
    fn min_orthogonality_top_sheared_hex_matches_analytic() {
        // Hex A is the unit cube, z in [0,1]. Hex B sits on top
        // with its BOTTOM face shared with A's top (nodes 4..7) but
        // its TOP face shifted +x by 0.5 (nodes 8..11). The shared
        // face normal is +z. Centroids: A = (0.5, 0.5, 0.5),
        // B = (0.75, 0.5, 1.5) (top half pulls +x). Cell-to-cell
        // vector d = (0.25, 0, 1), |d| = sqrt(17/16). Face normal
        // (0,0,1) so d.n = 1. Orthogonality = |1/|d|| = 4/sqrt(17).
        let mut m = Mesh::new("top-sheared-hexes");
        m.nodes = vec![
            // Hex A bottom (z=0)
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            // Shared layer (z=1)
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            // Hex B top (z=2, sheared +x by 0.5)
            Vector3::new(0.5, 0.0, 2.0),
            Vector3::new(1.5, 0.0, 2.0),
            Vector3::new(1.5, 1.0, 2.0),
            Vector3::new(0.5, 1.0, 2.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11];
        m.element_blocks.push(block);
        let o = min_orthogonality(&m).expect("interior face exists");
        let expected = 4.0 / (17.0_f64).sqrt();
        assert!((o - expected).abs() < 1e-12, "expected {expected}, got {o}");
    }

    #[test]
    fn quality_report_display_includes_each_populated_field() {
        // A right-isoceles Tri3 produces: 1 element, AR=sqrt(2),
        // skew=0.25, no orthogonality (no interior face), 0
        // inverted. Display should render all populated values
        // and skip unpopulated ones (no "min orthogonality" line).
        let mut m = Mesh::new("right-isoceles");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let r = report(&m);
        let s = format!("{r}");
        assert!(s.contains("elements: 1"), "missing element count: {s}");
        assert!(s.contains("max aspect"), "missing aspect line: {s}");
        assert!(s.contains("max skew"), "missing skew line: {s}");
        // No interior faces -> orthogonality should be skipped.
        assert!(
            !s.contains("min orthogonality"),
            "orthogonality should be omitted when None: {s}"
        );
        // Healthy mesh -> the inverted-count line says "0 inverted".
        assert!(s.contains("inverted: 0"), "missing inverted line: {s}");
    }

    #[test]
    fn quality_report_display_includes_orthogonality_when_present() {
        // Two stacked unit cubes give orthogonality = 1.0; should
        // surface in the Display output.
        let mut m = Mesh::new("stacked");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(0.0, 0.0, 2.0),
            Vector3::new(1.0, 0.0, 2.0),
            Vector3::new(1.0, 1.0, 2.0),
            Vector3::new(0.0, 1.0, 2.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11];
        m.element_blocks.push(block);
        let r = report(&m);
        let s = format!("{r}");
        assert!(s.contains("min orthogonality"), "missing orth line: {s}");
        assert!(s.contains("1.000"), "expected 1.000 orth value: {s}");
    }

    #[test]
    fn quality_report_display_empty_mesh_is_minimal() {
        // Empty mesh -> "elements: 0" plus essentially nothing else.
        // Useful for `println!("{report}")` on a freshly-loaded mesh
        // before any quality pass has run.
        let r = QualityReport::default();
        let s = format!("{r}");
        assert!(s.contains("elements: 0"));
        assert!(!s.contains("max aspect"));
        assert!(!s.contains("max skew"));
        assert!(!s.contains("min orthogonality"));
    }

    #[test]
    fn min_orthogonality_two_axis_aligned_hexes_is_one() {
        // Two unit cubes stacked along z. The shared face is z=1
        // with normal (0,0,1); centroids are (0.5, 0.5, 0.5) and
        // (0.5, 0.5, 1.5), so the cell-to-cell vector is also
        // (0,0,1). Orthogonality = |dot(n, d)/(|n||d|)| = 1.0.
        let mut m = Mesh::new("two-stacked-hexes");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(0.0, 0.0, 2.0),
            Vector3::new(1.0, 0.0, 2.0),
            Vector3::new(1.0, 1.0, 2.0),
            Vector3::new(0.0, 1.0, 2.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11];
        m.element_blocks.push(block);
        let o = min_orthogonality(&m).expect("should have an interior face");
        assert!((o - 1.0).abs() < 1e-12, "expected 1.0, got {o}");
    }
}
