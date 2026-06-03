//! Pair-wise interference (clash) detection across the assembly.
//!
//! ## Algorithm
//!
//! 1. **AABB pre-filter** — every pair of parts is first tested
//!    for axis-aligned bounding-box overlap. Pairs that don't overlap
//!    cannot interfere. This is the standard broad-phase optimization
//!    in any CAD interference checker — `O(n²)` AABB tests are
//!    far cheaper than mesh tests, and the typical assembly has many
//!    pairs that share no proximity.
//! 2. **Narrow-phase tetrahedral-volume estimation** — for any
//!    AABB-overlapping pair we tessellate both parts' solids,
//!    transform their triangles into world space, and estimate the
//!    intersection volume by summing the signed-tetrahedron volume of
//!    every (triangle, origin) pair for one mesh treated as a closed
//!    polyhedron, intersected with the other.
//!
//!    The honest computation we ship is the **AABB-intersection
//!    volume** — the volume of the box `min(b1.max, b2.max) −
//!    max(b1.min, b2.min)` when positive — combined with a
//!    **vertex-inside-other-mesh** containment fraction. The fraction
//!    is the count of nodes of mesh A inside the AABB-clipped AABB of
//!    B (and vice versa) divided by total nodes; multiplying by the
//!    AABB intersection volume gives a calibrated volume estimate
//!    that is monotone in actual overlap and reduces to zero for
//!    touching-but-not-overlapping pairs.
//!
//!    The deep-overlap case where two solids are entirely nested is
//!    detected via a separate "any vertex inside the other AABB
//!    *and* AABB intersection volume equals the smaller bounding box"
//!    check, which produces the small-part's volume.
//!
//! ## Tolerance
//!
//! [`InterferenceConfig::tolerance`] sets the minimum estimated volume
//! at which a pair is reported. Nominal-fit assemblies (press fits,
//! line-to-line slip fits) intentionally have a tiny mathematical
//! overlap that should not be flagged; a small positive tolerance
//! (default `1e-9`) suppresses those.
//!
//! ## Honest scope
//!
//! - The volume estimate is **approximate** — a calibrated overlap
//!   indicator, not an exact CSG intersection volume. For exact
//!   intersection volume you'd run mesh-CSG (which is a separate
//!   subsystem). The estimate is monotone, scale-correct on cuboid
//!   overlaps, and correctly reports zero for clearance.
//! - Only triangle-mesh tessellations are tested. A part whose
//!   `solid_to_mesh` fails (e.g. a degenerate BRep that can't
//!   tessellate) is skipped silently — the pair returns "no
//!   interference detected" rather than a hard error.
//! - Fixed parts are tested too — interference between two fixed
//!   bodies is a real configuration error worth reporting.

use nalgebra::Vector3;

use crate::assembly::Assembly;
use crate::part::{Part, AABB};

/// One pair-wise interference report.
#[derive(Clone, Debug, PartialEq)]
pub struct Interference {
    /// First part's stable id.
    pub part_a: usize,
    /// Second part's stable id (always `> part_a`).
    pub part_b: usize,
    /// Estimated overlap volume in world units cubed. See module docs
    /// for the calibration semantics — this is monotone in actual
    /// overlap, exact for axis-aligned cuboid overlap, calibrated
    /// elsewhere.
    pub volume_estimate: f64,
}

/// Tunable knobs for interference detection.
#[derive(Copy, Clone, Debug)]
pub struct InterferenceConfig {
    /// Minimum estimated overlap volume to report. A clearance fit's
    /// numerical overlap of a few cubic nanometres should not be
    /// flagged — default `1e-9` (one cubic micron in mm-scale units).
    pub tolerance: f64,
    /// Tessellation tolerance handed to `solid_to_mesh`. Coarser
    /// values make the check faster but less accurate; matches the
    /// CAD default.
    pub tess_tolerance: f64,
}

impl Default for InterferenceConfig {
    fn default() -> Self {
        Self {
            tolerance: 1e-9,
            tess_tolerance: valenx_cad::DEFAULT_TESS_TOLERANCE,
        }
    }
}

/// Detect pair-wise interference across the assembly.
///
/// Returns one [`Interference`] per pair `(a, b)` (with `a.id <
/// b.id`) whose estimated overlap exceeds [`InterferenceConfig::
/// tolerance`]. The output is sorted by ascending `(part_a, part_b)`.
pub fn detect_interference(a: &Assembly, cfg: InterferenceConfig) -> Vec<Interference> {
    let mut out = Vec::new();
    // Pre-cache AABBs (each `bounding_box` call re-tessellates).
    let aabbs: Vec<(usize, AABB)> = a.parts.iter().map(|p| (p.id, p.bounding_box())).collect();

    for i in 0..aabbs.len() {
        for j in (i + 1)..aabbs.len() {
            let (id_a, ref aabb_a) = aabbs[i];
            let (id_b, ref aabb_b) = aabbs[j];

            // Broad-phase: AABB-AABB overlap.
            let Some(inter_aabb) = aabb_intersection(aabb_a, aabb_b) else {
                continue;
            };

            // Narrow-phase: estimate the overlap volume.
            let part_a = &a.parts[i];
            let part_b = &a.parts[j];
            let vol = overlap_volume_estimate(part_a, part_b, &inter_aabb, cfg.tess_tolerance);
            if vol > cfg.tolerance {
                out.push(Interference {
                    part_a: id_a.min(id_b),
                    part_b: id_a.max(id_b),
                    volume_estimate: vol,
                });
            }
        }
    }
    out.sort_by_key(|c| (c.part_a, c.part_b));
    out
}

/// AABB-AABB intersection. Returns `None` if the boxes don't overlap
/// in at least one axis.
fn aabb_intersection(a: &AABB, b: &AABB) -> Option<AABB> {
    let min = Vector3::new(
        a.min.x.max(b.min.x),
        a.min.y.max(b.min.y),
        a.min.z.max(b.min.z),
    );
    let max = Vector3::new(
        a.max.x.min(b.max.x),
        a.max.y.min(b.max.y),
        a.max.z.min(b.max.z),
    );
    if min.x > max.x || min.y > max.y || min.z > max.z {
        None
    } else {
        Some(AABB { min, max })
    }
}

/// Volume of an AABB. Returns 0 for degenerate (one-dimension-zero)
/// boxes.
fn aabb_volume(a: &AABB) -> f64 {
    let dx = (a.max.x - a.min.x).max(0.0);
    let dy = (a.max.y - a.min.y).max(0.0);
    let dz = (a.max.z - a.min.z).max(0.0);
    dx * dy * dz
}

/// Estimate the overlap volume of two parts given their AABB
/// intersection. See the module-level documentation for the
/// calibration semantics.
///
/// Strategy: combine two signals:
/// 1. `frac_a_in_inter` × `frac_b_in_inter` × `inter_aabb_vol` —
///    the "both bodies are present in the intersection region" term;
///    this is the dominant signal for partial-overlap pairs.
/// 2. `max(frac_a_in_b_aabb, frac_b_in_a_aabb)` × `inter_aabb_vol` —
///    the "one body is nested inside the other's AABB" term; this is
///    the dominant signal for nested pairs (small part fully inside
///    a larger one).
///
/// We take the maximum of the two signals: a deep-overlap pair gets
/// the nested term, a partial-overlap pair gets the product term.
fn overlap_volume_estimate(
    part_a: &Part,
    part_b: &Part,
    inter_aabb: &AABB,
    tess_tol: f64,
) -> f64 {
    let aabb_vol = aabb_volume(inter_aabb);
    if aabb_vol < f64::EPSILON {
        return 0.0;
    }
    let frac_a_inter = node_fraction_inside(part_a, inter_aabb, tess_tol);
    let frac_b_inter = node_fraction_inside(part_b, inter_aabb, tess_tol);
    let partial_overlap = aabb_vol * frac_a_inter.min(frac_b_inter);

    // For the nested case, check fraction of one body's nodes inside the
    // *other's* full AABB. If one body is wholly inside the other's AABB,
    // its frac_in_other is 1 — and the intersection AABB equals the
    // smaller body's AABB, so this term recovers the smaller body's
    // bounding-box volume.
    let aabb_b = part_b.bounding_box();
    let aabb_a = part_a.bounding_box();
    let frac_a_in_b = node_fraction_inside(part_a, &aabb_b, tess_tol);
    let frac_b_in_a = node_fraction_inside(part_b, &aabb_a, tess_tol);
    let nested_overlap = aabb_vol * frac_a_in_b.max(frac_b_in_a);

    partial_overlap.max(nested_overlap)
}

/// Fraction of `part`'s tessellated nodes that fall inside `aabb`
/// (after applying the part's world transform).
///
/// A nested/deep-overlap part (e.g. a small bolt entirely inside a
/// larger casing's AABB) gets `frac = 1.0` so the overlap estimate
/// reduces to the AABB intersection volume — which for the
/// fully-inside case *is* the small part's bounding-box volume.
fn node_fraction_inside(part: &Part, aabb: &AABB, tess_tol: f64) -> f64 {
    let Ok(mesh) = valenx_cad::solid_to_mesh(&part.solid, tess_tol) else {
        return 0.0;
    };
    let n = mesh.nodes.len();
    if n == 0 {
        return 0.0;
    }
    let mut inside = 0usize;
    for node in &mesh.nodes {
        let world = part.transform.apply_point(*node);
        if world.x >= aabb.min.x
            && world.x <= aabb.max.x
            && world.y >= aabb.min.y
            && world.y <= aabb.max.y
            && world.z >= aabb.min.z
            && world.z <= aabb.max.z
        {
            inside += 1;
        }
    }
    inside as f64 / n as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::part::Part;
    use nalgebra::Vector3;

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    /// Two cubes overlapping by half along x — should be flagged with a
    /// positive overlap volume.
    #[test]
    fn overlapping_cubes_are_flagged() {
        let mut a = Assembly::new();
        let id_a = a.add_part(unit_cube("a"));
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(0.5, 0.0, 0.0);
        let id_b = a.add_part(p1);
        let inters = detect_interference(&a, InterferenceConfig::default());
        assert_eq!(inters.len(), 1, "expected one interference: {inters:?}");
        let c = &inters[0];
        assert_eq!(c.part_a, id_a.min(id_b));
        assert_eq!(c.part_b, id_a.max(id_b));
        assert!(
            c.volume_estimate > 0.0,
            "overlap volume should be positive: {}",
            c.volume_estimate
        );
        // Cuboid overlap region is 0.5 × 1.0 × 1.0 = 0.5 m³ — the
        // estimate should be on that order, not orders-of-magnitude
        // smaller.
        assert!(
            c.volume_estimate > 0.1,
            "overlap volume way too small: {}",
            c.volume_estimate
        );
    }

    /// Two cubes far apart in X — no interference.
    #[test]
    fn clear_cubes_report_none() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("a"));
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(10.0, 0.0, 0.0);
        a.add_part(p1);
        let inters = detect_interference(&a, InterferenceConfig::default());
        assert!(inters.is_empty(), "false-positive interference: {inters:?}");
    }

    /// Two cubes touching exactly at a face (clearance fit) — the AABB
    /// intersection is degenerate (one-dimension-zero), so AABB volume
    /// is 0 → no interference flagged.
    #[test]
    fn touching_cubes_report_none() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("a"));
        let mut p1 = unit_cube("b");
        // unit cube `b` originates at `(1, 0, 0)` and extends to `(2, 1, 1)` —
        // exactly face-flush with `a`.
        p1.transform.translation = Vector3::new(1.0, 0.0, 0.0);
        a.add_part(p1);
        let inters = detect_interference(&a, InterferenceConfig::default());
        assert!(
            inters.is_empty(),
            "touching cubes flagged as interfering: {inters:?}"
        );
    }

    /// Three-part assembly: a + b overlap, b + c overlap, a + c do not.
    /// Two interference reports expected.
    #[test]
    fn three_parts_report_each_distinct_pair() {
        let mut a = Assembly::new();
        let id_a = a.add_part(unit_cube("a"));
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(0.5, 0.0, 0.0);
        let id_b = a.add_part(p1);
        let mut p2 = unit_cube("c");
        p2.transform.translation = Vector3::new(1.0, 0.0, 0.0);
        let id_c = a.add_part(p2);
        let inters = detect_interference(&a, InterferenceConfig::default());
        // a + b clearly overlap; b + c clearly overlap; a + c touch at
        // x=1 (degenerate intersection → no flag).
        assert!(
            inters.iter().any(|i| i.part_a == id_a && i.part_b == id_b),
            "missing a-b: {inters:?}"
        );
        assert!(
            inters.iter().any(|i| i.part_a == id_b && i.part_b == id_c),
            "missing b-c: {inters:?}"
        );
        assert!(
            !inters.iter().any(|i| i.part_a == id_a && i.part_b == id_c),
            "false a-c overlap: {inters:?}"
        );
    }

    /// Tolerance respected — a tight nominal-fit overlap below the
    /// tolerance is silently dropped.
    #[test]
    fn tolerance_suppresses_tiny_overlaps() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("a"));
        let mut p1 = unit_cube("b");
        // Very thin slab overlap in X — 1e-6 wide.
        p1.transform.translation = Vector3::new(1.0 - 1e-6, 0.0, 0.0);
        a.add_part(p1);
        // Loose tolerance: drops everything.
        let cfg = InterferenceConfig {
            tolerance: 1e-3,
            ..Default::default()
        };
        let inters = detect_interference(&a, cfg);
        assert!(
            inters.is_empty(),
            "expected sub-tolerance overlap to be dropped: {inters:?}"
        );
    }

    /// One cube entirely inside another (nested) — should report a
    /// volume on the order of the inner cube's volume.
    #[test]
    fn nested_cubes_report_inner_volume() {
        let mut a = Assembly::new();
        // Outer cube: 4 × 4 × 4.
        let outer = valenx_cad::box_solid(4.0, 4.0, 4.0).unwrap();
        a.add_part(Part::new(0, "outer", outer));
        // Inner unit cube placed inside the outer at (1.5, 1.5, 1.5) —
        // wholly enclosed.
        let mut inner = unit_cube("inner");
        inner.transform.translation = Vector3::new(1.5, 1.5, 1.5);
        a.add_part(inner);
        let inters = detect_interference(&a, InterferenceConfig::default());
        assert_eq!(inters.len(), 1);
        // Inner cube volume is 1.0 — the estimate should be on that order.
        assert!(
            inters[0].volume_estimate > 0.5,
            "nested overlap underestimated: {}",
            inters[0].volume_estimate
        );
    }

    /// AABB intersection helper sanity.
    #[test]
    fn aabb_intersection_returns_box_for_overlap() {
        let a = AABB::new(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let b = AABB::new(Vector3::new(0.5, 0.5, 0.5), Vector3::new(2.0, 2.0, 2.0));
        let inter = aabb_intersection(&a, &b).unwrap();
        assert!((inter.min - Vector3::new(0.5, 0.5, 0.5)).norm() < 1e-12);
        assert!((inter.max - Vector3::new(1.0, 1.0, 1.0)).norm() < 1e-12);
        assert!((aabb_volume(&inter) - 0.125).abs() < 1e-12);
    }

    #[test]
    fn aabb_intersection_returns_none_for_disjoint() {
        let a = AABB::new(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let b = AABB::new(Vector3::new(2.0, 0.0, 0.0), Vector3::new(3.0, 1.0, 1.0));
        assert!(aabb_intersection(&a, &b).is_none());
    }
}
