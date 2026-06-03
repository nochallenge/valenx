//! Phase 144 — `ShapeAnalysis_Edge::CheckSameRange` — analyze
//! edge-to-vertex distance ratios to catch near-degenerate geometry.
//!
//! ## What OCCT does
//!
//! Computes the ratio `distance(edge_endpoint, vertex) / edge_length`
//! for every edge in a shape. When the ratio is large (> 1.0) the
//! edge "overshoots" its expected vertex, signalling a tolerance
//! mismatch between the geometric edge representation and the
//! topological vertex it claims to terminate at. Real OCCT uses this
//! to flag edges that round-tripped through a tolerant boolean
//! operation with too-coarse precision.
//!
//! ## v1 status
//!
//! **Honest v1** on a polyline-as-shape representation: caller
//! provides a list of (segment_a, segment_b) pairs giving an edge's
//! endpoints and the expected vertex coordinates; we compute the
//! distance ratio per pair and flag anything exceeding `max_ratio`.
//! The full BRep variant (consuming `truck_modeling::Edge`/`Vertex`)
//! is Phase 144.5 — depends on Phase 140.5's orientation-walk
//! infrastructure for the edge-vertex pairing.

use crate::error::OcctAdvancedError;

/// One (endpoint, claimed_vertex) pair to check.
///
/// `endpoint` is the geometric edge's actual endpoint coordinates;
/// `vertex` is the topological vertex it claims to terminate at;
/// `edge_length` is the edge's arc length (used as the ratio
/// denominator).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct EdgeVertexPair {
    /// Geometric edge endpoint coordinates.
    pub endpoint: [f64; 3],
    /// Topological vertex coordinates.
    pub vertex: [f64; 3],
    /// Edge arc length (used as the ratio denominator).
    pub edge_length: f64,
}

/// Report returned for a shape whose edges all pass the ratio check.
#[derive(Clone, Debug, PartialEq)]
pub struct DistRatioReport {
    /// Number of edge-vertex pairs analyzed.
    pub pairs: usize,
    /// Maximum ratio observed.
    pub max_ratio: f64,
    /// Mean ratio across all pairs.
    pub mean_ratio: f64,
}

/// Check that `pairs` all have distance ratios below `max_allowed`.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for empty input or non-positive
///   `max_allowed`.
/// - [`OcctAdvancedError::Defect`] when any pair exceeds the
///   threshold. Locus is the pair index.
pub fn shape_analysis_check_dist_ratio(
    pairs: &[EdgeVertexPair],
    max_allowed: f64,
) -> Result<DistRatioReport, OcctAdvancedError> {
    if pairs.is_empty() {
        return Err(OcctAdvancedError::bad_input(
            "pairs",
            "need at least one edge-vertex pair",
        ));
    }
    if !max_allowed.is_finite() || max_allowed <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "max_allowed",
            "must be positive finite",
        ));
    }

    let mut max_r = 0.0_f64;
    let mut sum = 0.0_f64;
    for (i, p) in pairs.iter().enumerate() {
        if !p.edge_length.is_finite() || p.edge_length <= 0.0 {
            return Err(OcctAdvancedError::bad_input(
                "edge_length",
                format!("pairs[{i}].edge_length must be positive finite"),
            ));
        }
        let dx = p.endpoint[0] - p.vertex[0];
        let dy = p.endpoint[1] - p.vertex[1];
        let dz = p.endpoint[2] - p.vertex[2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        let r = d / p.edge_length;
        if r > max_allowed {
            return Err(OcctAdvancedError::defect(
                format!("pairs[{i}]"),
                format!("distance ratio {r:.3e} > max {max_allowed:.3e}"),
            ));
        }
        max_r = max_r.max(r);
        sum += r;
    }

    Ok(DistRatioReport {
        pairs: pairs.len(),
        max_ratio: max_r,
        mean_ratio: sum / pairs.len() as f64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(endpoint: [f64; 3], vertex: [f64; 3], length: f64) -> EdgeVertexPair {
        EdgeVertexPair {
            endpoint,
            vertex,
            edge_length: length,
        }
    }

    #[test]
    fn rejects_empty_pairs() {
        let err = shape_analysis_check_dist_ratio(&[], 0.1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_max_allowed() {
        let p = [pair([0.0; 3], [0.0; 3], 1.0)];
        let err = shape_analysis_check_dist_ratio(&p, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_edge_length() {
        let p = [pair([0.0; 3], [0.0; 3], 0.0)];
        let err = shape_analysis_check_dist_ratio(&p, 0.1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn matched_endpoints_pass() {
        // endpoint == vertex => ratio 0.
        let p = [pair([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)];
        let r = shape_analysis_check_dist_ratio(&p, 0.1).unwrap();
        assert_eq!(r.pairs, 1);
        assert_eq!(r.max_ratio, 0.0);
    }

    #[test]
    fn small_drift_within_threshold() {
        // Drift 0.05 over edge length 1.0 → ratio 0.05 ≤ 0.1.
        let p = [pair([1.05, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)];
        let r = shape_analysis_check_dist_ratio(&p, 0.1).unwrap();
        assert!((r.max_ratio - 0.05).abs() < 1e-9);
    }

    #[test]
    fn large_drift_flagged() {
        // Drift 0.5 over edge length 1.0 → ratio 0.5 > 0.1.
        let p = [pair([1.5, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)];
        let err = shape_analysis_check_dist_ratio(&p, 0.1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.defect");
    }
}
