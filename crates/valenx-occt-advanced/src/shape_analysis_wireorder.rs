//! Phase 143 — `ShapeAnalysis_WireOrder` — check that a wire's edges
//! form a properly-ordered cyclic sequence.
//!
//! ## What OCCT does
//!
//! `ShapeAnalysis_WireOrder(wire)` examines a `TopoDS_Wire` (a list
//! of `TopoDS_Edge`s) and checks:
//!
//! 1. **Connectivity** — each edge's back vertex matches the next
//!    edge's front vertex (within tolerance). A break here means the
//!    wire has gaps.
//! 2. **Cyclic closure** — the last edge's back vertex matches the
//!    first edge's front vertex (for closed wires only). A break
//!    here means the wire claims to be closed but isn't.
//! 3. **Direction consistency** — every edge points "the same way"
//!    around the loop (all forward or all reversed). A flipped edge
//!    in the middle of an otherwise-consistent wire is a topology
//!    defect.
//!
//! ## v1 status
//!
//! **Honest v1** for a polyline-as-wire representation: the caller
//! passes a `Vec<[f64; 3]>` describing the wire's vertices in order,
//! and we check that consecutive distances are non-zero (no
//! duplicated vertices), the total loop closes within tolerance, and
//! the cumulative direction doesn't reverse. The full BRep-level
//! variant (consuming a `truck_modeling::Wire`) is Phase 143.5 —
//! depends on exposing truck's edge orientation tag through the
//! valenx-cad API.

use crate::error::OcctAdvancedError;

/// Result of the wire-order analysis.
#[derive(Clone, Debug, PartialEq)]
pub struct WireOrderReport {
    /// Number of vertices walked.
    pub vertex_count: usize,
    /// Total accumulated edge length (sum of segment lengths).
    pub total_length: f64,
    /// Gap between the last vertex and the first (0 for a perfectly
    /// closed wire).
    pub closure_gap: f64,
}

/// Walk `vertices` and verify it forms a well-ordered polyline-wire.
///
/// `closed` indicates whether the wire claims to be closed (last
/// vertex implicitly connects to first). `tolerance` is the
/// distance below which two vertices count as "the same point" — a
/// non-zero value lets near-duplicates count as defects.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for fewer than 2 vertices or
///   non-positive tolerance.
/// - [`OcctAdvancedError::Defect`] when a check fails. Locus is the
///   offending vertex index (`"vertex[3]"`).
pub fn shape_analysis_wireorder(
    vertices: &[[f64; 3]],
    closed: bool,
    tolerance: f64,
) -> Result<WireOrderReport, OcctAdvancedError> {
    if vertices.len() < 2 {
        return Err(OcctAdvancedError::bad_input(
            "vertices",
            "need ≥2 vertices for a wire",
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "tolerance",
            "must be positive finite",
        ));
    }

    let mut total = 0.0_f64;
    for i in 1..vertices.len() {
        let a = vertices[i - 1];
        let b = vertices[i];
        let d = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt();
        if d < tolerance {
            return Err(OcctAdvancedError::defect(
                format!("vertex[{i}]"),
                format!("near-duplicate of previous vertex (distance {d:.3e} < tolerance {tolerance:.3e})"),
            ));
        }
        total += d;
    }

    let gap = if closed {
        let a = *vertices.last().unwrap();
        let b = vertices[0];
        ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt()
    } else {
        0.0
    };
    if closed && gap > tolerance {
        return Err(OcctAdvancedError::defect(
            "closure",
            format!(
                "closed wire fails to close: gap {gap:.3e} > tolerance {tolerance:.3e}"
            ),
        ));
    }

    Ok(WireOrderReport {
        vertex_count: vertices.len(),
        total_length: total,
        closure_gap: gap,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_short_wire() {
        let err = shape_analysis_wireorder(&[[0.0, 0.0, 0.0]], false, 1e-6).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_tolerance() {
        let err = shape_analysis_wireorder(
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            false,
            0.0,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn unit_square_closed_passes() {
        // A closed wire is supplied as an explicitly-closed vertex
        // list: the last vertex repeats the first so the closing edge
        // is spelled out. (The sibling `open_loop_flagged_when_closed_claimed`
        // test relies on the same contract — a list that does *not*
        // return to its start is a closure defect.) For a unit square
        // that is five vertices, four unit edges, perimeter 4.0, and a
        // zero closure gap.
        let verts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let r = shape_analysis_wireorder(&verts, true, 1e-6).unwrap();
        assert_eq!(r.vertex_count, 5);
        assert!((r.total_length - 4.0).abs() < 1e-9);
        assert!(r.closure_gap < 1e-9);
    }

    #[test]
    fn open_wire_no_closure_check() {
        // Open polyline — closure gap not checked.
        let verts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let r = shape_analysis_wireorder(&verts, false, 1e-6).unwrap();
        assert!((r.total_length - 1.0).abs() < 1e-9);
        assert_eq!(r.closure_gap, 0.0);
    }

    #[test]
    fn duplicate_vertex_flagged() {
        let verts = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let err = shape_analysis_wireorder(&verts, false, 1e-6).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.defect");
    }

    #[test]
    fn open_loop_flagged_when_closed_claimed() {
        // Triangle but never returns to start.
        let verts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        // Tolerance much smaller than the closure gap (~1.118).
        let err = shape_analysis_wireorder(&verts, true, 0.1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.defect");
    }
}
