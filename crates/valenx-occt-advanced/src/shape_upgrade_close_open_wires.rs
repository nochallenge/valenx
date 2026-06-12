//! Phase 151 — `ShapeFix_Wire::FixClosed` — close open wires by
//! adding edges.
//!
//! ## What OCCT does
//!
//! `ShapeFix_Wire(wire, ...).FixClosed(tolerance)` walks an open
//! polyline wire and synthesises new edges to close it. The closure
//! strategy:
//!
//! 1. If the wire's endpoints are within `tolerance`, fuse them
//!    (merge the front vertex of the first edge with the back vertex
//!    of the last edge).
//! 2. Otherwise insert a new straight edge from `back` to `front`.
//! 3. Re-orient any edges whose direction is inconsistent with the
//!    walk-around-the-loop direction.
//!
//! Required cleanup after a partial boolean leaves the boundary of a
//! face as an almost-closed wire missing the connecting bit.
//!
//! ## v1 status
//!
//! **Honest v1** for the polyline-vertex representation: caller
//! provides `Vec<[f64; 3]>` and we synthesise the closing point if
//! needed (within tolerance: snap last to first; otherwise append a
//! new vertex equal to the first). The full BRep-level variant
//! (consuming `truck_modeling::Wire`, returning a closed wire) is
//! Phase 151.5 — depends on truck exposing wire mutation.
//!
//! ## Phase 151.5 — arc closure variant
//!
//! [`shape_upgrade_close_open_wires_arc`] is the graduated arc
//! variant: instead of a straight bridge, the open wire is closed
//! with a **circular arc** fitted through three points — the wire's
//! last vertex, a caller-supplied through-point, and the wire's
//! first vertex. The unique circle through those three points is
//! computed (the 3D circumcircle), and the arc from `last` to
//! `first` passing through the through-point is sampled into the
//! requested number of straight segments. This is the closure CAD
//! sketchers want when the boundary should rejoin with a fillet-like
//! curve rather than a hard corner.

use crate::error::OcctAdvancedError;

/// Result of the close-open-wire op.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedWire {
    /// Closed polyline (first and last vertex are the same point
    /// within tolerance — or exactly equal after snap).
    pub vertices: Vec<[f64; 3]>,
    /// Action taken: `Snap` when the wire was already nearly closed,
    /// `Bridge` when a new edge was synthesised.
    pub action: CloseAction,
}

/// What [`shape_upgrade_close_open_wires`] did to close the wire.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum CloseAction {
    /// The wire was already closed within tolerance — last vertex
    /// was snapped to match the first exactly.
    Snap,
    /// The wire's endpoints were too far apart for snap — a new edge
    /// (last → first) was synthesised by appending the first vertex
    /// to the end.
    Bridge,
}

/// Close an open polyline wire.
///
/// `vertices` — the open wire's vertex sequence.
/// `tolerance` — distance below which the endpoints count as "the
/// same point" and get snapped (vs bridged).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for fewer than 2 vertices or
///   non-positive tolerance.
pub fn shape_upgrade_close_open_wires(
    vertices: &[[f64; 3]],
    tolerance: f64,
) -> Result<ClosedWire, OcctAdvancedError> {
    if vertices.len() < 2 {
        return Err(OcctAdvancedError::bad_input(
            "vertices",
            "need ≥2 vertices to close",
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "tolerance",
            "must be positive finite",
        ));
    }

    let first = vertices[0];
    let last = *vertices.last().unwrap();
    let dx = last[0] - first[0];
    let dy = last[1] - first[1];
    let dz = last[2] - first[2];
    let gap = (dx * dx + dy * dy + dz * dz).sqrt();

    let mut out = vertices.to_vec();
    let action = if gap < tolerance {
        // Snap: replace the last vertex with first.
        *out.last_mut().unwrap() = first;
        CloseAction::Snap
    } else {
        // Bridge: append first to close the loop.
        out.push(first);
        CloseAction::Bridge
    };

    Ok(ClosedWire {
        vertices: out,
        action,
    })
}

/// Close an open polyline wire with a **circular arc** (Phase 151.5).
///
/// The closing arc is the unique circle through three points — the
/// wire's last vertex, the caller-supplied `arc_through` point, and
/// the wire's first vertex — sampled from `last` to `first` (passing
/// through `arc_through`) into `segments` straight chords. The arc's
/// interior sample points are appended after the original vertices;
/// the final appended point equals the wire's first vertex, closing
/// the loop.
///
/// `segments` is the number of chord segments the arc is sampled
/// into (`segments >= 1`; `1` degenerates to a straight bridge).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for fewer than 2 vertices,
///   `segments == 0`, a non-finite `arc_through`, or three points
///   that are collinear (no unique circle exists — use the straight
///   [`shape_upgrade_close_open_wires`] for that case).
///
/// # Example
///
/// ```
/// use valenx_occt_advanced::shape_upgrade_close_open_wires::shape_upgrade_close_open_wires_arc;
/// // Open wire on the X axis; close it with an arc bulging through
/// // (0.5, 1, 0).
/// let wire = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
/// let closed = shape_upgrade_close_open_wires_arc(&wire, [0.5, 1.0, 0.0], 8).unwrap();
/// assert_eq!(closed.vertices.first().unwrap(), &[0.0, 0.0, 0.0]);
/// assert_eq!(closed.vertices.last().unwrap(), &[0.0, 0.0, 0.0]);
/// ```
pub fn shape_upgrade_close_open_wires_arc(
    vertices: &[[f64; 3]],
    arc_through: [f64; 3],
    segments: usize,
) -> Result<ClosedWire, OcctAdvancedError> {
    if vertices.len() < 2 {
        return Err(OcctAdvancedError::bad_input(
            "vertices",
            "need ≥2 vertices to close",
        ));
    }
    if segments == 0 {
        return Err(OcctAdvancedError::bad_input("segments", "must be ≥1"));
    }
    if arc_through.iter().any(|c| !c.is_finite()) {
        return Err(OcctAdvancedError::bad_input(
            "arc_through",
            "must be finite",
        ));
    }
    let first = vertices[0];
    let last = *vertices.last().unwrap();

    // Sample the circular arc last -> arc_through -> first.
    let arc = sample_arc(last, arc_through, first, segments).ok_or_else(|| {
        OcctAdvancedError::bad_input(
            "arc_through",
            "last, arc_through and first are collinear — no unique arc",
        )
    })?;

    // arc[0] == last (already in `vertices`); append arc[1..] which
    // ends exactly at `first`, closing the loop.
    let mut out = vertices.to_vec();
    out.extend_from_slice(&arc[1..]);

    Ok(ClosedWire {
        vertices: out,
        // The arc closure always synthesises new geometry — it is a
        // (curved) bridge.
        action: CloseAction::Bridge,
    })
}

/// Sample the circular arc through three points `p0 -> p1 -> p2`,
/// returning `segments + 1` points from `p0` to `p2` (inclusive).
/// Returns `None` when the three points are collinear.
fn sample_arc(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], segments: usize) -> Option<Vec<[f64; 3]>> {
    // Plane of the three points.
    let v01 = sub(p1, p0);
    let v02 = sub(p2, p0);
    let normal = cross(v01, v02);
    let n_len = norm(normal);
    if n_len < 1e-12 {
        return None; // collinear
    }
    let n = scale(normal, 1.0 / n_len);

    // Circumcentre: intersection of the perpendicular bisectors of
    // the chords p0-p1 and p0-p2, solved in the triangle's plane.
    // Centre = p0 + α·v01 + β·v02 with α, β from the 2x2 system
    //   [v01·v01  v01·v02] [α]   [0.5·v01·v01]
    //   [v01·v02  v02·v02] [β] = [0.5·v02·v02]
    let a = dot(v01, v01);
    let b = dot(v01, v02);
    let c = dot(v02, v02);
    let det = a * c - b * b;
    if det.abs() < 1e-18 {
        return None;
    }
    let rhs0 = 0.5 * a;
    let rhs1 = 0.5 * c;
    let alpha = (rhs0 * c - rhs1 * b) / det;
    let beta = (a * rhs1 - b * rhs0) / det;
    let centre = add(p0, add(scale(v01, alpha), scale(v02, beta)));
    let radius = norm(sub(p0, centre));
    if radius < 1e-12 {
        return None;
    }

    // In-plane basis with `ex` toward p0.
    let ex = scale(sub(p0, centre), 1.0 / radius);
    let ey = cross(n, ex); // unit (n ⟂ ex, both unit)

    // Angle of each point about the centre in (ex, ey).
    let ang = |p: [f64; 3]| -> f64 {
        let d = sub(p, centre);
        dot(d, ey).atan2(dot(d, ex))
    };
    let a0 = ang(p0); // == 0 by construction
    let mut a1 = ang(p1);
    let mut a2 = ang(p2);
    // Unwrap so the sweep p0 -> p1 -> p2 is monotonic. Walk in the
    // direction that passes through p1.
    let two_pi = std::f64::consts::TAU;
    // Normalise relative to a0.
    while a1 - a0 <= 0.0 {
        a1 += two_pi;
    }
    while a2 - a0 <= 0.0 {
        a2 += two_pi;
    }
    // The arc must contain p1: if p2 comes "before" p1, the sweep
    // actually goes the long way — add a turn so p2 follows p1.
    if a2 < a1 {
        a2 += two_pi;
    }

    let mut out = Vec::with_capacity(segments + 1);
    for s in 0..=segments {
        let t = s as f64 / segments as f64;
        let theta = a0 + (a2 - a0) * t;
        out.push(add(
            centre,
            add(
                scale(ex, radius * theta.cos()),
                scale(ey, radius * theta.sin()),
            ),
        ));
    }
    // Pin the endpoints exactly (guard against float drift).
    out[0] = p0;
    let last = out.len() - 1;
    out[last] = p2;
    Some(out)
}

// --- vector helpers ---
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_single_vertex() {
        let err = shape_upgrade_close_open_wires(&[[0.0; 3]], 1e-6).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_tolerance() {
        let v = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let err = shape_upgrade_close_open_wires(&v, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn snap_when_almost_closed() {
        // Almost-closed triangle — last vertex within tolerance of
        // first.
        let v = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1e-9, 1e-9, 0.0],
        ];
        let r = shape_upgrade_close_open_wires(&v, 1e-6).unwrap();
        assert_eq!(r.action, CloseAction::Snap);
        assert_eq!(r.vertices.len(), 4);
        assert_eq!(r.vertices.last().unwrap(), &[0.0, 0.0, 0.0]);
    }

    #[test]
    fn bridge_when_endpoints_far_apart() {
        // Open triangle, large gap → bridge by appending first.
        let v = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let r = shape_upgrade_close_open_wires(&v, 1e-6).unwrap();
        assert_eq!(r.action, CloseAction::Bridge);
        assert_eq!(r.vertices.len(), 4);
        assert_eq!(r.vertices.last().unwrap(), &v[0]);
    }

    // --- Phase 151.5 arc-closure tests ---

    #[test]
    fn arc_close_appends_arc_segments_and_closes_loop() {
        // Open wire on the X axis; close with an arc bulging up.
        let wire = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let r = shape_upgrade_close_open_wires_arc(&wire, [1.0, 1.0, 0.0], 8).unwrap();
        assert_eq!(r.action, CloseAction::Bridge);
        // 2 original + 8 arc samples (arc[1..] of 9) = 10.
        assert_eq!(r.vertices.len(), 10);
        // Loop closes exactly at the wire's first vertex.
        assert_eq!(r.vertices.first().unwrap(), &[0.0, 0.0, 0.0]);
        assert_eq!(r.vertices.last().unwrap(), &[0.0, 0.0, 0.0]);
    }

    #[test]
    fn arc_close_samples_lie_on_the_fitted_circle() {
        // last=(2,0,0), through=(1,1,0), first=(0,0,0): the unique
        // circle is centred at (1,0,0) with radius 1. Every arc
        // sample must be radius 1 from that centre.
        let wire = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let r = shape_upgrade_close_open_wires_arc(&wire, [1.0, 1.0, 0.0], 16).unwrap();
        let centre = [1.0, 0.0, 0.0];
        for p in &r.vertices[1..] {
            let d = (((p[0] - centre[0]).powi(2)
                + (p[1] - centre[1]).powi(2)
                + (p[2] - centre[2]).powi(2))
            .sqrt()
                - 1.0)
                .abs();
            assert!(d < 1e-9, "arc sample {p:?} is not on the unit circle");
        }
        // The arc bulges through y>0 (passes near the through-point).
        assert!(r.vertices.iter().any(|p| p[1] > 0.9));
    }

    #[test]
    fn arc_close_rejects_collinear_through_point() {
        // last, through, first all on the X axis → no unique circle.
        let wire = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let err = shape_upgrade_close_open_wires_arc(&wire, [1.0, 0.0, 0.0], 8).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn arc_close_rejects_zero_segments() {
        let wire = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let err = shape_upgrade_close_open_wires_arc(&wire, [1.0, 1.0, 0.0], 0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn arc_close_rejects_short_wire() {
        let err = shape_upgrade_close_open_wires_arc(&[[0.0; 3]], [1.0, 1.0, 0.0], 4).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }
}
