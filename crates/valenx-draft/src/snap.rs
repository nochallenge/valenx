//! Snap helpers for the Draft workbench.
//!
//! Pure-data routines the UI layer calls to:
//!
//! - collect all snap-relevant vertices in a document
//!   ([`endpoints`], [`midpoints`]),
//! - pick the closest snap candidate to a target cursor position
//!   ([`nearest`]),
//! - round a free cursor position onto a grid ([`grid_snap`]).
//!
//! All helpers operate on 2D points in the local frame of the
//! document's [`super::WorkingPlane`].

use crate::document::DraftDocument;
use crate::entity::DraftEntity;

/// Collect every endpoint visible in the document.
///
/// Includes:
/// - line start + end,
/// - polyline vertices,
/// - arc start / end positions (centre + radius·(cos, sin) on each angle),
/// - rectangle corners,
/// - polygon vertices,
/// - linear-dimension `from` / `to`,
/// - circle centre + 4 cardinal perimeter points (top / bottom / left / right).
pub fn endpoints(doc: &DraftDocument) -> Vec<[f64; 2]> {
    let mut out = Vec::new();
    for e in &doc.entities {
        match e {
            DraftEntity::Line { start, end } => {
                out.push(*start);
                out.push(*end);
            }
            DraftEntity::Polyline { points, .. } => {
                out.extend_from_slice(points);
            }
            DraftEntity::Arc {
                center,
                radius,
                start_angle,
                end_angle,
            } => {
                out.push([
                    center[0] + radius * start_angle.cos(),
                    center[1] + radius * start_angle.sin(),
                ]);
                out.push([
                    center[0] + radius * end_angle.cos(),
                    center[1] + radius * end_angle.sin(),
                ]);
                out.push(*center);
            }
            DraftEntity::Circle { center, radius } => {
                out.push(*center);
                out.push([center[0] + radius, center[1]]);
                out.push([center[0] - radius, center[1]]);
                out.push([center[0], center[1] + radius]);
                out.push([center[0], center[1] - radius]);
            }
            DraftEntity::Rectangle { min, max } => {
                out.push(*min);
                out.push([max[0], min[1]]);
                out.push(*max);
                out.push([min[0], max[1]]);
            }
            DraftEntity::Polygon {
                center,
                radius,
                sides,
            } => {
                let n = *sides as i32;
                if n >= 1 {
                    let two_pi = std::f64::consts::TAU;
                    for i in 0..n {
                        let theta = two_pi * (i as f64) / (n as f64);
                        out.push([
                            center[0] + radius * theta.cos(),
                            center[1] + radius * theta.sin(),
                        ]);
                    }
                }
            }
            DraftEntity::LinearDimension { from, to, .. } => {
                out.push(*from);
                out.push(*to);
            }
            DraftEntity::Text { position, .. } => {
                out.push(*position);
            }
        }
    }
    out
}

/// Collect every midpoint visible in the document.
///
/// Includes:
/// - line midpoint,
/// - midpoint of every segment in a polyline (closing segment too
///   when `closed`),
/// - midpoint of every rectangle edge,
/// - midpoint of every polygon edge.
///
/// Arcs / circles are not included — their geometric "midpoint" is
/// ambiguous (centre vs. point-on-arc); the perimeter snap handles
/// those separately if a future revision needs them.
pub fn midpoints(doc: &DraftDocument) -> Vec<[f64; 2]> {
    let mut out = Vec::new();
    let mid = |a: [f64; 2], b: [f64; 2]| [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    for e in &doc.entities {
        match e {
            DraftEntity::Line { start, end } => {
                out.push(mid(*start, *end));
            }
            DraftEntity::Polyline { points, closed } => {
                for w in points.windows(2) {
                    out.push(mid(w[0], w[1]));
                }
                if *closed && points.len() >= 2 {
                    out.push(mid(*points.last().unwrap(), points[0]));
                }
            }
            DraftEntity::Rectangle { min, max } => {
                let bl = *min;
                let br = [max[0], min[1]];
                let tr = *max;
                let tl = [min[0], max[1]];
                out.push(mid(bl, br));
                out.push(mid(br, tr));
                out.push(mid(tr, tl));
                out.push(mid(tl, bl));
            }
            DraftEntity::Polygon {
                center,
                radius,
                sides,
            } => {
                let n = *sides as i32;
                if n >= 2 {
                    let two_pi = std::f64::consts::TAU;
                    let v = |i: i32| -> [f64; 2] {
                        let theta = two_pi * (i as f64) / (n as f64);
                        [
                            center[0] + radius * theta.cos(),
                            center[1] + radius * theta.sin(),
                        ]
                    };
                    for i in 0..n {
                        out.push(mid(v(i), v((i + 1) % n)));
                    }
                }
            }
            _ => {}
        }
    }
    out
}

/// Pick the closest candidate to `target`, but only if its distance
/// is `<= threshold`. `None` when no candidate is within range.
pub fn nearest(target: [f64; 2], candidates: &[[f64; 2]], threshold: f64) -> Option<[f64; 2]> {
    let mut best: Option<([f64; 2], f64)> = None;
    for c in candidates {
        let dx = c[0] - target[0];
        let dy = c[1] - target[1];
        let d2 = dx * dx + dy * dy;
        match best {
            None => best = Some((*c, d2)),
            Some((_, bd2)) if d2 < bd2 => best = Some((*c, d2)),
            _ => {}
        }
    }
    let (p, d2) = best?;
    if d2.sqrt() <= threshold {
        Some(p)
    } else {
        None
    }
}

/// Round a free cursor coordinate onto a regular grid with the given
/// `spacing`. Each axis is independently snapped to the nearest
/// multiple of `spacing`. Returns `point` unchanged when
/// `spacing <= 0`.
pub fn grid_snap(point: [f64; 2], spacing: f64) -> [f64; 2] {
    if spacing <= 0.0 {
        return point;
    }
    [
        (point[0] / spacing).round() * spacing,
        (point[1] / spacing).round() * spacing,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plane::WorkingPlane;

    #[test]
    fn endpoints_collects_line_and_circle_points() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        d.add_entity(DraftEntity::Line {
            start: [0.0, 0.0],
            end: [10.0, 0.0],
        });
        d.add_entity(DraftEntity::Circle {
            center: [5.0, 5.0],
            radius: 2.0,
        });
        let pts = endpoints(&d);
        assert!(pts.contains(&[0.0, 0.0]));
        assert!(pts.contains(&[10.0, 0.0]));
        assert!(pts.contains(&[5.0, 5.0]));
        assert!(pts.contains(&[7.0, 5.0]));
        assert!(pts.contains(&[3.0, 5.0]));
    }

    #[test]
    fn endpoints_collects_polygon_vertices() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        d.add_entity(DraftEntity::Polygon {
            center: [0.0, 0.0],
            radius: 1.0,
            sides: 4,
        });
        let pts = endpoints(&d);
        // Square: vertices at (1,0), (0,1), (-1,0), (0,-1).
        assert!(pts
            .iter()
            .any(|p| (p[0] - 1.0).abs() < 1e-9 && p[1].abs() < 1e-9));
        assert!(pts
            .iter()
            .any(|p| p[0].abs() < 1e-9 && (p[1] - 1.0).abs() < 1e-9));
    }

    #[test]
    fn midpoints_line_is_average() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        d.add_entity(DraftEntity::Line {
            start: [0.0, 0.0],
            end: [4.0, 6.0],
        });
        assert_eq!(midpoints(&d), vec![[2.0, 3.0]]);
    }

    #[test]
    fn midpoints_polyline_closed_emits_closing_segment_midpoint() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        d.add_entity(DraftEntity::Polyline {
            points: vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0]],
            closed: true,
        });
        // open polyline mids: (1,0), (2,1); closing segment mid: (1,1).
        let mids = midpoints(&d);
        assert!(mids.contains(&[1.0, 0.0]));
        assert!(mids.contains(&[2.0, 1.0]));
        assert!(mids.contains(&[1.0, 1.0]));
    }

    #[test]
    fn midpoints_polyline_open_does_not_close() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        d.add_entity(DraftEntity::Polyline {
            points: vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0]],
            closed: false,
        });
        let mids = midpoints(&d);
        assert!(mids.contains(&[1.0, 0.0]));
        assert!(mids.contains(&[2.0, 1.0]));
        assert!(!mids.contains(&[1.0, 1.0]));
    }

    #[test]
    fn nearest_picks_closest_within_threshold() {
        let cands = vec![[0.0, 0.0], [10.0, 0.0], [3.0, 4.0]];
        // distance 0 → exact hit
        assert_eq!(nearest([0.0, 0.0], &cands, 0.1), Some([0.0, 0.0]));
        // closest is (3,4) at distance 5; threshold 6 → hit.
        assert_eq!(nearest([3.0, 4.0], &cands, 6.0), Some([3.0, 4.0]));
        // closest is (3,4) at distance ~1; threshold 0.5 → miss.
        assert_eq!(nearest([3.5, 4.0], &cands, 0.2), None);
    }

    #[test]
    fn nearest_returns_none_on_empty_candidates() {
        assert_eq!(nearest([0.0, 0.0], &[], 1.0), None);
    }

    #[test]
    fn grid_snap_rounds_each_axis() {
        assert_eq!(grid_snap([1.2, 3.7], 1.0), [1.0, 4.0]);
        assert_eq!(grid_snap([0.49, 0.51], 0.5), [0.5, 0.5]);
        // spacing <= 0 returns input unchanged
        assert_eq!(grid_snap([1.234, 5.678], 0.0), [1.234, 5.678]);
        assert_eq!(grid_snap([1.234, 5.678], -1.0), [1.234, 5.678]);
    }
}
