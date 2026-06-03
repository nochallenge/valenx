//! 2D Voronoi via the dual of a Delaunay triangulation (v1).
//!
//! CGAL ships a Fortune's-algorithm implementation, but a Delaunay
//! dual produces equivalent cell topology for our v1 needs — every
//! Voronoi vertex is the circumcentre of a Delaunay triangle, and
//! every cell's edges connect circumcentres of triangles sharing an
//! edge.

use crate::delaunay::triangulate_2d;
use crate::error::CgalError;

/// A single Voronoi cell — site point + ordered cell vertices.
///
/// `vertices` may be partial / open along the convex hull (cells on
/// the boundary of the input have unbounded edges in the true Voronoi
/// diagram). v1 just returns the finite circumcentre vertices it
/// found; callers needing bounded regions should clip against the
/// input AABB.
#[derive(Clone, Debug)]
pub struct VoronoiCell {
    /// Site coordinate.
    pub site: [f64; 2],
    /// Ordered cell vertices (best-effort, may be partial).
    pub vertices: Vec<[f64; 2]>,
}

/// Compute the 2D Voronoi diagram (dual of Delaunay).
pub fn voronoi_2d(points: &[[f64; 2]]) -> Result<Vec<VoronoiCell>, CgalError> {
    let tris = triangulate_2d(points)?;
    // For each site (input point), collect circumcentres of incident
    // triangles.
    let mut cells: Vec<VoronoiCell> = points
        .iter()
        .map(|p| VoronoiCell {
            site: *p,
            vertices: Vec::new(),
        })
        .collect();
    for t in &tris {
        let cc = circumcentre(points[t.0], points[t.1], points[t.2]);
        if let Some(cc) = cc {
            cells[t.0].vertices.push(cc);
            cells[t.1].vertices.push(cc);
            cells[t.2].vertices.push(cc);
        }
    }
    // Sort each cell's vertices CCW around the site so callers can
    // emit a polygon.
    for cell in &mut cells {
        let cx = cell.site[0];
        let cy = cell.site[1];
        cell.vertices.sort_by(|a, b| {
            let aa = (a[1] - cy).atan2(a[0] - cx);
            let bb = (b[1] - cy).atan2(b[0] - cx);
            aa.partial_cmp(&bb).unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    Ok(cells)
}

fn circumcentre(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<[f64; 2]> {
    let ax = b[0] - a[0];
    let ay = b[1] - a[1];
    let bx = c[0] - a[0];
    let by = c[1] - a[1];
    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-18 {
        return None;
    }
    let ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    let uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    Some([a[0] + ux, a[1] + uy])
}
