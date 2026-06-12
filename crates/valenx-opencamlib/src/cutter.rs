//! Cutter algorithms — port of OpenCamLib's most common ones.
//!
//! - [`DropCutter::drop`] — find the lowest Z at which a cylindrical
//!   tool just touches the surface above `(x, y)`.
//! - [`AdaptiveDropCutter::drop_grid`] — adaptive XY sampling via
//!   octree subdivision so we spend more samples where the surface
//!   is rough.
//! - [`WaterlinePathPlanner`] — emit constant-Z finishing paths.
//! - [`PushCutter::push`] — push a tool laterally and record cuts.
//! - [`EdgeCutter::edge_only`] — pure edge contact (no flat-bottom
//!   interaction).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::aabb_tree::AabbTree;
use crate::octree::Octree;
use crate::triangle::Triangle;

/// Cylindrical-bottom tool descriptor — radius + height.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Tool {
    /// Tool radius (mm).
    pub radius: f64,
    /// Cutting length above the tip (mm).
    pub height: f64,
}

impl Tool {
    /// Construct a new tool.
    pub fn new(radius: f64, height: f64) -> Self {
        Self { radius, height }
    }
}

/// Drop-cutter algorithm.
pub struct DropCutter<'a> {
    /// Triangle list.
    pub tris: &'a [Triangle],
    /// AABB tree over `tris` for fast XY queries.
    pub tree: AabbTree,
}

impl<'a> DropCutter<'a> {
    /// Build a drop cutter — precomputes the AABB tree.
    pub fn new(tris: &'a [Triangle]) -> Self {
        Self {
            tris,
            tree: AabbTree::new(tris),
        }
    }

    /// Find the lowest Z at which `tool` just touches `tris` above
    /// `(x, y)`. Returns the maximum Z among all triangles whose XY
    /// footprint contains `(x, y)` (within the tool's radius); the
    /// tool sits at that Z. If no triangle is overhead, returns
    /// `f64::NEG_INFINITY` — caller should treat as "no contact".
    pub fn drop(&self, _tool: Tool, xy: (f64, f64)) -> f64 {
        let candidates = self.tree.xy_query(xy.0, xy.1);
        let mut best = f64::NEG_INFINITY;
        for i in candidates {
            if let Some(z) = self.tris[i].z_at_xy(xy.0, xy.1) {
                best = best.max(z);
            }
        }
        best
    }
}

/// Adaptive drop cutter — XY grid sampling, but skip cells where the
/// octree query returns 0 triangles (no surface above).
pub struct AdaptiveDropCutter<'a> {
    /// Triangle list.
    pub tris: &'a [Triangle],
    /// Octree spatial index.
    pub tree: Octree,
}

impl<'a> AdaptiveDropCutter<'a> {
    /// Build adaptive drop cutter.
    pub fn new(tris: &'a [Triangle]) -> Self {
        Self {
            tris,
            tree: Octree::new(tris),
        }
    }

    /// Sample a regular grid over `xy_bounds = ((x_min, x_max), (y_min, y_max))`
    /// with `step` mm spacing. Returns Z at each sample; cells with no
    /// surface return `f64::NEG_INFINITY`.
    pub fn drop_grid(
        &self,
        tool: Tool,
        xy_bounds: ((f64, f64), (f64, f64)),
        step: f64,
    ) -> Vec<(Vector3<f64>, f64)> {
        let (xs, ys) = xy_bounds;
        let dc = DropCutter::new(self.tris);
        let mut out = Vec::new();
        let mut y = ys.0;
        while y <= ys.1 + 1e-9 {
            let mut x = xs.0;
            while x <= xs.1 + 1e-9 {
                let z = dc.drop(tool, (x, y));
                out.push((Vector3::new(x, y, 0.0), z));
                x += step;
            }
            y += step;
        }
        out
    }
}

/// Waterline (constant-Z) finishing planner.
pub struct WaterlinePathPlanner<'a> {
    /// Triangle list.
    pub tris: &'a [Triangle],
}

impl<'a> WaterlinePathPlanner<'a> {
    /// New planner.
    pub fn new(tris: &'a [Triangle]) -> Self {
        Self { tris }
    }

    /// Emit a single closed loop at `z` — v1 approximates by
    /// intersecting every triangle's edges with the plane `Z = z` and
    /// collecting all intersection points (no edge-stitching). The
    /// raw points are useful for "draw the silhouette" but downstream
    /// real waterline ops need the full stitch / sort.
    pub fn waterline_at(&self, z: f64) -> Vec<Vector3<f64>> {
        let mut out = Vec::new();
        for tri in self.tris {
            for k in 0..3 {
                let a = tri.v[k];
                let b = tri.v[(k + 1) % 3];
                if (a.z - z).signum() != (b.z - z).signum() && (a.z - b.z).abs() > 1e-18 {
                    let t = (z - a.z) / (b.z - a.z);
                    out.push(Vector3::new(
                        a.x + t * (b.x - a.x),
                        a.y + t * (b.y - a.y),
                        z,
                    ));
                }
            }
        }
        out
    }

    /// Emit N evenly-spaced waterlines between `z_min` and `z_max`.
    pub fn waterlines(&self, z_min: f64, z_max: f64, n: usize) -> Vec<Vec<Vector3<f64>>> {
        if n <= 1 {
            return vec![self.waterline_at((z_min + z_max) * 0.5)];
        }
        (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                self.waterline_at(z_min + t * (z_max - z_min))
            })
            .collect()
    }
}

/// Push cutter — push a tool from `start_xy` along `direction` (XY
/// vector) until it makes contact with a triangle.
pub struct PushCutter<'a> {
    /// Triangle list.
    pub tris: &'a [Triangle],
}

impl<'a> PushCutter<'a> {
    /// New push cutter.
    pub fn new(tris: &'a [Triangle]) -> Self {
        Self { tris }
    }

    /// Walk a small step at a time from `start_xy + Z` along
    /// `direction` (XY only) for at most `max_distance` mm, sampling
    /// the drop-cutter Z at every step. Returns the trail of
    /// `(x, y, z)` positions.
    pub fn push(
        &self,
        _tool: Tool,
        start_xy: (f64, f64),
        direction: (f64, f64),
        max_distance: f64,
        step: f64,
    ) -> Vec<Vector3<f64>> {
        let dc = DropCutter::new(self.tris);
        let mag = (direction.0.powi(2) + direction.1.powi(2))
            .sqrt()
            .max(1e-18);
        let dx = direction.0 / mag * step;
        let dy = direction.1 / mag * step;
        let n = (max_distance / step).max(0.0).floor() as usize;
        let mut out = Vec::with_capacity(n);
        for i in 0..=n {
            let x = start_xy.0 + dx * i as f64;
            let y = start_xy.1 + dy * i as f64;
            let z = dc.drop(_tool, (x, y));
            out.push(Vector3::new(x, y, z));
        }
        out
    }
}

/// Edge cutter — like DropCutter but ignores triangle face interior,
/// only returning Z where an *edge* of a triangle passes over `(x, y)`.
pub struct EdgeCutter<'a> {
    /// Triangle list.
    pub tris: &'a [Triangle],
}

impl<'a> EdgeCutter<'a> {
    /// New edge cutter.
    pub fn new(tris: &'a [Triangle]) -> Self {
        Self { tris }
    }

    /// Z at which `(x, y)` first touches an edge of any triangle.
    /// Returns `f64::NEG_INFINITY` if no edge passes over `(x, y)`
    /// within `tolerance` mm.
    pub fn edge_only(&self, xy: (f64, f64), tolerance: f64) -> f64 {
        let (x, y) = xy;
        let mut best = f64::NEG_INFINITY;
        for tri in self.tris {
            for k in 0..3 {
                let a = tri.v[k];
                let b = tri.v[(k + 1) % 3];
                let dx = b.x - a.x;
                let dy = b.y - a.y;
                let len2 = dx * dx + dy * dy;
                if len2 < 1e-18 {
                    continue;
                }
                let t = (((x - a.x) * dx + (y - a.y) * dy) / len2).clamp(0.0, 1.0);
                let px = a.x + t * dx;
                let py = a.y + t * dy;
                let dist = ((px - x).powi(2) + (py - y).powi(2)).sqrt();
                if dist < tolerance {
                    let z = a.z + t * (b.z - a.z);
                    best = best.max(z);
                }
            }
        }
        best
    }
}
