//! Inter-atomic distances and a grid-accelerated neighbour search.
//!
//! Pairwise distance is trivial; the value here is [`NeighborGrid`], a
//! uniform spatial hash that turns "every atom within `r` of a query
//! point" from an `O(n)` scan into an `O(1)`-expected lookup. Contact
//! maps, clash detection and SASA all build on it.

use crate::error::{BiostructError, Result};
use crate::structure::Model;
use nalgebra::Point3;
use std::collections::HashMap;

/// Euclidean distance between two points, in ångström.
pub fn distance(a: &Point3<f64>, b: &Point3<f64>) -> f64 {
    (a - b).norm()
}

/// Squared distance — cheaper when only a comparison is needed.
pub fn distance_sq(a: &Point3<f64>, b: &Point3<f64>) -> f64 {
    (a - b).norm_squared()
}

/// A uniform-grid spatial hash over a fixed point set.
///
/// The grid cell size equals the search cutoff, so a neighbour query
/// only ever inspects the 27 cells of a 3×3×3 block around the query
/// point. Build once, query many times.
#[derive(Debug, Clone)]
pub struct NeighborGrid {
    /// All indexed points, in insertion order.
    points: Vec<Point3<f64>>,
    /// `cell -> point indices`.
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
    /// Cell edge length (== the cutoff the grid was built for).
    cell_size: f64,
    /// Lower corner of the bounding box (grid origin).
    origin: Point3<f64>,
}

impl NeighborGrid {
    /// Build a grid over `points` sized for queries up to `cutoff`
    /// ångström. `cutoff` must be strictly positive.
    pub fn new(points: &[Point3<f64>], cutoff: f64) -> Result<NeighborGrid> {
        if cutoff <= 0.0 || cutoff.is_nan() {
            return Err(BiostructError::invalid(
                "cutoff",
                "neighbour-grid cutoff must be positive",
            ));
        }
        let origin = if points.is_empty() {
            Point3::origin()
        } else {
            let mut lo = points[0];
            for p in points {
                lo.x = lo.x.min(p.x);
                lo.y = lo.y.min(p.y);
                lo.z = lo.z.min(p.z);
            }
            lo
        };
        let mut cells: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for (i, p) in points.iter().enumerate() {
            let key = cell_key(p, &origin, cutoff);
            cells.entry(key).or_default().push(i);
        }
        Ok(NeighborGrid {
            points: points.to_vec(),
            cells,
            cell_size: cutoff,
            origin,
        })
    }

    /// Build a grid over every atom of `model`'s first-model-style
    /// coordinate set. Index `i` corresponds to the `i`-th atom in
    /// [`Model::atoms`] order.
    pub fn from_model(model: &Model, cutoff: f64) -> Result<NeighborGrid> {
        let pts: Vec<Point3<f64>> = model.atoms().map(|a| a.coord).collect();
        NeighborGrid::new(&pts, cutoff)
    }

    /// Number of indexed points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the grid holds no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Indices of every point within `radius` of `query`. `radius`
    /// must not exceed the cutoff the grid was built for, or the
    /// 27-cell window can miss a neighbour — the method clamps it and
    /// reports nothing beyond `cell_size`.
    pub fn within(&self, query: &Point3<f64>, radius: f64) -> Vec<usize> {
        let r = radius.min(self.cell_size);
        let r2 = r * r;
        let (cx, cy, cz) = cell_key(query, &self.origin, self.cell_size);
        let mut out = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(indices) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &i in indices {
                            if distance_sq(query, &self.points[i]) <= r2 {
                                out.push(i);
                            }
                        }
                    }
                }
            }
        }
        out
    }

    /// Like [`within`](Self::within) but excludes the point whose
    /// index equals `self_index` — for "neighbours of atom *i*"
    /// queries that should not return *i* itself.
    pub fn within_excluding(
        &self,
        query: &Point3<f64>,
        radius: f64,
        self_index: usize,
    ) -> Vec<usize> {
        self.within(query, radius)
            .into_iter()
            .filter(|&i| i != self_index)
            .collect()
    }

    /// Every unordered index pair `(i, j)` with `i < j` closer than
    /// `cutoff`. The result is the contact graph used by clash
    /// detection and contact maps.
    pub fn all_pairs_within(&self, cutoff: f64) -> Vec<(usize, usize, f64)> {
        let r = cutoff.min(self.cell_size);
        let r2 = r * r;
        let mut out = Vec::new();
        for (i, p) in self.points.iter().enumerate() {
            for j in self.within(p, r) {
                if j > i && distance_sq(p, &self.points[j]) <= r2 {
                    out.push((i, j, distance(p, &self.points[j])));
                }
            }
        }
        out
    }
}

/// Map a point to its integer cell coordinate.
fn cell_key(p: &Point3<f64>, origin: &Point3<f64>, size: f64) -> (i64, i64, i64) {
    (
        ((p.x - origin.x) / size).floor() as i64,
        ((p.y - origin.y) / size).floor() as i64,
        ((p.z - origin.z) / size).floor() as i64,
    )
}

/// Brute-force "all atoms within `radius` of `query`" — the reference
/// the grid is validated against; `O(n)` but allocation-light.
pub fn within_brute(points: &[Point3<f64>], query: &Point3<f64>, radius: f64) -> Vec<usize> {
    let r2 = radius * radius;
    points
        .iter()
        .enumerate()
        .filter(|(_, p)| distance_sq(query, p) <= r2)
        .map(|(i, _)| i)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lattice() -> Vec<Point3<f64>> {
        let mut v = Vec::new();
        for x in 0..5 {
            for y in 0..5 {
                for z in 0..5 {
                    v.push(Point3::new(x as f64, y as f64, z as f64));
                }
            }
        }
        v
    }

    #[test]
    fn pairwise_distance() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(3.0, 4.0, 0.0);
        assert!((distance(&a, &b) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn grid_matches_brute_force() {
        let pts = lattice();
        let grid = NeighborGrid::new(&pts, 1.5).unwrap();
        let q = Point3::new(2.0, 2.0, 2.0);
        let mut g = grid.within(&q, 1.5);
        let mut b = within_brute(&pts, &q, 1.5);
        g.sort_unstable();
        b.sort_unstable();
        assert_eq!(g, b);
    }

    #[test]
    fn grid_self_exclusion() {
        let pts = lattice();
        let grid = NeighborGrid::new(&pts, 1.2).unwrap();
        let n = grid.within_excluding(&pts[62], 1.2, 62);
        assert!(!n.contains(&62));
        // a point well inside the lattice has 6 axis neighbours at d=1.
        assert_eq!(n.len(), 6);
    }

    #[test]
    fn all_pairs_count() {
        // 3 collinear points 1 Å apart: pairs within 1.5 -> (0,1),(1,2)
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let grid = NeighborGrid::new(&pts, 1.5).unwrap();
        assert_eq!(grid.all_pairs_within(1.5).len(), 2);
    }

    #[test]
    fn rejects_bad_cutoff() {
        assert!(NeighborGrid::new(&[], 0.0).is_err());
        assert!(NeighborGrid::new(&[], -1.0).is_err());
    }

    #[test]
    fn empty_grid_is_empty() {
        let grid = NeighborGrid::new(&[], 1.0).unwrap();
        assert!(grid.is_empty());
        assert!(grid.within(&Point3::origin(), 1.0).is_empty());
    }
}
