//! O(log N) spatial queries over geometry, backed by an R*-tree.
//!
//! Linear scans over point sets are the usual bottleneck once a model
//! grows past a few thousand vertices: a nearest-neighbour or "what is
//! inside this box" lookup is O(N) per query, so an N-vs-N pass is
//! O(N²). This module wraps the [`rstar`](https://docs.rs/rstar) R*-tree
//! so the same queries run in roughly O(log N).
//!
//! The unit of indexing is [`SpatialPoint`]: a 3-D position plus an
//! arbitrary payload `T` (typically a vertex index, a node id, or a
//! handle into some other structure). Queries hand back references to
//! the stored points, so the caller recovers their own payload — the
//! index never has to own the heavyweight geometry.
//!
//! ```
//! use valenx_geo::spatial_index::SpatialIndex;
//! use nalgebra::Vector3;
//!
//! // Index four labelled points.
//! let idx = SpatialIndex::from_points([
//!     (Vector3::new(0.0, 0.0, 0.0), "origin"),
//!     (Vector3::new(10.0, 0.0, 0.0), "east"),
//!     (Vector3::new(0.0, 10.0, 0.0), "north"),
//!     (Vector3::new(10.0, 10.0, 0.0), "ne"),
//! ]);
//!
//! // Nearest neighbour to a query position.
//! let near = idx.nearest(Vector3::new(9.0, 1.0, 0.0)).unwrap();
//! assert_eq!(near.data, "east");
//! ```

use crate::bounding_box::BoundingBox;
use nalgebra::Vector3;
use rstar::{PointDistance, RTree, RTreeObject, AABB};

/// rstar's primitive point type for this index: a 3-D `f64` coordinate.
type Point3 = [f64; 3];

/// A point stored in a [`SpatialIndex`]: a position plus a payload.
///
/// `T` rides along with the coordinate so query results carry whatever
/// the caller needs (a vertex index, an id, a label, …) without the
/// index owning the original geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialPoint<T> {
    /// World-space position of the point.
    pub pos: Point3,
    /// Caller-supplied payload travelling with the point.
    pub data: T,
}

impl<T> SpatialPoint<T> {
    /// Construct from an [`nalgebra`] vector and a payload.
    pub fn new(pos: Vector3<f64>, data: T) -> Self {
        Self {
            pos: [pos.x, pos.y, pos.z],
            data,
        }
    }

    /// Position as an [`nalgebra`] vector.
    pub fn position(&self) -> Vector3<f64> {
        Vector3::new(self.pos[0], self.pos[1], self.pos[2])
    }
}

impl<T> RTreeObject for SpatialPoint<T> {
    type Envelope = AABB<Point3>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(self.pos)
    }
}

impl<T> PointDistance for SpatialPoint<T> {
    /// Squared distance from this point to `query`. Required by the
    /// nearest-neighbour and `locate_within_distance` queries.
    fn distance_2(&self, query: &Point3) -> f64 {
        dist_sq(&self.pos, query)
    }
}

/// Squared Euclidean distance between two 3-D points.
///
/// Squared distances avoid a `sqrt` and preserve ordering, which is all
/// the nearest-neighbour comparisons and radius checks need.
#[inline]
fn dist_sq(a: &Point3, b: &Point3) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// An R*-tree spatial index over [`SpatialPoint`]s.
///
/// Build it once with [`SpatialIndex::from_points`] (bulk loading packs
/// the tree for good query performance) and then run as many
/// [`nearest`](Self::nearest) / [`k_nearest`](Self::k_nearest) /
/// [`points_in_box`](Self::points_in_box) queries as you like.
pub struct SpatialIndex<T> {
    tree: RTree<SpatialPoint<T>>,
}

impl<T> SpatialIndex<T>
where
    T: Clone,
{
    /// Build an index from `(position, payload)` pairs.
    ///
    /// Uses rstar's bulk loader, which is both faster to construct and
    /// produces a better-balanced tree than repeated insertion.
    pub fn from_points<I>(points: I) -> Self
    where
        I: IntoIterator<Item = (Vector3<f64>, T)>,
    {
        let pts: Vec<SpatialPoint<T>> = points
            .into_iter()
            .map(|(pos, data)| SpatialPoint::new(pos, data))
            .collect();
        Self {
            tree: RTree::bulk_load(pts),
        }
    }

    /// Number of points in the index.
    pub fn len(&self) -> usize {
        self.tree.size()
    }

    /// `true` if the index holds no points.
    pub fn is_empty(&self) -> bool {
        self.tree.size() == 0
    }

    /// Nearest point to `query`, or `None` if the index is empty.
    ///
    /// Runs in ~O(log N).
    pub fn nearest(&self, query: Vector3<f64>) -> Option<&SpatialPoint<T>> {
        let q: Point3 = [query.x, query.y, query.z];
        self.tree.nearest_neighbor(q)
    }

    /// The `k` nearest points to `query`, ordered nearest-first.
    ///
    /// Returns fewer than `k` only when the index holds fewer than `k`
    /// points. `k == 0` yields an empty vector.
    pub fn k_nearest(&self, query: Vector3<f64>, k: usize) -> Vec<&SpatialPoint<T>> {
        if k == 0 {
            return Vec::new();
        }
        let q: Point3 = [query.x, query.y, query.z];
        // `nearest_neighbor_iter` yields points in increasing distance;
        // the tree walk is lazy, so taking `k` stays ~O(k log N).
        self.tree.nearest_neighbor_iter(q).take(k).collect()
    }

    /// Every point whose position lies within `bbox` (inclusive of the
    /// faces), in arbitrary order.
    ///
    /// This is the range / window query: it touches only the tree nodes
    /// overlapping `bbox`, so it is ~O(log N + m) for `m` hits rather
    /// than scanning all N points.
    pub fn points_in_box(&self, bbox: &BoundingBox) -> Vec<&SpatialPoint<T>> {
        let envelope = AABB::from_corners(
            [bbox.min.x, bbox.min.y, bbox.min.z],
            [bbox.max.x, bbox.max.y, bbox.max.z],
        );
        self.tree.locate_in_envelope(envelope).collect()
    }

    /// Every point within Euclidean `radius` of `query`, in arbitrary
    /// order.
    ///
    /// A spherical range query, complementing the axis-aligned
    /// [`points_in_box`](Self::points_in_box).
    pub fn points_within_radius(&self, query: Vector3<f64>, radius: f64) -> Vec<&SpatialPoint<T>> {
        let q: Point3 = [query.x, query.y, query.z];
        let r2 = radius * radius;
        self.tree
            .locate_within_distance(q, r2)
            .filter(|p| dist_sq(&p.pos, &q) <= r2)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Five points on a line at x = 0, 1, 2, 3, 4 (y = z = 0), each
    /// labelled by its index.
    fn line_index() -> SpatialIndex<usize> {
        SpatialIndex::from_points((0..5).map(|i| (Vector3::new(i as f64, 0.0, 0.0), i)))
    }

    #[test]
    fn nearest_returns_closest_point() {
        let idx = line_index();
        // 2.4 is closest to the point at x = 2.
        let n = idx.nearest(Vector3::new(2.4, 0.0, 0.0)).unwrap();
        assert_eq!(n.data, 2);
        // 3.6 is closest to the point at x = 4.
        let n = idx.nearest(Vector3::new(3.6, 0.0, 0.0)).unwrap();
        assert_eq!(n.data, 4);
    }

    #[test]
    fn nearest_on_empty_is_none() {
        let idx = SpatialIndex::<usize>::from_points(std::iter::empty());
        assert!(idx.is_empty());
        assert_eq!(idx.len(), 0);
        assert!(idx.nearest(Vector3::new(0.0, 0.0, 0.0)).is_none());
    }

    #[test]
    fn nearest_matches_brute_force_3d() {
        // A scattered 3-D cloud; the R*-tree result must equal the
        // linear-scan nearest neighbour for every query.
        let pts: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(-4.0, 0.5, 7.0),
            Vector3::new(2.0, -2.0, -1.0),
            Vector3::new(9.0, 9.0, 9.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(-3.0, -3.0, 2.0),
        ];
        let idx = SpatialIndex::from_points(pts.iter().enumerate().map(|(i, p)| (*p, i)));

        let queries = [
            Vector3::new(1.1, 1.9, 3.2),
            Vector3::new(-3.5, -3.0, 2.1),
            Vector3::new(8.0, 8.0, 8.0),
            Vector3::new(0.1, 0.1, -0.1),
        ];
        for q in queries {
            let brute = pts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    (*a - q).norm_squared().total_cmp(&(*b - q).norm_squared())
                })
                .map(|(i, _)| i)
                .unwrap();
            let got = idx.nearest(q).unwrap().data;
            assert_eq!(got, brute, "query {q:?}");
        }
    }

    #[test]
    fn k_nearest_returns_k_in_order() {
        let idx = line_index();
        // Around x = 1.1 the order by distance is 1, 0/2 (tie ~), then 3, 4.
        let got: Vec<usize> = idx
            .k_nearest(Vector3::new(1.1, 0.0, 0.0), 3)
            .iter()
            .map(|p| p.data)
            .collect();
        assert_eq!(got.len(), 3);
        // Nearest is the point at x = 1.
        assert_eq!(got[0], 1);
        // The next two are x = 0 and x = 2 (distances 1.1 and 0.9).
        let mut rest = [got[1], got[2]];
        rest.sort_unstable();
        assert_eq!(rest, [0, 2]);

        // Strict distance ordering check on an asymmetric query.
        let ordered: Vec<usize> = idx
            .k_nearest(Vector3::new(0.2, 0.0, 0.0), 5)
            .iter()
            .map(|p| p.data)
            .collect();
        assert_eq!(ordered, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn k_nearest_clamps_and_zero() {
        let idx = line_index();
        // Asking for more than exist returns them all.
        assert_eq!(idx.k_nearest(Vector3::new(0.0, 0.0, 0.0), 99).len(), 5);
        // k = 0 is empty.
        assert!(idx.k_nearest(Vector3::new(0.0, 0.0, 0.0), 0).is_empty());
    }

    #[test]
    fn bbox_query_returns_exactly_inside_points() {
        // 3x3x1 grid of points at integer x,y in [0,2].
        let mut pts = Vec::new();
        let mut id = 0;
        for x in 0..3 {
            for y in 0..3 {
                pts.push((Vector3::new(x as f64, y as f64, 0.0), id));
                id += 1;
            }
        }
        let idx = SpatialIndex::from_points(pts);

        // Box covering x,y in [0.5, 2.5] -> the 2x2 block at x,y in {1,2}:
        // ids for (1,1),(1,2),(2,1),(2,2).
        let bbox = BoundingBox::new(Vector3::new(0.5, 0.5, -0.5), Vector3::new(2.5, 2.5, 0.5));
        let mut got: Vec<i32> = idx.points_in_box(&bbox).iter().map(|p| p.data).collect();
        got.sort_unstable();

        // (x,y)->id mapping is id = x*3 + y.
        // reason: keep the `x * 3 + y` literal form so the mapping is self-documenting.
        #[allow(clippy::identity_op)]
        let expected = vec![1 * 3 + 1, 1 * 3 + 2, 2 * 3 + 1, 2 * 3 + 2]; // 4,5,7,8
        assert_eq!(got, expected);
        assert_eq!(got.len(), 4);
    }

    #[test]
    fn bbox_query_includes_boundary() {
        let idx = line_index();
        // Box [1,3] on x must include x = 1, 2, 3 (faces inclusive).
        let bbox = BoundingBox::new(Vector3::new(1.0, -0.1, -0.1), Vector3::new(3.0, 0.1, 0.1));
        let mut got: Vec<usize> = idx.points_in_box(&bbox).iter().map(|p| p.data).collect();
        got.sort_unstable();
        assert_eq!(got, vec![1, 2, 3]);
    }

    #[test]
    fn radius_query_returns_inside_sphere() {
        let idx = line_index();
        // Radius 1.5 around x = 2 -> points at x = 1, 2, 3 (dist 1,0,1);
        // x = 0 and 4 are at distance 2 > 1.5, excluded.
        let mut got: Vec<usize> = idx
            .points_within_radius(Vector3::new(2.0, 0.0, 0.0), 1.5)
            .iter()
            .map(|p| p.data)
            .collect();
        got.sort_unstable();
        assert_eq!(got, vec![1, 2, 3]);
    }
}
