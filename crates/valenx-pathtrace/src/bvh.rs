//! Bounding-volume hierarchy — the spatial acceleration structure that
//! turns the renderer's `O(triangles)` brute-force ray cast into an
//! `O(log triangles)` traversal.
//!
//! # What it is
//!
//! A binary tree of axis-aligned bounding boxes. Each leaf holds a
//! small run of triangles; each interior node bounds its two children.
//! A ray descends the tree, skipping any subtree whose box it misses,
//! so a ray that hits nothing touches a handful of boxes instead of
//! every triangle.
//!
//! # How it is built
//!
//! Top-down recursive subdivision with a **binned surface-area
//! heuristic** (SAH). At each node the triangles' centroids are binned
//! along the longest axis; for every bin boundary the SAH cost
//! `area(L)·count(L) + area(R)·count(R)` is evaluated and the cheapest
//! split is taken. The SAH is the standard production split metric — it
//! adapts to the geometry's density far better than a plain
//! spatial-median cut. If no binned split beats leaving the node whole
//! (or the node is already small), a leaf is emitted.
//!
//! The tree is stored as a **flat `Vec<BvhNode>`** — index 0 is the
//! root, an interior node stores the index of its right child (the left
//! child is always the next slot), and a leaf stores a `(start, count)`
//! range into a reordered triangle-index array. A flat array keeps the
//! traversal cache-friendly and the whole structure trivially
//! `Clone`able with no pointer chasing.

use crate::geometry::{Aabb, Hit, Ray, Triangle};

/// One node of the flattened BVH.
///
/// A leaf and an interior node share the struct; `tri_count > 0`
/// distinguishes a leaf. For a leaf, `payload` is the start index into
/// [`Bvh::tri_indices`]; for an interior node, `payload` is the array
/// index of the **right** child (the left child is `self_index + 1`).
#[derive(Clone, Copy, Debug)]
struct BvhNode {
    /// The node's bounding box.
    bounds: Aabb,
    /// Leaf: start offset into `tri_indices`. Interior: right-child
    /// node index.
    payload: u32,
    /// Number of triangles in this leaf, or 0 for an interior node.
    tri_count: u32,
}

/// A built bounding-volume hierarchy over a triangle list.
///
/// Construct one with [`Bvh::build`]; query it with
/// [`Bvh::intersect`] (nearest hit) and [`Bvh::occluded`] (any-hit, for
/// shadow rays).
#[derive(Clone, Debug)]
pub struct Bvh {
    /// The flattened node array; `nodes[0]` is the root.
    nodes: Vec<BvhNode>,
    /// Triangle indices, reordered so each leaf's triangles form a
    /// contiguous run. `intersect` indexes the caller's triangle slice
    /// through this.
    tri_indices: Vec<u32>,
}

/// Triangles below this count in a node become a leaf without further
/// subdivision — a small linear scan beats the traversal overhead.
const MAX_LEAF_TRIS: usize = 4;

/// Number of SAH bins evaluated per split axis. 12 is the usual
/// production sweet spot — enough resolution to find a good split,
/// cheap enough that the build stays fast.
const SAH_BINS: usize = 12;

/// Working record for one triangle during the build — its bounds and
/// centroid, cached so the recursion never recomputes them.
#[derive(Clone, Copy)]
struct BuildPrim {
    bounds: Aabb,
    centroid: crate::math::Vec3,
    index: u32,
}

impl Bvh {
    /// Build a BVH over `triangles`.
    ///
    /// Returns an empty hierarchy (a single empty leaf) for an empty
    /// triangle list, so [`Bvh::intersect`] is always safe to call.
    pub fn build(triangles: &[Triangle]) -> Bvh {
        if triangles.is_empty() {
            return Bvh {
                nodes: vec![BvhNode {
                    bounds: Aabb::empty(),
                    payload: 0,
                    tri_count: 0,
                }],
                tri_indices: Vec::new(),
            };
        }
        // Cache each triangle's bounds + centroid up front.
        let mut prims: Vec<BuildPrim> = triangles
            .iter()
            .enumerate()
            .map(|(i, t)| {
                let bounds = t.bounds();
                BuildPrim {
                    bounds,
                    centroid: t.centroid(),
                    index: i as u32,
                }
            })
            .collect();

        let mut nodes: Vec<BvhNode> = Vec::with_capacity(2 * triangles.len());
        let mut tri_indices: Vec<u32> = Vec::with_capacity(triangles.len());
        // Reserve the root slot; subdivide fills it.
        nodes.push(BvhNode {
            bounds: Aabb::empty(),
            payload: 0,
            tri_count: 0,
        });
        Self::subdivide(&mut prims, 0, &mut nodes, &mut tri_indices, 0);

        Bvh { nodes, tri_indices }
    }

    /// Recursively subdivide `prims[..]` into the node at `node_index`.
    ///
    /// `prims` is partitioned in place; `nodes` and `tri_indices` grow
    /// as the tree is filled. `depth` guards against pathological deep
    /// recursion on degenerate input.
    fn subdivide(
        prims: &mut [BuildPrim],
        node_index: usize,
        nodes: &mut Vec<BvhNode>,
        tri_indices: &mut Vec<u32>,
        depth: u32,
    ) {
        // Bounds of every triangle in this node.
        let mut bounds = Aabb::empty();
        let mut centroid_bounds = Aabb::empty();
        for p in prims.iter() {
            bounds.expand_box(&p.bounds);
            centroid_bounds.expand_point(p.centroid);
        }
        nodes[node_index].bounds = bounds;

        // Leaf cut-offs: too few triangles, too deep, or every centroid
        // collapsed to a point (no axis to split on).
        let make_leaf = |prims: &[BuildPrim],
                         node_index: usize,
                         nodes: &mut Vec<BvhNode>,
                         tri_indices: &mut Vec<u32>| {
            let start = tri_indices.len() as u32;
            for p in prims {
                tri_indices.push(p.index);
            }
            nodes[node_index].payload = start;
            nodes[node_index].tri_count = prims.len() as u32;
        };

        let extent = centroid_bounds.extent();
        let axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        let axis_extent = extent.axis(axis);
        if prims.len() <= MAX_LEAF_TRIS || depth >= 48 || axis_extent < 1e-12 {
            make_leaf(prims, node_index, nodes, tri_indices);
            return;
        }

        // --- binned SAH ---
        // Bucket each triangle by its centroid along the split axis.
        #[derive(Clone, Copy)]
        struct Bin {
            count: u32,
            bounds: Aabb,
        }
        let mut bins = [Bin {
            count: 0,
            bounds: Aabb::empty(),
        }; SAH_BINS];
        let axis_min = centroid_bounds.min.axis(axis);
        let scale = SAH_BINS as f32 / axis_extent;
        for p in prims.iter() {
            let mut b = ((p.centroid.axis(axis) - axis_min) * scale) as usize;
            b = b.min(SAH_BINS - 1);
            bins[b].count += 1;
            bins[b].bounds.expand_box(&p.bounds);
        }

        // Sweep the SAH_BINS-1 candidate split planes. For each plane
        // accumulate the left running cost and the right running cost.
        let mut left_area = [0.0f32; SAH_BINS - 1];
        let mut left_count = [0u32; SAH_BINS - 1];
        {
            let mut box_acc = Aabb::empty();
            let mut count_acc = 0u32;
            for i in 0..SAH_BINS - 1 {
                count_acc += bins[i].count;
                box_acc.expand_box(&bins[i].bounds);
                left_count[i] = count_acc;
                left_area[i] = box_acc.surface_area();
            }
        }
        let mut right_area = [0.0f32; SAH_BINS - 1];
        let mut right_count = [0u32; SAH_BINS - 1];
        {
            let mut box_acc = Aabb::empty();
            let mut count_acc = 0u32;
            for i in (0..SAH_BINS - 1).rev() {
                count_acc += bins[i + 1].count;
                box_acc.expand_box(&bins[i + 1].bounds);
                right_count[i] = count_acc;
                right_area[i] = box_acc.surface_area();
            }
        }

        // Pick the cheapest split: cost = area_L·count_L + area_R·count_R.
        let mut best_cost = f32::INFINITY;
        let mut best_split = usize::MAX;
        for i in 0..SAH_BINS - 1 {
            if left_count[i] == 0 || right_count[i] == 0 {
                continue;
            }
            let cost = left_area[i] * left_count[i] as f32 + right_area[i] * right_count[i] as f32;
            if cost < best_cost {
                best_cost = cost;
                best_split = i;
            }
        }

        // The cost of *not* splitting (leaving this node a leaf).
        let leaf_cost = bounds.surface_area() * prims.len() as f32;
        if best_split == usize::MAX || (best_cost >= leaf_cost && prims.len() <= 16) {
            // No worthwhile split — emit a leaf. (The `<= 16` guard
            // still forces a split for big nodes even if the SAH ties,
            // so a degenerate-but-large node never becomes a slow leaf.)
            make_leaf(prims, node_index, nodes, tri_indices);
            return;
        }

        // Partition `prims` in place: everything in a bin ≤ best_split
        // goes left, the rest right.
        let split_bin = best_split + 1;
        let mid = itertools_partition(prims, |p| {
            let mut b = ((p.centroid.axis(axis) - axis_min) * scale) as usize;
            b = b.min(SAH_BINS - 1);
            b < split_bin
        });
        // A pathological case: every centroid landed on one side of the
        // chosen plane. Fall back to a median split so we still make
        // progress.
        let mid = if mid == 0 || mid == prims.len() {
            prims.len() / 2
        } else {
            mid
        };

        // The flat-array invariant the traversal relies on is "left
        // child = node_index + 1, and the whole left subtree is laid
        // out before the right child". So the left child must be
        // allocated *and fully recursed* before the right child slot is
        // taken — otherwise the left subtree's nodes would land between
        // the parent and its right child. (Allocating both children
        // up-front and only then recursing breaks the invariant for
        // every node below the root.)
        let (left_prims, right_prims) = prims.split_at_mut(mid);

        // Left child: the slot immediately after this node.
        let left_index = nodes.len();
        debug_assert_eq!(left_index, node_index + 1);
        nodes.push(BvhNode {
            bounds: Aabb::empty(),
            payload: 0,
            tri_count: 0,
        });
        Self::subdivide(left_prims, left_index, nodes, tri_indices, depth + 1);

        // Right child: the slot after the entire left subtree.
        let right_index = nodes.len();
        nodes.push(BvhNode {
            bounds: Aabb::empty(),
            payload: 0,
            tri_count: 0,
        });
        nodes[node_index].payload = right_index as u32;
        nodes[node_index].tri_count = 0; // interior
        Self::subdivide(right_prims, right_index, nodes, tri_indices, depth + 1);
    }

    /// The number of nodes in the tree (interior + leaf). Exposed for
    /// tests and diagnostics.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Find the nearest triangle the ray strikes in `(t_min, t_max)`.
    ///
    /// `triangles` must be the exact slice the BVH was [`Bvh::build`]'d
    /// over — the hierarchy stores indices into it. Returns the closest
    /// [`Hit`], or `None` if the ray escapes the scene.
    ///
    /// Traversal is an explicit stack (no recursion): the box of each
    /// node is slab-tested against the ray's *current* `[t_min,
    /// t_max]`; interior nodes push both children, leaves linearly test
    /// their triangle run, and every hit tightens `t_max` so later
    /// subtrees beyond the closest hit are skipped.
    pub fn intersect(
        &self,
        triangles: &[Triangle],
        ray: &Ray,
        t_min: f32,
        t_max: f32,
    ) -> Option<Hit> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut closest: Option<Hit> = None;
        let mut t_far = t_max;
        // A 64-deep stack covers any tree the 48-depth build cap can
        // produce with comfortable headroom.
        let mut stack: [u32; 64] = [0; 64];
        let mut sp: usize = 0;
        stack[sp] = 0;
        sp += 1;

        while sp > 0 {
            sp -= 1;
            let this = stack[sp];
            let node = self.nodes[this as usize];
            if !node.bounds.hit(ray, t_min, t_far) {
                continue;
            }
            if node.tri_count > 0 {
                // Leaf — test its triangle run.
                let start = node.payload as usize;
                for &ti in &self.tri_indices[start..start + node.tri_count as usize] {
                    if let Some(hit) = triangles[ti as usize].intersect(ray, t_min, t_far) {
                        t_far = hit.t;
                        closest = Some(hit);
                    }
                }
            } else {
                // Interior — push both children. The left child is the
                // next array slot; the right child is `payload`.
                let left_child = this + 1;
                let right_child = node.payload;
                if sp + 2 <= stack.len() {
                    stack[sp] = right_child;
                    sp += 1;
                    stack[sp] = left_child;
                    sp += 1;
                }
            }
        }
        closest
    }

    /// Test whether *any* triangle blocks the ray within `(t_min,
    /// t_max)` — the shadow-ray query.
    ///
    /// Returns `true` on the first hit found (no need to find the
    /// nearest), so it is cheaper than [`Bvh::intersect`]. Used to test
    /// occlusion between a surface point and a light sample.
    pub fn occluded(&self, triangles: &[Triangle], ray: &Ray, t_min: f32, t_max: f32) -> bool {
        if self.nodes.is_empty() {
            return false;
        }
        let mut stack: [u32; 64] = [0; 64];
        let mut sp: usize = 0;
        stack[sp] = 0;
        sp += 1;
        while sp > 0 {
            sp -= 1;
            let idx = stack[sp];
            let node = self.nodes[idx as usize];
            if !node.bounds.hit(ray, t_min, t_max) {
                continue;
            }
            if node.tri_count > 0 {
                let start = node.payload as usize;
                for &ti in &self.tri_indices[start..start + node.tri_count as usize] {
                    if triangles[ti as usize]
                        .intersect(ray, t_min, t_max)
                        .is_some()
                    {
                        return true;
                    }
                }
            } else {
                let left_child = idx + 1;
                let right_child = node.payload;
                if sp + 2 <= stack.len() {
                    stack[sp] = right_child;
                    sp += 1;
                    stack[sp] = left_child;
                    sp += 1;
                }
            }
        }
        false
    }
}

/// In-place partition: move every element for which `pred` is true to
/// the front, returning the index of the first false element. A
/// hand-rolled `Vec::partition`-style two-pointer sweep (the std
/// `[T]::iter_mut().partition` allocates; this does not).
fn itertools_partition<T, F: Fn(&T) -> bool>(slice: &mut [T], pred: F) -> usize {
    let mut i = 0;
    let mut j = slice.len();
    loop {
        while i < j && pred(&slice[i]) {
            i += 1;
        }
        while i < j && !pred(&slice[j - 1]) {
            j -= 1;
        }
        if i >= j {
            return i;
        }
        slice.swap(i, j - 1);
        i += 1;
        j -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    /// A grid of `n × n` unit triangles in the z = 0 plane.
    fn triangle_grid(n: usize) -> Vec<Triangle> {
        let mut tris = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let x = i as f32;
                let y = j as f32;
                tris.push(Triangle::flat(
                    [
                        vec3(x, y, 0.0),
                        vec3(x + 0.9, y, 0.0),
                        vec3(x, y + 0.9, 0.0),
                    ],
                    0,
                ));
            }
        }
        tris
    }

    #[test]
    fn empty_bvh_intersects_nothing() {
        let bvh = Bvh::build(&[]);
        let ray = Ray::new(vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, -1.0));
        assert!(bvh.intersect(&[], &ray, 1e-4, 1e30).is_none());
        assert!(!bvh.occluded(&[], &ray, 1e-4, 1e30));
    }

    #[test]
    fn bvh_finds_the_same_hit_as_brute_force() {
        // For a randomised-ish set of rays the BVH must agree with a
        // linear scan — this is the correctness contract.
        let tris = triangle_grid(8);
        let bvh = Bvh::build(&tris);
        // March a few rays straight down at known triangle interiors.
        for j in 0..8 {
            for i in 0..8 {
                let ox = i as f32 + 0.2;
                let oy = j as f32 + 0.2;
                let ray = Ray::new(vec3(ox, oy, 5.0), vec3(0.0, 0.0, -1.0));
                let bvh_hit = bvh.intersect(&tris, &ray, 1e-4, 1e30);
                // Brute force.
                let mut brute: Option<Hit> = None;
                let mut t_far = 1e30;
                for t in &tris {
                    if let Some(h) = t.intersect(&ray, 1e-4, t_far) {
                        t_far = h.t;
                        brute = Some(h);
                    }
                }
                match (bvh_hit, brute) {
                    (Some(a), Some(b)) => {
                        assert!((a.t - b.t).abs() < 1e-4, "BVH t {} vs brute t {}", a.t, b.t);
                    }
                    (None, None) => {}
                    _ => panic!("BVH / brute-force disagree at ({i}, {j})"),
                }
            }
        }
    }

    #[test]
    fn bvh_nearest_hit_is_the_closest_triangle() {
        // Two parallel triangles at different depths; the ray must
        // return the nearer one.
        let near = Triangle::flat(
            [
                vec3(-1.0, -1.0, 1.0),
                vec3(1.0, -1.0, 1.0),
                vec3(0.0, 1.0, 1.0),
            ],
            0,
        );
        let far = Triangle::flat(
            [
                vec3(-1.0, -1.0, -1.0),
                vec3(1.0, -1.0, -1.0),
                vec3(0.0, 1.0, -1.0),
            ],
            1,
        );
        let tris = vec![far, near]; // deliberately far-first
        let bvh = Bvh::build(&tris);
        let ray = Ray::new(vec3(0.0, 0.0, 5.0), vec3(0.0, 0.0, -1.0));
        let hit = bvh.intersect(&tris, &ray, 1e-4, 1e30).unwrap();
        // The near triangle is at z = 1 → t ≈ 4; material index 0.
        assert!((hit.t - 4.0).abs() < 1e-4, "should hit the near triangle");
        assert_eq!(hit.material, 0, "near triangle has material 0");
    }

    #[test]
    fn occluded_detects_a_blocker_and_clears_past_it() {
        let tris = triangle_grid(4);
        let bvh = Bvh::build(&tris);
        // A ray through a triangle, with t_max past it → occluded.
        let blocked = Ray::new(vec3(0.2, 0.2, 5.0), vec3(0.0, 0.0, -1.0));
        assert!(bvh.occluded(&tris, &blocked, 1e-4, 1e30));
        // The same ray but t_max stops short of the triangle → clear.
        assert!(!bvh.occluded(&tris, &blocked, 1e-4, 1.0));
        // A ray through empty space between triangles → clear.
        let gap = Ray::new(vec3(0.95, 0.95, 5.0), vec3(0.0, 0.0, -1.0));
        assert!(!bvh.occluded(&tris, &gap, 1e-4, 1e30));
    }

    #[test]
    fn bvh_over_many_triangles_builds_a_real_tree() {
        // A 16×16 grid is 256 triangles; the SAH build must produce a
        // multi-node hierarchy, not one giant leaf.
        let tris = triangle_grid(16);
        let bvh = Bvh::build(&tris);
        assert!(
            bvh.node_count() > 1,
            "256 triangles should yield an interior tree, got {} nodes",
            bvh.node_count()
        );
    }

    #[test]
    fn partition_moves_matching_elements_to_the_front() {
        let mut v = vec![1, 4, 2, 5, 3, 6];
        let mid = itertools_partition(&mut v, |&x| x <= 3);
        assert_eq!(mid, 3);
        // Front three are all ≤ 3, back three all > 3 (order within
        // each side is unspecified).
        assert!(v[..3].iter().all(|&x| x <= 3));
        assert!(v[3..].iter().all(|&x| x > 3));
    }
}
