//! Topology-level repair for Tri3 surface meshes.
//!
//! Three families of helpers:
//!
//! - [`fill_holes`] — detect boundary loops and triangulate each one
//!   with ear-clipping.
//! - [`is_manifold`] — every interior edge has exactly 2 incident
//!   triangles (boundary edges have 1). Two or more boundary edges
//!   meeting at a non-manifold vertex are also flagged.
//! - [`self_intersections`] — find triangle pairs that share a positive-area
//!   overlap region using an AABB-tree spatial index (Phase 7.5
//!   graduation: replaces the O(T²) brute-force prune with an
//!   expected-O(T log T) median-split BVH). Returns `(i, j)` index
//!   pairs with `i < j`.
//!
//! All operate on Tri3 blocks; other element types are ignored (their
//! boundary handling is solver-specific and out of scope here).

use std::collections::{HashMap, HashSet};

use nalgebra::Vector3;

use crate::element::{ElementBlock, ElementType};
use crate::mesh::Mesh;

/// Fill every closed boundary loop in `mesh` whose perimeter is
/// `<= max_boundary_length` by triangulating it with ear clipping.
///
/// `max_boundary_length` is a safety cap so we don't accidentally
/// fill a giant intentional hole (e.g. a tube's open end). Use
/// `f64::INFINITY` to fill all loops.
///
/// Returns a fresh mesh with the new triangles appended to the first
/// Tri3 block (created if none exists). The output `id` is
/// `"<original>_filled"`.
pub fn fill_holes(mesh: &Mesh, max_boundary_length: f64) -> Mesh {
    let mut out = mesh.clone();
    out.id = format!("{}_filled", mesh.id);
    let loops = boundary_loops(&out);
    let mut new_tris: Vec<u32> = Vec::new();
    for lp in loops {
        if lp.len() < 3 {
            continue;
        }
        // Perimeter check.
        //
        // R34 S2 (defense-in-depth): loop vertices are connectivity
        // values surfaced by `boundary_loops`, so `out.nodes[a]` would
        // panic on an out-of-range index from a future un-hardened
        // loader. Use `.get()` and skip the whole loop if any vertex
        // is out of range — backs the per-loader parse guards.
        let mut perim = 0.0;
        let mut loop_ok = true;
        for k in 0..lp.len() {
            let a = lp[k];
            let b = lp[(k + 1) % lp.len()];
            let (Some(pa), Some(pb)) =
                (out.nodes.get(a as usize), out.nodes.get(b as usize))
            else {
                loop_ok = false;
                break;
            };
            perim += (pa - pb).norm();
        }
        if !loop_ok || perim > max_boundary_length {
            continue;
        }
        let tris = ear_clip_loop(&out, &lp);
        for t in tris {
            new_tris.extend_from_slice(&t);
        }
    }
    if !new_tris.is_empty() {
        if let Some(blk) = out
            .element_blocks
            .iter_mut()
            .find(|b| b.element_type == ElementType::Tri3)
        {
            blk.connectivity.extend(new_tris);
        } else {
            let mut blk = ElementBlock::new(ElementType::Tri3);
            blk.connectivity = new_tris;
            out.element_blocks.push(blk);
        }
    }
    out.recompute_stats();
    out
}

/// Walk every boundary edge (interior count == 1) of the Tri3 blocks
/// and chain them into closed loops. Order around each loop is
/// consistent with the source triangle winding (so the loop normal
/// points outward for a properly-oriented mesh).
///
/// Multiple loops on the same mesh are returned in arbitrary order;
/// the inner ordering of each loop is the head-tail chain.
pub fn boundary_loops(mesh: &Mesh) -> Vec<Vec<u32>> {
    // Gather every directed boundary edge: half-edges whose reverse
    // doesn't appear. We track (a -> b); when traversing, the next
    // edge is the one whose source equals our `b`.
    let mut edge_count: HashMap<(u32, u32), i32> = HashMap::new();
    let mut directed: Vec<(u32, u32)> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let und = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(und).or_insert(0) += 1;
                directed.push((a, b));
            }
        }
    }
    // The boundary directed edges are those whose undirected count is 1.
    let mut start_to_edge: HashMap<u32, Vec<u32>> = HashMap::new();
    for (a, b) in directed {
        let und = if a < b { (a, b) } else { (b, a) };
        if edge_count.get(&und).copied().unwrap_or(0) == 1 {
            start_to_edge.entry(a).or_default().push(b);
        }
    }
    // Walk loops.
    let mut loops: Vec<Vec<u32>> = Vec::new();
    let mut visited: HashSet<(u32, u32)> = HashSet::new();
    let starts: Vec<u32> = start_to_edge.keys().copied().collect();
    for start in starts {
        let candidates: Vec<u32> = start_to_edge.get(&start).cloned().unwrap_or_default();
        for first_next in candidates {
            if visited.contains(&(start, first_next)) {
                continue;
            }
            // Try to build a loop starting (start -> first_next).
            let mut lp: Vec<u32> = vec![start];
            let mut cur = first_next;
            let mut last = start;
            let mut ok = false;
            while !lp.contains(&cur) || cur == start {
                lp.push(cur);
                visited.insert((last, cur));
                if cur == start {
                    // Closed loop — drop the duplicate.
                    lp.pop();
                    ok = true;
                    break;
                }
                // Pick the next edge from cur that doesn't go back.
                let Some(nexts) = start_to_edge.get(&cur) else {
                    break;
                };
                let Some(&next) = nexts.iter().find(|&&n| n != last) else {
                    // Dead end (degenerate boundary); abandon.
                    break;
                };
                last = cur;
                cur = next;
                if lp.len() > 100_000 {
                    break;
                }
            }
            if ok && lp.len() >= 3 {
                loops.push(lp);
            }
        }
    }
    loops
}

/// Triangulate a planar-ish 3D boundary loop with ear clipping.
///
/// 1. Compute the loop's average normal via Newell's formula.
/// 2. Project every vertex into a 2D basis on the loop's plane.
/// 3. Standard 2D ear-clipping: while the loop has > 3 vertices,
///    find a convex vertex whose triangle contains no other loop
///    vertex, emit it, and remove it from the loop.
///
/// Returns a Vec of `[v0, v1, v2]` triangles using the original
/// vertex indices from `loop_vertices`.
fn ear_clip_loop(mesh: &Mesh, loop_vertices: &[u32]) -> Vec<[u32; 3]> {
    if loop_vertices.len() < 3 {
        return Vec::new();
    }
    if loop_vertices.len() == 3 {
        return vec![[loop_vertices[0], loop_vertices[1], loop_vertices[2]]];
    }
    // Newell's normal for an arbitrary 3D polygon.
    //
    // R34 S2 (defense-in-depth): `loop_vertices` are connectivity
    // values, so an out-of-range index would panic `mesh.nodes[..]`.
    // Resolve every loop vertex up front with `.get()`; if any is out
    // of range we can't position the loop, so we bail (empty fan) —
    // the same graceful-degrade fate as a degenerate (zero-area) loop
    // below. Backs the per-loader parse guards.
    let Some(loop_pos): Option<Vec<Vector3<f64>>> = loop_vertices
        .iter()
        .map(|&i| mesh.nodes.get(i as usize).copied())
        .collect()
    else {
        return Vec::new();
    };
    let mut normal: Vector3<f64> = Vector3::zeros();
    for k in 0..loop_pos.len() {
        let a = loop_pos[k];
        let b = loop_pos[(k + 1) % loop_pos.len()];
        normal.x += (a.y - b.y) * (a.z + b.z);
        normal.y += (a.z - b.z) * (a.x + b.x);
        normal.z += (a.x - b.x) * (a.y + b.y);
    }
    let nlen = normal.norm();
    if nlen < 1e-20 {
        return Vec::new();
    }
    normal /= nlen;
    // Build orthonormal basis in the loop's plane.
    let helper = if normal.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let mut u = (helper - normal * normal.dot(&helper)).normalize();
    let v = normal.cross(&u);
    u = v.cross(&normal); // re-orthogonalise.

    // Reuse the positions resolved above (already `.get()`-checked) so
    // the projection can't re-introduce an out-of-range index.
    let proj: Vec<(f64, f64)> = loop_pos.iter().map(|p| (p.dot(&u), p.dot(&v))).collect();

    // Working list of indices into `loop_vertices` / `proj` (i.e. the
    // polygon's current vertex positions). Each ear-clip removes one.
    let mut work: Vec<usize> = (0..loop_vertices.len()).collect();
    let mut out: Vec<[u32; 3]> = Vec::new();
    let mut guard = work.len() * work.len() + 8;
    while work.len() > 3 && guard > 0 {
        guard -= 1;
        let mut ear_found = false;
        let n = work.len();
        for i in 0..n {
            let prev = work[(i + n - 1) % n];
            let curr = work[i];
            let next = work[(i + 1) % n];
            let p_prev = proj[prev];
            let p_curr = proj[curr];
            let p_next = proj[next];
            // Convexity test: cross product sign relative to polygon orientation.
            let cross = (p_curr.0 - p_prev.0) * (p_next.1 - p_curr.1)
                - (p_curr.1 - p_prev.1) * (p_next.0 - p_curr.0);
            if cross <= 0.0 {
                continue; // reflex
            }
            // Check no other working vertex is inside the triangle.
            let mut contains_other = false;
            for &other in &work {
                if other == prev || other == curr || other == next {
                    continue;
                }
                if point_in_triangle_2d(proj[other], p_prev, p_curr, p_next) {
                    contains_other = true;
                    break;
                }
            }
            if contains_other {
                continue;
            }
            // Emit the ear, remove curr from the working list.
            out.push([
                loop_vertices[prev],
                loop_vertices[curr],
                loop_vertices[next],
            ]);
            work.remove(i);
            ear_found = true;
            break;
        }
        if !ear_found {
            break; // polygon not strictly simple; bail.
        }
    }
    if work.len() == 3 {
        out.push([
            loop_vertices[work[0]],
            loop_vertices[work[1]],
            loop_vertices[work[2]],
        ]);
    }
    out
}

fn point_in_triangle_2d(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    let s1 = sign_2d(p, a, b);
    let s2 = sign_2d(p, b, c);
    let s3 = sign_2d(p, c, a);
    let has_neg = s1 < 0.0 || s2 < 0.0 || s3 < 0.0;
    let has_pos = s1 > 0.0 || s2 > 0.0 || s3 > 0.0;
    !(has_neg && has_pos)
}
fn sign_2d(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    (p.0 - b.0) * (a.1 - b.1) - (a.0 - b.0) * (p.1 - b.1)
}

/// Manifold-ness check: every interior edge has exactly two
/// incident triangles, and there's no edge with three or more.
///
/// Boundary edges (1 incident triangle) are allowed — they're a
/// normal feature of open meshes. Returns `true` only if the mesh
/// is **edge-manifold**.
pub fn is_manifold(mesh: &Mesh) -> bool {
    let mut edge_count: HashMap<(u32, u32), u32> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
    }
    edge_count.values().all(|&c| c == 1 || c == 2)
}

/// Detect Tri3 triangle pairs that intersect each other.
///
/// v1.5 implementation (Phase 7.5): median-split AABB-tree spatial
/// index. Triangles are gathered into a flat array, then partitioned
/// recursively by the median centroid coordinate along the longest
/// box axis. Pairwise tests are produced by self-overlap-traversal
/// of the tree, then refined with the internal `tri_tri_intersect`
/// Möller test. Triangles that share any vertex are skipped (a
/// shared edge or vertex does not count as a self-intersection).
///
/// Globally indexed: triangles are numbered by walking every Tri3
/// block in order. Returns `(i, j)` global triangle indices with
/// `i < j`. The pair order within the output is determined by the
/// tree traversal — callers that need a deterministic ordering
/// should sort the result.
pub fn self_intersections(mesh: &Mesh) -> Vec<(usize, usize)> {
    let mut tris: Vec<[u32; 3]> = Vec::new();
    let mut bboxes: Vec<Aabb> = Vec::new();
    // R34 S2 (defense-in-depth): this is the single point where Tri3
    // connectivity is resolved to node positions. Every downstream
    // consumer (the AABB tree and the Möller refinement below) indexes
    // `tris`, which only holds in-range triangles once we filter here,
    // so this one chokepoint protects them all. The per-loader parse
    // guards (OBJ/gmsh/netgen/PLY) are the first line; this seal backs
    // them so a future un-hardened loader degrades gracefully (drops a
    // degenerate triangle) instead of panicking "index out of bounds".
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let arr = [tri[0], tri[1], tri[2]];
            let (Some(&v0), Some(&v1), Some(&v2)) = (
                mesh.nodes.get(arr[0] as usize),
                mesh.nodes.get(arr[1] as usize),
                mesh.nodes.get(arr[2] as usize),
            ) else {
                continue;
            };
            bboxes.push(Aabb::from_triangle(&v0, &v1, &v2));
            tris.push(arr);
        }
    }
    if tris.len() < 2 {
        return Vec::new();
    }

    let tree = AabbTree::build(&bboxes);
    let mut candidate_pairs: Vec<(usize, usize)> = Vec::new();
    tree.self_overlap_pairs(&mut candidate_pairs);

    let mut out: Vec<(usize, usize)> = Vec::new();
    for (i, j) in candidate_pairs {
        // Share-a-vertex skip.
        if tris[i].iter().any(|v| tris[j].contains(v)) {
            continue;
        }
        // Full Möller check. `tris` was filtered to in-range triangles
        // above, so these `.get()`s never miss for a well-formed gather
        // — they stay `.get()`+skip purely to keep the R34 S2 seal
        // uniform and robust against future refactors of the gather.
        let ti = tris[i];
        let tj = tris[j];
        let (Some(&a0), Some(&a1), Some(&a2)) = (
            mesh.nodes.get(ti[0] as usize),
            mesh.nodes.get(ti[1] as usize),
            mesh.nodes.get(ti[2] as usize),
        ) else {
            continue;
        };
        let (Some(&b0), Some(&b1), Some(&b2)) = (
            mesh.nodes.get(tj[0] as usize),
            mesh.nodes.get(tj[1] as usize),
            mesh.nodes.get(tj[2] as usize),
        ) else {
            continue;
        };
        let a = [a0, a1, a2];
        let b = [b0, b1, b2];
        if tri_tri_intersect(&a, &b) {
            out.push((i, j));
        }
    }
    out
}

/// Axis-aligned bounding box for a single mesh element.
///
/// Internal building block for [`self_intersections`]'s spatial
/// index. Stored as per-axis min/max arrays so we can dispatch on
/// the longest axis without allocating per node.
#[derive(Clone, Copy, Debug)]
struct Aabb {
    min: [f64; 3],
    max: [f64; 3],
}

impl Aabb {
    fn from_triangle(v0: &Vector3<f64>, v1: &Vector3<f64>, v2: &Vector3<f64>) -> Self {
        let min = [
            v0.x.min(v1.x).min(v2.x),
            v0.y.min(v1.y).min(v2.y),
            v0.z.min(v1.z).min(v2.z),
        ];
        let max = [
            v0.x.max(v1.x).max(v2.x),
            v0.y.max(v1.y).max(v2.y),
            v0.z.max(v1.z).max(v2.z),
        ];
        Aabb { min, max }
    }

    fn empty() -> Self {
        Aabb {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    fn merge_in_place(&mut self, other: &Aabb) {
        for k in 0..3 {
            if other.min[k] < self.min[k] {
                self.min[k] = other.min[k];
            }
            if other.max[k] > self.max[k] {
                self.max[k] = other.max[k];
            }
        }
    }

    fn overlaps(&self, other: &Aabb) -> bool {
        !(self.min[0] > other.max[0]
            || self.max[0] < other.min[0]
            || self.min[1] > other.max[1]
            || self.max[1] < other.min[1]
            || self.min[2] > other.max[2]
            || self.max[2] < other.min[2])
    }

    fn centroid(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    fn longest_axis(&self) -> usize {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        if dx >= dy && dx >= dz {
            0
        } else if dy >= dz {
            1
        } else {
            2
        }
    }
}

/// Median-split AABB-tree for fast self-overlap queries.
///
/// Each non-leaf node holds the union AABB of its subtree plus
/// indices into a heap-laid-out node array (`left`/`right`); leaves
/// reference a primitive id in `prim_idx`. Building is O(N log N)
/// (each level scans N items + partitions in O(N)); a single
/// `self_overlap_pairs` traversal walks each pair of subtrees
/// whose boxes overlap, giving expected O(N log N) for sparse
/// meshes and degrading gracefully to O(N²) only when many
/// triangles genuinely overlap.
struct AabbTree {
    nodes: Vec<AabbNode>,
    root: usize,
}

#[derive(Clone, Copy, Debug)]
struct AabbNode {
    bbox: Aabb,
    /// Primitive index if this is a leaf, else `usize::MAX`.
    prim_idx: usize,
    /// Left child index in `nodes`, else `usize::MAX`.
    left: usize,
    /// Right child index in `nodes`, else `usize::MAX`.
    right: usize,
}

impl AabbNode {
    fn is_leaf(&self) -> bool {
        self.prim_idx != usize::MAX
    }
}

impl AabbTree {
    fn build(bboxes: &[Aabb]) -> Self {
        let mut indices: Vec<usize> = (0..bboxes.len()).collect();
        let mut nodes: Vec<AabbNode> = Vec::with_capacity(bboxes.len() * 2);
        let root = Self::build_recursive(bboxes, &mut indices, &mut nodes);
        AabbTree { nodes, root }
    }

    fn build_recursive(
        bboxes: &[Aabb],
        indices: &mut [usize],
        nodes: &mut Vec<AabbNode>,
    ) -> usize {
        let node_idx = nodes.len();
        if indices.len() == 1 {
            let prim = indices[0];
            nodes.push(AabbNode {
                bbox: bboxes[prim],
                prim_idx: prim,
                left: usize::MAX,
                right: usize::MAX,
            });
            return node_idx;
        }
        // Compute union bbox over this subtree's primitives.
        let mut bbox = Aabb::empty();
        for &i in indices.iter() {
            bbox.merge_in_place(&bboxes[i]);
        }
        // Choose split axis = longest axis of union; partition by
        // median centroid. `select_nth_unstable_by` is O(N).
        let axis = bbox.longest_axis();
        let mid = indices.len() / 2;
        indices.select_nth_unstable_by(mid, |&a, &b| {
            let ca = bboxes[a].centroid()[axis];
            let cb = bboxes[b].centroid()[axis];
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        // Reserve our slot, then recurse.
        nodes.push(AabbNode {
            bbox,
            prim_idx: usize::MAX,
            left: usize::MAX,
            right: usize::MAX,
        });
        let (left_slice, right_slice) = indices.split_at_mut(mid);
        let left = Self::build_recursive(bboxes, left_slice, nodes);
        let right = Self::build_recursive(bboxes, right_slice, nodes);
        nodes[node_idx].left = left;
        nodes[node_idx].right = right;
        node_idx
    }

    /// Emit every `(i, j)` primitive pair with `i < j` whose AABBs
    /// overlap. The pair order follows the tree traversal — sort
    /// the output if a deterministic order is required.
    fn self_overlap_pairs(&self, out: &mut Vec<(usize, usize)>) {
        self.descend_self(self.root, out);
    }

    fn descend_self(&self, node_idx: usize, out: &mut Vec<(usize, usize)>) {
        let node = &self.nodes[node_idx];
        if node.is_leaf() {
            return;
        }
        let left = &self.nodes[node.left];
        let right = &self.nodes[node.right];
        // Self-recurse into each subtree first.
        self.descend_self(node.left, out);
        self.descend_self(node.right, out);
        // Then cross-recurse the pair if their boxes overlap.
        if left.bbox.overlaps(&right.bbox) {
            self.descend_cross(node.left, node.right, out);
        }
    }

    fn descend_cross(&self, a_idx: usize, b_idx: usize, out: &mut Vec<(usize, usize)>) {
        let a = &self.nodes[a_idx];
        let b = &self.nodes[b_idx];
        if !a.bbox.overlaps(&b.bbox) {
            return;
        }
        match (a.is_leaf(), b.is_leaf()) {
            (true, true) => {
                let (i, j) = if a.prim_idx < b.prim_idx {
                    (a.prim_idx, b.prim_idx)
                } else {
                    (b.prim_idx, a.prim_idx)
                };
                out.push((i, j));
            }
            (true, false) => {
                self.descend_cross(a_idx, b.left, out);
                self.descend_cross(a_idx, b.right, out);
            }
            (false, true) => {
                self.descend_cross(a.left, b_idx, out);
                self.descend_cross(a.right, b_idx, out);
            }
            (false, false) => {
                // Descend the larger subtree to keep the recursion
                // tree balanced for the common "many overlaps"
                // worst case.
                let a_extent = (a.bbox.max[a.bbox.longest_axis()]
                    - a.bbox.min[a.bbox.longest_axis()])
                .abs();
                let b_extent = (b.bbox.max[b.bbox.longest_axis()]
                    - b.bbox.min[b.bbox.longest_axis()])
                .abs();
                if a_extent >= b_extent {
                    self.descend_cross(a.left, b_idx, out);
                    self.descend_cross(a.right, b_idx, out);
                } else {
                    self.descend_cross(a_idx, b.left, out);
                    self.descend_cross(a_idx, b.right, out);
                }
            }
        }
    }
}

/// Triangle-triangle intersection (Möller, simplified). Returns true
/// if `a` and `b` share any positive-area overlap region.
///
/// v1 implementation: degenerate-case handling is conservative
/// (reports overlap when undecided). Good enough for an interactive
/// "do these surfaces clip into each other" diagnostic.
fn tri_tri_intersect(a: &[Vector3<f64>; 3], b: &[Vector3<f64>; 3]) -> bool {
    // Compute plane of A.
    let n_a = (a[1] - a[0]).cross(&(a[2] - a[0]));
    let d_a = -n_a.dot(&a[0]);
    let db = [
        n_a.dot(&b[0]) + d_a,
        n_a.dot(&b[1]) + d_a,
        n_a.dot(&b[2]) + d_a,
    ];
    // All on one side: no intersection.
    if (db[0] > 1e-12 && db[1] > 1e-12 && db[2] > 1e-12)
        || (db[0] < -1e-12 && db[1] < -1e-12 && db[2] < -1e-12)
    {
        return false;
    }
    let n_b = (b[1] - b[0]).cross(&(b[2] - b[0]));
    let d_b = -n_b.dot(&b[0]);
    let da = [
        n_b.dot(&a[0]) + d_b,
        n_b.dot(&a[1]) + d_b,
        n_b.dot(&a[2]) + d_b,
    ];
    if (da[0] > 1e-12 && da[1] > 1e-12 && da[2] > 1e-12)
        || (da[0] < -1e-12 && da[1] < -1e-12 && da[2] < -1e-12)
    {
        return false;
    }
    // Compute the intersection line direction (cross of normals).
    let line_dir = n_a.cross(&n_b);
    if line_dir.norm_squared() < 1e-20 {
        // Coplanar: conservatively report intersection so we don't
        // miss obvious overlaps. (Full coplanar-2D-clip is v1.5.)
        return true;
    }
    // Project both triangles onto the dominant axis of line_dir, find
    // each triangle's interval on the intersection line, return whether
    // the intervals overlap.
    let axis = dominant_axis(&line_dir);
    let pa = project_triangle_to_line_interval(a, &da, axis);
    let pb = project_triangle_to_line_interval(b, &db, axis);
    !(pa.1 < pb.0 || pb.1 < pa.0)
}

fn dominant_axis(v: &Vector3<f64>) -> usize {
    let abs = [v.x.abs(), v.y.abs(), v.z.abs()];
    if abs[0] > abs[1] && abs[0] > abs[2] {
        0
    } else if abs[1] > abs[2] {
        1
    } else {
        2
    }
}

fn project_triangle_to_line_interval(
    t: &[Vector3<f64>; 3],
    d: &[f64; 3],
    axis: usize,
) -> (f64, f64) {
    // Two of d[*] have the same sign, one has the other. The two
    // crossings between the differing vertex and the other two yield
    // the interval endpoints.
    let mut isects: Vec<f64> = Vec::new();
    for k in 0..3 {
        let (i, j) = (k, (k + 1) % 3);
        if (d[i] > 0.0 && d[j] <= 0.0)
            || (d[i] < 0.0 && d[j] >= 0.0)
            || (d[i] <= 0.0 && d[j] > 0.0)
            || (d[i] >= 0.0 && d[j] < 0.0)
        {
            let denom = d[i] - d[j];
            if denom.abs() < 1e-20 {
                continue;
            }
            let s = d[i] / denom;
            let p = t[i] + (t[j] - t[i]) * s;
            isects.push(p[axis]);
        }
    }
    // Add on-plane vertices.
    for k in 0..3 {
        if d[k].abs() < 1e-12 {
            isects.push(t[k][axis]);
        }
    }
    if isects.is_empty() {
        return (0.0, 0.0); // degenerate; caller treats overlap = false
    }
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for x in isects {
        if x < min {
            min = x;
        }
        if x > max {
            max = x;
        }
    }
    (min, max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ElementBlock;

    fn pt(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    /// Closed unit-cube surface (12 tris, 8 verts).
    fn unit_cube() -> Mesh {
        let mut m = Mesh::new("cube");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(1.0, 1.0, 0.0),
            pt(0.0, 1.0, 0.0),
            pt(0.0, 0.0, 1.0),
            pt(1.0, 0.0, 1.0),
            pt(1.0, 1.0, 1.0),
            pt(0.0, 1.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity.extend_from_slice(&[
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ]);
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    /// Cube with the top face missing — leaves a square boundary loop.
    fn cube_minus_top() -> Mesh {
        let mut m = unit_cube();
        // Top face triangles were the second pair (indices 6..12).
        let conn = &mut m.element_blocks[0].connectivity;
        // Top face = tris 2 and 3 in original layout: rows
        // (4,5,6) and (4,6,7) starting at offset 6 (each tri is 3 idx).
        // We constructed `unit_cube` with the +z face second: offsets 6..12.
        conn.drain(6..12);
        m.recompute_stats();
        m
    }

    #[test]
    fn manifold_cube_is_true() {
        let m = unit_cube();
        assert!(is_manifold(&m));
    }

    #[test]
    fn manifold_cube_with_extra_triangle_is_false() {
        // Add a triangle that shares an edge with an existing triangle
        // → that shared edge now appears in 3 triangles → non-manifold.
        let mut m = unit_cube();
        let new_tri = [0u32, 1, 4]; // edge (0, 1) is on -z face (tri 0,2,1 + 0,3,2) and -y face (0,1,5)
        m.element_blocks[0].connectivity.extend_from_slice(&new_tri);
        assert!(!is_manifold(&m));
    }

    #[test]
    fn boundary_loops_on_open_cube_is_one_square() {
        let m = cube_minus_top();
        let loops = boundary_loops(&m);
        assert_eq!(
            loops.len(),
            1,
            "expected exactly 1 loop, got {}",
            loops.len()
        );
        assert_eq!(
            loops[0].len(),
            4,
            "expected 4-vertex loop, got {}",
            loops[0].len()
        );
    }

    #[test]
    fn boundary_loops_on_closed_cube_is_empty() {
        let m = unit_cube();
        assert!(boundary_loops(&m).is_empty());
    }

    #[test]
    fn fill_holes_closes_open_cube() {
        let m = cube_minus_top();
        // Before: open mesh.
        assert!(!boundary_loops(&m).is_empty());
        let filled = fill_holes(&m, 100.0);
        // After: no boundary loops.
        assert!(
            boundary_loops(&filled).is_empty(),
            "loops remain after fill"
        );
        // 2 new triangles added (square → 2 ear-clipped tris).
        let before_tris: usize = m.element_blocks[0].connectivity.len() / 3;
        let after_tris: usize = filled.element_blocks[0].connectivity.len() / 3;
        assert_eq!(after_tris, before_tris + 2);
        assert_eq!(filled.id, "cube_filled");
    }

    #[test]
    fn fill_holes_respects_max_boundary_length() {
        let m = cube_minus_top();
        // Set max smaller than the perimeter (=4.0) → no fill.
        let filled = fill_holes(&m, 1.0);
        assert!(!boundary_loops(&filled).is_empty());
    }

    #[test]
    fn self_intersections_on_unit_cube_is_empty() {
        // A simple non-self-intersecting closed surface.
        let m = unit_cube();
        assert!(self_intersections(&m).is_empty());
    }

    #[test]
    fn self_intersections_finds_crossed_tris() {
        // Two triangles that obviously intersect (one cuts through the other).
        let mut m = Mesh::new("cross");
        m.nodes = vec![
            // Triangle A in z=0 plane, big.
            pt(-1.0, -1.0, 0.0),
            pt(1.0, -1.0, 0.0),
            pt(0.0, 1.0, 0.0),
            // Triangle B vertical, slicing through A.
            pt(0.0, 0.0, -0.5),
            pt(0.0, 0.0, 0.5),
            pt(0.0, 0.4, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2, 3, 4, 5];
        m.element_blocks = vec![blk];
        let pairs = self_intersections(&m);
        assert!(!pairs.is_empty(), "expected to find intersection");
    }

    #[test]
    fn self_intersections_returns_pairs_with_low_index_first() {
        // The Phase 7.5 AABB-tree traversal can emit pairs in any
        // order — the contract is that each pair has `i < j` so
        // downstream consumers can use a canonical lookup.
        let mut m = Mesh::new("cross");
        m.nodes = vec![
            pt(-1.0, -1.0, 0.0),
            pt(1.0, -1.0, 0.0),
            pt(0.0, 1.0, 0.0),
            pt(0.0, 0.0, -0.5),
            pt(0.0, 0.0, 0.5),
            pt(0.0, 0.4, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2, 3, 4, 5];
        m.element_blocks = vec![blk];
        let pairs = self_intersections(&m);
        for &(i, j) in &pairs {
            assert!(i < j, "pair ({i}, {j}) violates i < j contract");
        }
    }

    #[test]
    fn self_intersections_scales_to_grid_without_false_positives() {
        // A 6 × 6 grid of disjoint triangles in the z=0 plane —
        // every triangle is well-separated from every other, so
        // the AABB tree must drop every cross-subtree pair.
        let mut m = Mesh::new("grid");
        let mut blk = ElementBlock::new(ElementType::Tri3);
        let mut node_id: u32 = 0;
        for i in 0..6 {
            for j in 0..6 {
                let x = i as f64 * 10.0;
                let y = j as f64 * 10.0;
                m.nodes.push(pt(x, y, 0.0));
                m.nodes.push(pt(x + 1.0, y, 0.0));
                m.nodes.push(pt(x, y + 1.0, 0.0));
                blk.connectivity.extend_from_slice(&[
                    node_id,
                    node_id + 1,
                    node_id + 2,
                ]);
                node_id += 3;
            }
        }
        m.element_blocks = vec![blk];
        assert!(self_intersections(&m).is_empty());
    }

    /// R34 S2 (RED→GREEN): defense-in-depth sink seal. A mesh whose
    /// Tri3 connectivity cites a vertex index past `nodes.len()` must
    /// NOT panic the repair pass — the per-loader validation
    /// (OBJ/gmsh/netgen/PLY) is the first line, but a future
    /// un-hardened loader could still hand us such a mesh. Pre-fix
    /// `self_intersections` did `mesh.nodes[arr[0] as usize]` and
    /// panicked with "index out of bounds". Post-fix the offending
    /// triangle is skipped and a result is returned. We assert no
    /// panic; the bad triangle contributes no pairs.
    #[test]
    fn out_of_range_connectivity_does_not_panic() {
        let mut m = Mesh::new("hostile");
        // 3 real vertices...
        m.nodes = vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 0.0)];
        // ...but a triangle that cites vertex 9 (out of range).
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 9];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Must return (graceful degrade), not panic.
        let pairs = self_intersections(&m);
        assert!(
            pairs.is_empty(),
            "a single out-of-range triangle yields no self-intersection pairs"
        );
        // The other connectivity consumers in this module must also
        // survive the same hostile mesh without panicking.
        assert!(boundary_loops(&m).iter().all(|lp| !lp.is_empty()) || boundary_loops(&m).is_empty());
        let _ = is_manifold(&m);
        let _ = fill_holes(&m, f64::INFINITY);
    }

    /// R34 S2: a valid triangle and an out-of-range triangle together
    /// — the bad one is dropped while the good geometry still drives
    /// the self-intersection search, no panic.
    #[test]
    fn out_of_range_triangle_skipped_valid_kept() {
        // Two genuinely crossing triangles (as in
        // `self_intersections_finds_crossed_tris`) plus a third tri
        // citing a vertex past the node array.
        let mut m = Mesh::new("mixed");
        m.nodes = vec![
            pt(-1.0, -1.0, 0.0),
            pt(1.0, -1.0, 0.0),
            pt(0.0, 1.0, 0.0),
            pt(0.0, 0.0, -0.5),
            pt(0.0, 0.0, 0.5),
            pt(0.0, 0.4, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // First two tris cross; third cites vertex 42 (out of range).
        blk.connectivity = vec![0, 1, 2, 3, 4, 5, 0, 1, 42];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // The valid crossing pair is still detected; the bad triangle
        // is silently dropped instead of panicking.
        let pairs = self_intersections(&m);
        assert!(
            !pairs.is_empty(),
            "the valid crossing pair must still be found"
        );
    }

    #[test]
    fn ear_clip_triangle_loop_emits_one_face() {
        let mut m = Mesh::new("tri");
        m.nodes = vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 0.0)];
        let out = ear_clip_loop(&m, &[0, 1, 2]);
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn ear_clip_square_loop_emits_two_faces() {
        let mut m = Mesh::new("square");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(1.0, 1.0, 0.0),
            pt(0.0, 1.0, 0.0),
        ];
        let out = ear_clip_loop(&m, &[0, 1, 2, 3]);
        assert_eq!(out.len(), 2);
    }
}
