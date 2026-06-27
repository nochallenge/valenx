//! GPU-render-oriented mesh optimization passes (vertex-cache,
//! overdraw, vertex-fetch reordering) and LOD simplification.
//!
//! These wrap the C++ [`meshoptimizer`] library (via the `meshopt`
//! crate) over the canonical [`Mesh`]. They target the **render**
//! representation of a surface mesh — the order triangle indices are
//! submitted to the GPU and the order vertices are laid out in
//! memory — and a fast appearance-preserving LOD decimator.
//!
//! [`meshoptimizer`]: https://github.com/zeux/meshoptimizer
//!
//! ## Relationship to [`crate::decimate`]
//!
//! [`crate::decimate::quadric_error_decimate`] is valenx's in-house
//! Garland–Heckbert QEM decimator — it is the right tool when you
//! want a precise, vertex-budget-driven reduction with valenx-owned
//! semantics. The [`simplify`] here is a *different* algorithm
//! (meshopt's error-bounded edge collapse, tuned for real-time LOD
//! chains) and is complemented by the cache / overdraw / fetch
//! reordering passes, which the QEM decimator does not provide.
//! Both are exposed; callers pick by use case.
//!
//! ## Scope
//!
//! All passes operate on **`Tri3`** surface blocks. Triangles whose
//! indices reference a vertex past `nodes.len()` are dropped
//! (defensive — matches [`crate::decimate`]). Non-`Tri3` element
//! blocks are preserved with their connectivity remapped through the
//! survivor table when vertices are reordered/compacted, and left
//! untouched when only indices are reordered. `regions` and
//! `boundaries` index into the original element arrays and are not
//! reconstructed across a reorder/decimation, so they are dropped on
//! the returned mesh (same policy as the QEM decimator).

// NB: `meshopt::optimize_vertex_cache` is intentionally NOT imported
// unqualified — this module exposes its own `optimize_vertex_cache`
// over `&Mesh`, so the meshoptimizer free function is always called
// fully-qualified to avoid the name clash.
use meshopt::{
    optimize_overdraw_in_place, optimize_vertex_fetch, simplify as mo_simplify, SimplifyOptions,
    VertexDataAdapter,
};
use nalgebra::Vector3;

use crate::element::{ElementBlock, ElementType};
use crate::mesh::Mesh;

/// Default overdraw threshold passed to [`optimize`].
///
/// `1.05` lets the overdraw optimizer degrade vertex-cache
/// efficiency by up to 5 % in exchange for reduced pixel overdraw —
/// the value recommended by meshoptimizer for the common case.
pub const DEFAULT_OVERDRAW_THRESHOLD: f32 = 1.05;

/// Collect every valid `Tri3` triangle as a flat index buffer, plus
/// the non-`Tri3` blocks (cloned) to splice back into the output.
///
/// Triangles citing an out-of-range vertex are skipped (defensive
/// seal — see module docs).
fn collect_tris(mesh: &Mesh) -> (Vec<u32>, Vec<ElementBlock>) {
    let node_count = mesh.nodes.len() as u32;
    let mut indices: Vec<u32> = Vec::new();
    let mut other_blocks: Vec<ElementBlock> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type == ElementType::Tri3 {
            for t in block.connectivity.chunks_exact(3) {
                if t[0] >= node_count || t[1] >= node_count || t[2] >= node_count {
                    continue;
                }
                indices.extend_from_slice(&[t[0], t[1], t[2]]);
            }
        } else {
            other_blocks.push(block.clone());
        }
    }
    (indices, other_blocks)
}

/// Pack node coordinates into the tightly-interleaved `f32` position
/// **byte** buffer meshoptimizer expects (3 × little-/native-endian
/// `f32` per vertex, stride 12, position offset 0). meshoptimizer is
/// single-precision internally, so the canonical `f64` nodes are
/// narrowed for the optimizer's quadric / cache analysis only — the
/// returned mesh keeps the original `f64` node coordinates untouched
/// (we only ever reorder or drop whole vertices, never edit their
/// values, except where [`simplify`]'s collapses move a vertex, which
/// it does in `f32`).
///
/// We materialize a `Vec<u8>` via `f32::to_ne_bytes` rather than
/// reinterpreting an `&[f32]` so the crate's `#![forbid(unsafe_code)]`
/// holds with zero `unsafe`. The extra allocation is negligible next
/// to the C++ optimizer's own work. `to_ne_bytes` matches the host
/// endianness `meshopt_*` reads the buffer with.
fn pack_position_bytes(nodes: &[Vector3<f64>]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(nodes.len() * 12);
    for n in nodes {
        buf.extend_from_slice(&(n.x as f32).to_ne_bytes());
        buf.extend_from_slice(&(n.y as f32).to_ne_bytes());
        buf.extend_from_slice(&(n.z as f32).to_ne_bytes());
    }
    buf
}

/// Build a `VertexDataAdapter` over a packed position byte buffer
/// (stride 12, offset 0).
///
/// Returns `None` for an empty buffer so callers can early-out (the
/// adapter requires `data.len() % stride == 0`, which always holds
/// for our fixed 12-byte stride, but an empty mesh has no positions
/// to analyze).
fn adapter(position_bytes: &[u8]) -> Option<VertexDataAdapter<'_>> {
    if position_bytes.is_empty() {
        return None;
    }
    VertexDataAdapter::new(position_bytes, 12, 0).ok()
}

/// Rebuild a fresh `Mesh` from optimized node + Tri3-index data,
/// splicing the (already index-remapped) non-`Tri3` blocks back in.
fn rebuild(
    src_id: &str,
    suffix: &str,
    nodes: Vec<Vector3<f64>>,
    tri_indices: Vec<u32>,
    other_blocks: Vec<ElementBlock>,
) -> Mesh {
    let mut out = Mesh::new(format!("{src_id}_{suffix}"));
    out.nodes = nodes;
    if !tri_indices.is_empty() {
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = tri_indices;
        out.element_blocks.push(block);
    }
    out.element_blocks.extend(other_blocks);
    out.recompute_stats();
    out
}

/// Reorder the `Tri3` index buffer to improve **GPU vertex-cache**
/// efficiency (meshoptimizer's `optimizeVertexCache`).
///
/// This is a pure permutation of the triangle index list: the
/// returned mesh has the **same vertices** (unchanged positions,
/// unchanged count, unchanged order) and the **same set of
/// triangles** — only the order indices are submitted changes, so
/// post-transform-cache GPU vertex-shader invocations drop. Geometry
/// and topology are bit-for-bit identical.
///
/// Non-`Tri3` blocks pass through untouched. The output `id` is
/// `"<original>_vcache"`.
pub fn optimize_vertex_cache(mesh: &Mesh) -> Mesh {
    let (indices, other_blocks) = collect_tris(mesh);
    if indices.is_empty() {
        return rebuild(
            &mesh.id,
            "vcache",
            mesh.nodes.clone(),
            Vec::new(),
            other_blocks,
        );
    }
    let reordered = meshopt::optimize_vertex_cache(&indices, mesh.nodes.len());
    rebuild(
        &mesh.id,
        "vcache",
        mesh.nodes.clone(),
        reordered,
        other_blocks,
    )
}

/// Full render-optimization pipeline over the `Tri3` surface:
/// vertex-cache → overdraw → vertex-fetch.
///
/// 1. **vertex-cache** reorders indices for the post-transform cache.
/// 2. **overdraw** reorders indices further to reduce pixel overdraw,
///    allowed to spend up to `overdraw_threshold` (e.g. `1.05` = 5 %)
///    of cache efficiency.
/// 3. **vertex-fetch** reorders *vertices* to match the index access
///    order and **compacts** them, dropping any vertex no triangle
///    references. Indices are rewritten to the new vertex order.
///
/// Because step 3 reorders and may drop vertices, the non-`Tri3`
/// blocks' connectivity is remapped through the same vertex
/// permutation; any non-`Tri3` element that referenced a vertex which
/// the fetch pass dropped (i.e. used by no surviving triangle) is
/// itself dropped. If you need volume blocks preserved verbatim, use
/// [`optimize_vertex_cache`] alone (which never touches vertices).
///
/// The output `id` is `"<original>_optimized"`.
pub fn optimize(mesh: &Mesh, overdraw_threshold: f32) -> Mesh {
    let (indices, other_blocks) = collect_tris(mesh);
    if indices.is_empty() {
        return rebuild(
            &mesh.id,
            "optimized",
            mesh.nodes.clone(),
            Vec::new(),
            other_blocks,
        );
    }

    let position_bytes = pack_position_bytes(&mesh.nodes);
    let Some(adapter) = adapter(&position_bytes) else {
        // No usable position buffer — fall back to a cache-only
        // reorder, vertices unchanged.
        let reordered = meshopt::optimize_vertex_cache(&indices, mesh.nodes.len());
        return rebuild(
            &mesh.id,
            "optimized",
            mesh.nodes.clone(),
            reordered,
            other_blocks,
        );
    };

    // 1. vertex cache.
    let mut idx = meshopt::optimize_vertex_cache(&indices, mesh.nodes.len());
    // 2. overdraw (requires cache-optimized input — satisfied above).
    optimize_overdraw_in_place(&mut idx, &adapter, overdraw_threshold);
    // 3. vertex fetch — reorder + compact the *original f64* nodes so
    //    they stay double precision. `optimize_vertex_fetch` rewrites
    //    `idx` in place to the compacted vertex order and returns the
    //    new node list.
    //
    // Build the old->new vertex remap so non-Tri3 blocks follow the
    // SAME permutation. `optimize_vertex_fetch_remap` and
    // `optimize_vertex_fetch` agree on the new vertex numbering only
    // when fed identical index buffers, so the remap is computed from
    // the post-cache/overdraw `idx` snapshot taken *before* the
    // in-place fetch rewrites it. `remap[old] = new` (`u32::MAX` for a
    // vertex no surviving triangle references).
    let remap = meshopt::optimize_vertex_fetch_remap(&idx, mesh.nodes.len());
    let new_nodes = optimize_vertex_fetch(&mut idx, &mesh.nodes);
    let remapped_other = remap_blocks(&other_blocks, &remap);

    rebuild(&mesh.id, "optimized", new_nodes, idx, remapped_other)
}

/// Apply a `remap[old] = new` table (with `u32::MAX` marking dropped
/// vertices) to non-`Tri3` blocks. Elements that reference a dropped
/// vertex are removed; survivors have their indices rewritten.
fn remap_blocks(blocks: &[ElementBlock], remap: &[u32]) -> Vec<ElementBlock> {
    let mut out = Vec::with_capacity(blocks.len());
    for block in blocks {
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let mut conn: Vec<u32> = Vec::with_capacity(block.connectivity.len());
        'elem: for elem in block.connectivity.chunks_exact(npe) {
            let mut mapped = [0u32; 32];
            for (k, &v) in elem.iter().enumerate() {
                let new = remap.get(v as usize).copied().unwrap_or(u32::MAX);
                if new == u32::MAX {
                    // References a dropped vertex — drop this element.
                    continue 'elem;
                }
                mapped[k] = new;
            }
            conn.extend_from_slice(&mapped[..npe]);
        }
        if !conn.is_empty() {
            out.push(ElementBlock {
                element_type: block.element_type,
                connectivity: conn,
            });
        }
    }
    out
}

/// **LOD simplification** — reduce the `Tri3` triangle count toward
/// `target_ratio × current` while preserving appearance, using
/// meshoptimizer's error-bounded edge-collapse simplifier.
///
/// `target_ratio` is clamped to `[0.0, 1.0]`: `0.5` aims for half the
/// triangles, `1.0` is (effectively) a no-op, `0.0` collapses as far
/// as the error bound allows. The simplifier respects an internal
/// relative error target (`1e-2` of the mesh extent) and will stop
/// short of the triangle target rather than introduce gross error, so
/// the result is *bounded-error* rather than an exact count.
///
/// The returned mesh's vertex buffer is **compacted** to only the
/// vertices the simplified triangles reference (via vertex-fetch), so
/// `nodes.len()` shrinks too. Non-`Tri3` blocks are remapped through
/// that same compaction and dropped if they referenced a removed
/// vertex (same policy as [`optimize`]). The output `id` is
/// `"<original>_lod"`.
pub fn simplify(mesh: &Mesh, target_ratio: f64) -> Mesh {
    let ratio = target_ratio.clamp(0.0, 1.0);
    let (indices, other_blocks) = collect_tris(mesh);
    if indices.is_empty() {
        return rebuild(
            &mesh.id,
            "lod",
            mesh.nodes.clone(),
            Vec::new(),
            other_blocks,
        );
    }

    let position_bytes = pack_position_bytes(&mesh.nodes);
    let Some(adapter) = adapter(&position_bytes) else {
        return rebuild(&mesh.id, "lod", mesh.nodes.clone(), indices, other_blocks);
    };

    // Target *index* count = triangles * 3. Keep at least one triangle.
    let target_index_count = {
        let tris = indices.len() / 3;
        let want = ((tris as f64) * ratio).round() as usize;
        want.max(1) * 3
    };

    // Relative error budget: allow collapses up to 1% of the mesh
    // extent. The simplifier returns early (fewer collapses) if it
    // cannot reach `target_index_count` within this budget — i.e.
    // error stays bounded even when the count target is not met.
    let target_error = 1e-2_f32;
    let mut result_error = 0.0_f32;
    let simplified = mo_simplify(
        &indices,
        &adapter,
        target_index_count,
        target_error,
        SimplifyOptions::None,
        Some(&mut result_error),
    );

    // Compact vertices to the simplified index set so the LOD mesh is
    // self-contained and small. The remap (`old -> new`, `u32::MAX`
    // for dropped) MUST be derived from the *same* index buffer the
    // fetch pass consumes — they share meshoptimizer's first-use
    // ordering only when computed over identical indices — so we
    // snapshot `simplified` before `optimize_vertex_fetch` rewrites it
    // in place, then remap the non-`Tri3` blocks through that table.
    let mut idx = simplified;
    let remap = meshopt::optimize_vertex_fetch_remap(&idx, mesh.nodes.len());
    let new_nodes = optimize_vertex_fetch(&mut idx, &mesh.nodes);
    let remapped_other = remap_blocks(&other_blocks, &remap);

    rebuild(&mesh.id, "lod", new_nodes, idx, remapped_other)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::{ElementBlock, ElementType};

    /// Build an `n × n`-cell subdivided unit plane in the XY plane as
    /// a `Tri3` mesh. `(n+1)²` vertices, `2 n²` triangles. A dense,
    /// flat mesh — ideal for exercising reduction (lots of coplanar
    /// triangles the simplifier can safely collapse) and for checking
    /// that reordering preserves geometry.
    fn subdivided_plane(n: usize) -> Mesh {
        let mut m = Mesh::new("plane");
        let step = 1.0 / n as f64;
        for j in 0..=n {
            for i in 0..=n {
                m.nodes
                    .push(Vector3::new(i as f64 * step, j as f64 * step, 0.0));
            }
        }
        let idx = |i: usize, j: usize| (j * (n + 1) + i) as u32;
        let mut block = ElementBlock::new(ElementType::Tri3);
        for j in 0..n {
            for i in 0..n {
                let a = idx(i, j);
                let b = idx(i + 1, j);
                let c = idx(i + 1, j + 1);
                let d = idx(i, j + 1);
                block.connectivity.extend_from_slice(&[a, b, c]);
                block.connectivity.extend_from_slice(&[a, c, d]);
            }
        }
        m.element_blocks.push(block);
        m.recompute_stats();
        m
    }

    /// Closed unit-cube surface (8 verts, 12 Tri3 tris) — a watertight
    /// manifold for sanity checks.
    fn unit_cube_surface() -> Mesh {
        let mut m = Mesh::new("cube");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        #[rustfmt::skip]
        let tris: [u32; 36] = [
            0, 2, 1, 0, 3, 2, // bottom
            4, 5, 6, 4, 6, 7, // top
            0, 1, 5, 0, 5, 4, // front
            1, 2, 6, 1, 6, 5, // right
            2, 3, 7, 2, 7, 6, // back
            3, 0, 4, 3, 4, 7, // left
        ];
        block.connectivity = tris.to_vec();
        m.element_blocks.push(block);
        m.recompute_stats();
        m
    }

    fn tri_count(m: &Mesh) -> usize {
        m.element_blocks
            .iter()
            .filter(|b| b.element_type == ElementType::Tri3)
            .map(|b| b.connectivity.len() / 3)
            .sum()
    }

    fn aabb(m: &Mesh) -> (Vector3<f64>, Vector3<f64>) {
        let mut min = m.nodes[0];
        let mut max = m.nodes[0];
        for n in &m.nodes {
            for k in 0..3 {
                if n[k] < min[k] {
                    min[k] = n[k];
                }
                if n[k] > max[k] {
                    max[k] = n[k];
                }
            }
        }
        (min, max)
    }

    #[test]
    fn empty_mesh_passthrough() {
        let m = Mesh::new("empty");
        for out in [
            optimize_vertex_cache(&m),
            optimize(&m, DEFAULT_OVERDRAW_THRESHOLD),
            simplify(&m, 0.5),
        ] {
            assert!(out.nodes.is_empty());
            assert_eq!(tri_count(&out), 0);
        }
    }

    #[test]
    fn vertex_cache_preserves_geometry_and_triangle_set() {
        // The cache pass is a pure index permutation: same vertices
        // (count, order, values) and the same *set* of triangles
        // (as sorted-vertex tuples), only index ORDER differs.
        let m = unit_cube_surface();
        let out = optimize_vertex_cache(&m);
        assert_eq!(out.id, "cube_vcache");

        // Vertices identical, in the same order.
        assert_eq!(out.nodes.len(), m.nodes.len());
        for (a, b) in out.nodes.iter().zip(m.nodes.iter()) {
            assert_eq!(a, b);
        }

        // Same number of triangles.
        assert_eq!(tri_count(&out), tri_count(&m));

        // Same multiset of triangles (each as a sorted index triple).
        let canon = |mesh: &Mesh| {
            let mut v: Vec<[u32; 3]> = mesh.element_blocks[0]
                .connectivity
                .chunks_exact(3)
                .map(|t| {
                    let mut s = [t[0], t[1], t[2]];
                    s.sort_unstable();
                    s
                })
                .collect();
            v.sort_unstable();
            v
        };
        assert_eq!(canon(&out), canon(&m));
    }

    #[test]
    fn vertex_cache_passes_through_non_tri3() {
        // A Tet4-only mesh has no Tri3 surface; nodes + the volume
        // block pass through unchanged.
        let mut m = Mesh::new("tet");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tet4);
        blk.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![blk];
        let out = optimize_vertex_cache(&m);
        assert_eq!(out.nodes.len(), 4);
        assert_eq!(out.element_blocks.len(), 1);
        assert_eq!(out.element_blocks[0].element_type, ElementType::Tet4);
        assert_eq!(out.element_blocks[0].connectivity, vec![0, 1, 2, 3]);
    }

    #[test]
    fn optimize_compacts_and_preserves_geometry() {
        // Full pipeline on the cube. Every cube vertex is used by a
        // triangle, so vertex-fetch keeps all 8 (reordered) and the
        // triangle count is unchanged. AABB must not drift (pure
        // reorder, no collapses).
        let m = unit_cube_surface();
        let (min0, max0) = aabb(&m);
        let out = optimize(&m, DEFAULT_OVERDRAW_THRESHOLD);
        assert_eq!(out.id, "cube_optimized");
        assert_eq!(out.nodes.len(), m.nodes.len());
        assert_eq!(tri_count(&out), tri_count(&m));
        let (min1, max1) = aabb(&out);
        for k in 0..3 {
            assert!((min0[k] - min1[k]).abs() < 1e-9);
            assert!((max0[k] - max1[k]).abs() < 1e-9);
        }
    }

    #[test]
    fn optimize_drops_unreferenced_vertices() {
        // Add an orphan vertex no triangle uses; vertex-fetch should
        // compact it away.
        let mut m = unit_cube_surface();
        m.nodes.push(Vector3::new(5.0, 5.0, 5.0)); // orphan
        m.recompute_stats();
        let out = optimize(&m, DEFAULT_OVERDRAW_THRESHOLD);
        // 8 referenced verts survive, the orphan is gone.
        assert_eq!(out.nodes.len(), 8);
        assert_eq!(tri_count(&out), 12);
    }

    #[test]
    fn simplify_reduces_triangle_count_toward_target() {
        // Dense 20x20 plane: 800 triangles, 441 vertices, all
        // coplanar -> the simplifier can collapse aggressively within
        // the error budget. Target 50% should noticeably reduce.
        let m = subdivided_plane(20);
        let before = tri_count(&m);
        assert_eq!(before, 800);
        let out = simplify(&m, 0.5);
        assert_eq!(out.id, "plane_lod");
        let after = tri_count(&out);
        assert!(
            after < before,
            "simplify must reduce triangles: before={before} after={after}"
        );
        // A flat plane is the easy case; it should reach at most the
        // requested count (never exceed it).
        assert!(
            after <= before / 2,
            "flat plane should reach the 50% target: before={before} after={after}"
        );
        // Vertices are compacted, so the count must not grow.
        assert!(out.nodes.len() <= m.nodes.len());
    }

    #[test]
    fn simplify_keeps_planar_geometry_in_bounds() {
        // Simplifying a flat plane must keep every surviving vertex on
        // the z=0 plane and inside the original XY bounds (bounded
        // error: a coplanar collapse introduces zero geometric error).
        let m = subdivided_plane(16);
        let (min0, max0) = aabb(&m);
        let out = simplify(&m, 0.3);
        assert!(!out.nodes.is_empty());
        for n in &out.nodes {
            assert!(n.z.abs() < 1e-6, "vertex left the z=0 plane: {n:?}");
            assert!(n.x >= min0.x - 1e-6 && n.x <= max0.x + 1e-6);
            assert!(n.y >= min0.y - 1e-6 && n.y <= max0.y + 1e-6);
        }
    }

    #[test]
    fn simplify_ratio_one_is_essentially_noop() {
        // target_ratio = 1.0 -> target index count == current, so the
        // simplifier has no reduction to do; triangle count is
        // preserved and geometry is unchanged.
        let m = subdivided_plane(8);
        let out = simplify(&m, 1.0);
        assert_eq!(tri_count(&out), tri_count(&m));
    }

    #[test]
    fn simplify_clamps_out_of_range_ratio() {
        // Ratios outside [0,1] are clamped, not panicking.
        let m = subdivided_plane(6);
        let hi = simplify(&m, 5.0); // clamps to 1.0 -> noop-ish
        assert_eq!(tri_count(&hi), tri_count(&m));
        let lo = simplify(&m, -1.0); // clamps to 0.0 -> max reduction
        assert!(tri_count(&lo) <= tri_count(&m));
    }
}
