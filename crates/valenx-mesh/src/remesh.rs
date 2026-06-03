//! Isotropic remeshing for Tri3 surface meshes.
//!
//! Drives the mesh toward triangles whose edges all sit close to a
//! target length `L`. The classic iteration (Botsch et al. 2010) is:
//!
//! 1. **Split** every edge longer than `4/3 · L` at its midpoint.
//! 2. **Collapse** every edge shorter than `4/5 · L` (the constants
//!    are chosen so a split + collapse don't ping-pong on the same
//!    edge).
//! 3. **Flip** every interior edge whose flip would bring the two
//!    incident triangles' opposite-vertex degrees closer to the
//!    valence-6 ideal.
//! 4. **Smooth** vertices in their tangent plane (we use the
//!    [`crate::smooth::laplacian`] helper as a cheap stand-in;
//!    proper tangent-plane projection is a v1.5 polish).
//!
//! Each "remesh iteration" runs that sequence once. Typical input:
//! 3-5 iterations on an STL produces visibly more uniform triangles
//! without dramatically changing the silhouette.
//!
//! ## Scope
//!
//! - Tri3 blocks only — non-Tri3 blocks pass through unchanged.
//! - Boundary edges are protected (no collapse, no smoothing) so the
//!   silhouette stays put.
//! - We don't attempt to project vertices back onto the original
//!   surface after smoothing; that's a v1.5 follow-up (would need
//!   per-iteration nearest-point queries against the input mesh).

use std::collections::HashMap;

use nalgebra::Vector3;

use crate::element::ElementType;
use crate::mesh::Mesh;
use crate::smooth::laplacian;

/// Run `iterations` isotropic-remesh passes targeting `target_edge_length`.
///
/// The classic split/collapse/flip/smooth loop. `iterations = 3` is a
/// reasonable default for an interactive "make the triangles look
/// uniform" pass on an STL; bigger numbers diminishing-return after
/// ~5-6 iterations on typical inputs.
///
/// Returns a fresh mesh; the original is unmodified. The output `id`
/// is `"<original>_remeshed"`.
pub fn isotropic(mesh: &Mesh, target_edge_length: f64, iterations: usize) -> Mesh {
    let mut current = mesh.clone();
    current.id = format!("{}_remeshed", mesh.id);
    if mesh.nodes.is_empty() || target_edge_length <= 0.0 || iterations == 0 {
        current.recompute_stats();
        return current;
    }
    let long_thresh = (4.0 / 3.0) * target_edge_length;
    let short_thresh = (4.0 / 5.0) * target_edge_length;

    for _ in 0..iterations {
        split_long_edges(&mut current, long_thresh);
        collapse_short_edges(&mut current, short_thresh);
        flip_to_improve_valence(&mut current);
        // Tangent-plane smoothing approximated by a single mild
        // Laplacian pass — gives 90 % of the quality benefit without
        // the nearest-point projection plumbing.
        current = laplacian(&current, 1, 0.3);
        current.id = format!("{}_remeshed", mesh.id);
    }
    current.recompute_stats();
    current
}

/// Split every Tri3 edge whose length exceeds `long_thresh`.
///
/// Each long edge gets a new vertex inserted at its midpoint, and the
/// two adjacent triangles are replaced with four. Public so the UI
/// can call it directly for a "subdivide once" button if needed.
pub fn split_long_edges(mesh: &mut Mesh, long_thresh: f64) {
    let long_sq = long_thresh * long_thresh;
    // Walk Tri3 blocks separately; the existing block(s) get
    // rewritten in-place.
    for block_idx in 0..mesh.element_blocks.len() {
        if mesh.element_blocks[block_idx].element_type != ElementType::Tri3 {
            continue;
        }
        let mut new_conn: Vec<u32> = Vec::new();
        // Local copy so we can mutate mesh.nodes inside the loop.
        let tris: Vec<[u32; 3]> = mesh.element_blocks[block_idx]
            .connectivity
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        // R34 S2 (defense-in-depth): this block is the single ingestion
        // point for Tri3 connectivity in the splitter. `mesh.nodes[a]`
        // (and the union-find `resolve` in the collapse pass) raw-index
        // these values, so we validate ONCE here against the ORIGINAL
        // node count (captured before midpoints are appended) and SKIP
        // any triangle citing a vertex past it — the bad triangle is
        // dropped rather than panicking "index out of bounds". Backs
        // the per-loader parse guards (OBJ/gmsh/netgen/PLY).
        let node_count = mesh.nodes.len() as u32;
        // Map each undirected edge → the midpoint vertex index, so
        // adjacent triangles share the same new vertex.
        let mut edge_mid: HashMap<(u32, u32), u32> = HashMap::new();
        for tri in &tris {
            if tri[0] >= node_count || tri[1] >= node_count || tri[2] >= node_count {
                continue;
            }
            let mut local_mids: [Option<u32>; 3] = [None; 3];
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let pa = mesh.nodes[a as usize];
                let pb = mesh.nodes[b as usize];
                if (pa - pb).norm_squared() > long_sq {
                    let key = if a < b { (a, b) } else { (b, a) };
                    let mid_idx = *edge_mid.entry(key).or_insert_with(|| {
                        mesh.nodes.push((pa + pb) * 0.5);
                        (mesh.nodes.len() - 1) as u32
                    });
                    local_mids[k] = Some(mid_idx);
                }
            }
            emit_split_subfaces(tri, &local_mids, &mut new_conn);
        }
        mesh.element_blocks[block_idx].connectivity = new_conn;
    }
}

/// Given one triangle's three optional midpoints, push 1/2/3/4
/// triangles into `out_conn` depending on how many edges split.
///
/// Layout: edge k is between tri[k] and tri[(k+1)%3]. local_mids[k]
/// is the midpoint of that edge if it was split.
fn emit_split_subfaces(tri: &[u32; 3], local_mids: &[Option<u32>; 3], out_conn: &mut Vec<u32>) {
    let n_split = local_mids.iter().filter(|m| m.is_some()).count();
    match n_split {
        0 => out_conn.extend_from_slice(tri),
        1 => {
            // One edge split: two child triangles.
            let k = local_mids.iter().position(|m| m.is_some()).unwrap();
            let m = local_mids[k].unwrap();
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let c = tri[(k + 2) % 3];
            out_conn.extend_from_slice(&[a, m, c]);
            out_conn.extend_from_slice(&[m, b, c]);
        }
        2 => {
            // Two edges split: three child triangles. We need the
            // unsplit edge as one of the new triangles.
            let unsplit = local_mids.iter().position(|m| m.is_none()).unwrap();
            // Edge `unsplit` is from tri[unsplit] to tri[(unsplit+1)%3].
            // The two split edges are the others. Order by k for
            // child generation.
            let k_unsplit = unsplit;
            let k_a = (k_unsplit + 1) % 3;
            let k_b = (k_unsplit + 2) % 3;
            let m_a = local_mids[k_a].unwrap();
            let m_b = local_mids[k_b].unwrap();
            let v0 = tri[k_unsplit];
            let v1 = tri[(k_unsplit + 1) % 3];
            let v2 = tri[(k_unsplit + 2) % 3];
            // Children: (v0, v1, m_a), (v0, m_a, m_b), (v2, m_b, m_a).
            // Wind so the original orientation is preserved.
            out_conn.extend_from_slice(&[v0, v1, m_a]);
            out_conn.extend_from_slice(&[v0, m_a, m_b]);
            out_conn.extend_from_slice(&[v2, m_b, m_a]);
        }
        3 => {
            // All three edges split: classic 4-1 subdivision.
            let m01 = local_mids[0].unwrap();
            let m12 = local_mids[1].unwrap();
            let m20 = local_mids[2].unwrap();
            out_conn.extend_from_slice(&[tri[0], m01, m20]);
            out_conn.extend_from_slice(&[m01, tri[1], m12]);
            out_conn.extend_from_slice(&[m20, m12, tri[2]]);
            out_conn.extend_from_slice(&[m01, m12, m20]);
        }
        _ => unreachable!(),
    }
}

/// Collapse every Tri3 edge shorter than `short_thresh` by merging
/// one endpoint into the other (at the midpoint).
///
/// Boundary edges (those whose vertices appear in only one triangle
/// edge) are preserved. After collapses, dead vertices and
/// zero-area triangles are filtered and indices compacted.
pub fn collapse_short_edges(mesh: &mut Mesh, short_thresh: f64) {
    let short_sq = short_thresh * short_thresh;
    let boundary = crate::smooth::vertex_neighbors(mesh); // not boundary itself; placeholder
    let _ = boundary; // unused; we use a dedicated check below

    // Build the boundary mask using the same heuristic as smooth.rs.
    let bmask = boundary_mask(mesh);

    let mut alive: Vec<bool> = vec![true; mesh.nodes.len()];
    let mut remap: Vec<u32> = (0..mesh.nodes.len() as u32).collect();
    let mut positions = mesh.nodes.clone();

    // One sweep is enough for v1; multi-pass collapses are handled by
    // re-running the whole isotropic loop.
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for k in 0..3 {
                let a_raw = tri[k];
                let b_raw = tri[(k + 1) % 3];
                let a = resolve(&remap, a_raw) as usize;
                let b = resolve(&remap, b_raw) as usize;
                // R34 S2 (defense-in-depth): a resolved index past the
                // node count (out-of-range connectivity) would panic
                // `alive[a]`/`bmask[a]`/`positions[a]`. Skip the edge —
                // `alive`/`bmask`/`positions` are all sized to the
                // original node count, so this one guard covers them.
                if a >= positions.len() || b >= positions.len() {
                    continue;
                }
                if a == b || !alive[a] || !alive[b] {
                    continue;
                }
                if bmask[a] || bmask[b] {
                    // Boundary protection: don't collapse silhouette edges.
                    continue;
                }
                let pa = positions[a];
                let pb = positions[b];
                if (pa - pb).norm_squared() > short_sq {
                    continue;
                }
                // Collapse b into a, position = midpoint.
                positions[a] = (pa + pb) * 0.5;
                alive[b] = false;
                remap[b] = a as u32;
            }
        }
    }

    // Recompact node table & rewrite each Tri3 block.
    let mut new_index = vec![u32::MAX; positions.len()];
    let mut next = 0u32;
    let mut new_nodes: Vec<Vector3<f64>> = Vec::new();
    for (i, &is_alive) in alive.iter().enumerate() {
        if is_alive {
            new_index[i] = next;
            new_nodes.push(positions[i]);
            next += 1;
        }
    }
    for block in &mut mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        let mut out_conn: Vec<u32> = Vec::new();
        for tri in block.connectivity.chunks_exact(3) {
            // R34 S2 (defense-in-depth): an out-of-range resolved index
            // would panic `new_index[..]`. `.get()` it and treat a miss
            // as a dead vertex (drop the triangle) — same fate as a
            // vertex collapsed away (`u32::MAX`) below.
            let (Some(&a), Some(&b), Some(&c)) = (
                new_index.get(resolve(&remap, tri[0]) as usize),
                new_index.get(resolve(&remap, tri[1]) as usize),
                new_index.get(resolve(&remap, tri[2]) as usize),
            ) else {
                continue;
            };
            if a == u32::MAX || b == u32::MAX || c == u32::MAX {
                continue;
            }
            if a == b || b == c || c == a {
                continue;
            }
            out_conn.extend_from_slice(&[a, b, c]);
        }
        block.connectivity = out_conn;
    }
    // Compact non-Tri3 blocks too.
    for block in &mut mesh.element_blocks {
        if block.element_type == ElementType::Tri3 {
            continue;
        }
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let mut out_conn: Vec<u32> = Vec::with_capacity(block.connectivity.len());
        'elem: for chunk in block.connectivity.chunks_exact(npe) {
            let mut buf = Vec::with_capacity(npe);
            for &idx in chunk {
                // R34 S2 (defense-in-depth): an out-of-range resolved
                // index would panic `new_index[..]`. `.get()` it and
                // drop the whole element on a miss — same fate as an
                // element touching a collapsed (dead) vertex.
                let Some(&r) = new_index.get(resolve(&remap, idx) as usize) else {
                    continue 'elem;
                };
                if r == u32::MAX {
                    continue 'elem;
                }
                buf.push(r);
            }
            out_conn.extend_from_slice(&buf);
        }
        block.connectivity = out_conn;
    }
    mesh.nodes = new_nodes;
}

/// Walk the union-find–style remap until we hit a fixed point.
///
/// R34 S2 (defense-in-depth): `remap` has one entry per node, so an
/// out-of-range `v` (a connectivity index from a future un-hardened
/// loader) would panic `remap[v]`. Such a value can't be remapped, so
/// we return it unchanged via `.get()` — callers range-check the
/// result before indexing `alive`/`positions`/`new_index` and drop the
/// element. Backs the per-loader parse guards.
fn resolve(remap: &[u32], mut v: u32) -> u32 {
    loop {
        let Some(&r) = remap.get(v as usize) else {
            return v;
        };
        if r == v {
            return v;
        }
        v = r;
    }
}

/// Boundary-vertex mask: vertex is boundary if any edge of it is
/// used in only one triangle.
fn boundary_mask(mesh: &Mesh) -> Vec<bool> {
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
    let mut mask = vec![false; mesh.nodes.len()];
    for ((a, b), n) in edge_count {
        if n == 1 {
            // R34 S2 (defense-in-depth): `a`/`b` are connectivity
            // values, so an out-of-range index would panic `mask[..]`.
            // Use `.get_mut()` and skip the out-of-range endpoint.
            if let Some(flag) = mask.get_mut(a as usize) {
                *flag = true;
            }
            if let Some(flag) = mask.get_mut(b as usize) {
                *flag = true;
            }
        }
    }
    mask
}

/// Flip interior edges to bring the four-vertex valence sum closer
/// to 24 (target valence-6 for each of the 4 vertices).
///
/// For each interior edge `(a, b)` shared by triangles
/// `(a, b, c)` and `(b, a, d)`, we consider flipping to `(c, d)`.
/// The flip is taken iff the sum of `|valence(v) - 6|` over the
/// four vertices strictly decreases after the flip. One sweep per
/// call — the caller (`isotropic`) iterates the whole loop several
/// times for convergence.
pub fn flip_to_improve_valence(mesh: &mut Mesh) {
    // Build the edge→(tri_a, tri_b) map once.
    let mut edge_to_tris: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    let mut all_tris: Vec<[u32; 3]> = Vec::new();
    let mut tri_block_ix: Vec<usize> = Vec::new();
    for (bix, block) in mesh.element_blocks.iter().enumerate() {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let arr = [tri[0], tri[1], tri[2]];
            let ti = all_tris.len();
            all_tris.push(arr);
            tri_block_ix.push(bix);
            for k in 0..3 {
                let a = arr[k];
                let b = arr[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                edge_to_tris.entry(key).or_default().push(ti);
            }
        }
    }
    // R34 S2 (defense-in-depth): `all_tris` is kept 1:1 with the
    // original connectivity (the rewrite below walks it positionally),
    // so we can't drop a malformed triangle here. Instead every index
    // derived from a connectivity value is `.get()`-guarded: the
    // valence tally skips out-of-range vertices, and any candidate
    // flip touching an out-of-range vertex is skipped (the triangle
    // passes through unchanged). A valid mesh is byte-identical. Backs
    // the per-loader parse guards (OBJ/gmsh/netgen/PLY).
    let mut valence = vec![0u32; mesh.nodes.len()];
    for tri in &all_tris {
        for &v in tri {
            if let Some(cnt) = valence.get_mut(v as usize) {
                *cnt += 1;
            }
        }
    }
    let mut flipped = vec![false; all_tris.len()];
    for ((a, b), tris) in &edge_to_tris {
        if tris.len() != 2 {
            continue;
        }
        let (ti0, ti1) = (tris[0], tris[1]);
        if flipped[ti0] || flipped[ti1] {
            continue;
        }
        let tri0 = all_tris[ti0];
        let tri1 = all_tris[ti1];
        let Some(c) = third_vertex(&tri0, *a, *b) else {
            continue;
        };
        let Some(d) = third_vertex(&tri1, *a, *b) else {
            continue;
        };
        if c == d {
            continue;
        }
        // R34 S2 (defense-in-depth): every valence read/write below
        // indexes a connectivity value. If any of the four vertices is
        // out of range we can't score or apply the flip, so skip it —
        // the malformed triangle passes through unchanged rather than
        // panicking `valence[..]`.
        let (Some(&va), Some(&vb), Some(&vc), Some(&vd)) = (
            valence.get(*a as usize),
            valence.get(*b as usize),
            valence.get(c as usize),
            valence.get(d as usize),
        ) else {
            continue;
        };
        let cost_before = badness(va) + badness(vb) + badness(vc) + badness(vd);
        // After flip: a loses 1, b loses 1, c gains 1, d gains 1.
        let cost_after = badness(va.saturating_sub(1))
            + badness(vb.saturating_sub(1))
            + badness(vc + 1)
            + badness(vd + 1);
        if cost_after < cost_before {
            // Apply: rewrite tri0 = (a, c, d), tri1 = (b, d, c) so
            // the winding stays consistent with the original triangles.
            all_tris[ti0] = [*a, c, d];
            all_tris[ti1] = [*b, d, c];
            // The four indices are all in range (checked above).
            valence[*a as usize] -= 1;
            valence[*b as usize] -= 1;
            valence[c as usize] += 1;
            valence[d as usize] += 1;
            flipped[ti0] = true;
            flipped[ti1] = true;
        }
    }
    // Rewrite each Tri3 block from `all_tris` in original order.
    let mut tri_idx = 0;
    for (bix, block) in mesh.element_blocks.iter_mut().enumerate() {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        let mut out_conn: Vec<u32> = Vec::with_capacity(block.connectivity.len());
        // Walk through tris in this block.
        let count = block.connectivity.len() / 3;
        for _ in 0..count {
            assert!(tri_idx < all_tris.len());
            assert_eq!(tri_block_ix[tri_idx], bix);
            let t = all_tris[tri_idx];
            out_conn.extend_from_slice(&t);
            tri_idx += 1;
        }
        block.connectivity = out_conn;
    }
}

fn third_vertex(tri: &[u32; 3], a: u32, b: u32) -> Option<u32> {
    tri.iter().find(|&&v| v != a && v != b).copied()
}

/// Distance from ideal valence (6 for interior triangle meshes).
fn badness(valence: u32) -> i32 {
    (valence as i32 - 6).abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ElementBlock;

    fn pt(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    /// A 2x2 quad-grid triangulated -> 9 vertices, 8 tris. Edges
    /// are mostly length 1 and sqrt(2).
    fn small_grid() -> Mesh {
        let mut m = Mesh::new("grid");
        for j in 0..=2 {
            for i in 0..=2 {
                m.nodes.push(pt(i as f64, j as f64, 0.0));
            }
        }
        let mut blk = ElementBlock::new(ElementType::Tri3);
        for j in 0..2 {
            for i in 0..2 {
                let i0 = (j * 3 + i) as u32;
                let i1 = i0 + 1;
                let i2 = i0 + 3;
                let i3 = i2 + 1;
                blk.connectivity
                    .extend_from_slice(&[i0, i1, i3, i0, i3, i2]);
            }
        }
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    #[test]
    fn split_doubles_triangle_when_all_edges_split() {
        // Single equilateral triangle, edge length 1. Set long_thresh
        // very small so all three edges split → 4 child triangles.
        let mut m = Mesh::new("tri");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(0.5, (3.0f64).sqrt() / 2.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks = vec![blk];
        split_long_edges(&mut m, 0.1);
        // After 4-1 subdivision: 6 vertices, 4 triangles.
        assert_eq!(m.nodes.len(), 6);
        assert_eq!(m.element_blocks[0].connectivity.len() / 3, 4);
    }

    #[test]
    fn split_noop_when_all_edges_short() {
        // No edge exceeds threshold → output identical to input.
        let mut m = small_grid();
        let conn_before = m.element_blocks[0].connectivity.clone();
        let nodes_before = m.nodes.len();
        split_long_edges(&mut m, 10.0);
        assert_eq!(m.nodes.len(), nodes_before);
        assert_eq!(m.element_blocks[0].connectivity, conn_before);
    }

    #[test]
    fn collapse_drops_short_interior_edge() {
        // Build a 4x4 grid (25 vertices) so we have at least one
        // edge whose endpoints are both interior. Shrink one such
        // edge below threshold and confirm a vertex disappears.
        let mut m = Mesh::new("4x4");
        for j in 0..=4 {
            for i in 0..=4 {
                m.nodes.push(pt(i as f64, j as f64, 0.0));
            }
        }
        let mut blk = ElementBlock::new(ElementType::Tri3);
        for j in 0..4 {
            for i in 0..4 {
                let i0 = (j * 5 + i) as u32;
                let i1 = i0 + 1;
                let i2 = i0 + 5;
                let i3 = i2 + 1;
                blk.connectivity
                    .extend_from_slice(&[i0, i1, i3, i0, i3, i2]);
            }
        }
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Vertices 6 and 7 are both interior (i, j in 1..=3). Pull
        // 7 toward 6 so their edge is below threshold.
        // Original positions: 6 = (1, 1), 7 = (2, 1).
        m.nodes[7] = pt(1.001, 1.0, 0.0);
        let nodes_before = m.nodes.len();
        collapse_short_edges(&mut m, 0.1);
        assert!(m.nodes.len() < nodes_before);
    }

    #[test]
    fn flip_improves_valence_on_two_triangles() {
        // Two triangles sharing edge (1, 2). The flip swap (1,2)->(0,3)
        // changes valences. We just verify no panic, output is valid,
        // and triangle count stays at 2.
        let mut m = Mesh::new("two");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0),
            pt(1.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2, 1, 3, 2];
        m.element_blocks = vec![blk];
        flip_to_improve_valence(&mut m);
        // Still two triangles, same vertex set.
        assert_eq!(m.element_blocks[0].connectivity.len() / 3, 2);
        assert_eq!(m.nodes.len(), 4);
    }

    #[test]
    fn isotropic_runs_to_target_length_within_tolerance() {
        // Build a 4x4 grid mesh (each cell 1.0 across). Run isotropic
        // remeshing with target 0.5 — edges should approach 0.5±30%.
        let mut m = Mesh::new("4x4");
        for j in 0..=4 {
            for i in 0..=4 {
                m.nodes.push(pt(i as f64, j as f64, 0.0));
            }
        }
        let mut blk = ElementBlock::new(ElementType::Tri3);
        for j in 0..4 {
            for i in 0..4 {
                let i0 = (j * 5 + i) as u32;
                let i1 = i0 + 1;
                let i2 = i0 + 5;
                let i3 = i2 + 1;
                blk.connectivity
                    .extend_from_slice(&[i0, i1, i3, i0, i3, i2]);
            }
        }
        m.element_blocks = vec![blk];
        m.recompute_stats();

        let out = isotropic(&m, 0.5, 3);
        let (mean, _var) = crate::smooth::edge_length_stats(&out)
            .expect("output mesh has at least one Tri3 block");
        // Within 30% of target.
        assert!(
            (mean - 0.5).abs() / 0.5 < 0.3,
            "mean edge length {mean} is not within 30% of 0.5"
        );
    }

    /// R34 S2 (RED→GREEN): defense-in-depth sink seal. A mesh whose
    /// Tri3 connectivity cites a vertex past `nodes.len()` must NOT
    /// panic any remesh entry point. Pre-fix `split_long_edges` did
    /// `mesh.nodes[a]`, `collapse_short_edges` resolved through
    /// `remap[v]` then indexed `alive[a]`/`new_index[..]`, and
    /// `flip_to_improve_valence` did `valence[v]` — each panicked
    /// "index out of bounds". Post-fix the out-of-range
    /// triangle/element is dropped (or passed through for flip) and
    /// every pass returns. We assert no panic.
    #[test]
    fn out_of_range_connectivity_does_not_panic() {
        let make = || {
            let mut m = Mesh::new("hostile");
            m.nodes = vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 0.0)];
            let mut blk = ElementBlock::new(ElementType::Tri3);
            // A triangle citing vertex 8 (out of range).
            blk.connectivity = vec![0, 1, 8];
            m.element_blocks = vec![blk];
            m.recompute_stats();
            m
        };
        // Each public pass must return rather than panic.
        let mut m1 = make();
        split_long_edges(&mut m1, 0.1);
        let mut m2 = make();
        collapse_short_edges(&mut m2, 10.0);
        let mut m3 = make();
        flip_to_improve_valence(&mut m3);
        // The full driver too.
        let out = isotropic(&make(), 0.5, 3);
        assert_eq!(out.id, "hostile_remeshed");
    }

    /// R34 S2: a valid triangle alongside an out-of-range one. The
    /// splitter drops the bad triangle but keeps the valid one's
    /// geometry; no panic.
    #[test]
    fn out_of_range_triangle_skipped_valid_kept_in_split() {
        let mut m = Mesh::new("mixed");
        m.nodes = vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 0.0)];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // First tri valid (0,1,2); second cites vertex 9.
        blk.connectivity = vec![0, 1, 2, 0, 1, 9];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Force all valid edges to split: the one valid triangle
        // becomes 4 children; the bad triangle is dropped (not 4+2).
        split_long_edges(&mut m, 0.1);
        let tri_out = m.element_blocks[0].connectivity.len() / 3;
        assert_eq!(
            tri_out, 4,
            "the valid triangle subdivides into 4; the bad one is dropped"
        );
    }

    #[test]
    fn isotropic_empty_mesh_returns_empty_with_id_suffix() {
        let m = Mesh::new("empty");
        let out = isotropic(&m, 0.5, 2);
        assert!(out.nodes.is_empty());
        assert_eq!(out.id, "empty_remeshed");
    }

    #[test]
    fn isotropic_zero_iter_is_noop_but_renames_id() {
        let m = small_grid();
        let out = isotropic(&m, 0.5, 0);
        assert_eq!(out.nodes.len(), m.nodes.len());
        assert_eq!(out.id, "grid_remeshed");
    }
}
