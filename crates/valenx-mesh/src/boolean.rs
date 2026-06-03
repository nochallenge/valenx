//! Boolean-style mesh operations.
//!
//! ## Scope
//!
//! This is the v0 scaffold per the workspace's "ship the simple
//! operation, name the gap" pattern. What's here:
//!
//! - **`concatenate`** — naive index-offset merge: stack two meshes'
//!   node arrays, re-index the second mesh's connectivity by the
//!   first mesh's node count, concatenate the element blocks. **Not**
//!   a topological union — overlapping volumes get double-counted, no
//!   surface trimming, no shared-edge merging.
//! - **`merge_coincident_nodes`** — O(N²) pass that snaps nodes
//!   within `tolerance` to a single index and rewrites every
//!   element's connectivity to match. Useful as a post-pass after
//!   concatenate when the inputs share boundary faces.
//! - **`union_concatenate`** — convenience wrapper:
//!   `concatenate` + `merge_coincident_nodes`.
//!
//! ## What's NOT here yet
//!
//! Real CSG (BSP-tree intersection, exact-arithmetic robustness,
//! coplanar-face handling) is its own follow-up — the standard
//! libraries that do this well (CGAL, OpenSCAD, libigl) are LGPL or
//! GPL, so a Rust-native implementation has to come before we can
//! drop a real intersect / subtract path in.
//!
//! For the workshop "load two STLs, view them together" workflow that
//! covers ~80% of user requests today, concatenate + dedup is enough.

use std::collections::HashMap;

use nalgebra::Vector3;

use crate::element::ElementBlock;
use crate::mesh::Mesh;

/// Concatenate two meshes into one. The second mesh's connectivity
/// is re-indexed by `a.nodes.len()` so element references stay
/// consistent with the merged node array.
///
/// `regions` and `boundaries` from `b` are dropped — they reference
/// positional element indices that don't survive the concatenation
/// (the element blocks get appended, but a region indexing
/// `block[3]` would now refer to the wrong block in the merged
/// mesh). Caller must rebuild any region metadata after the merge.
///
/// `id` of the output is `<a.id>+<b.id>` so audit logs / file dumps
/// retain a hint of the inputs.
pub fn concatenate(a: &Mesh, b: &Mesh) -> Mesh {
    let mut out = Mesh::new(format!("{}+{}", a.id, b.id));
    out.nodes.extend(a.nodes.iter().copied());
    out.nodes.extend(b.nodes.iter().copied());

    // First mesh's blocks land verbatim.
    for blk in &a.element_blocks {
        out.element_blocks.push(blk.clone());
    }
    // Second mesh's blocks need each connectivity index shifted by
    // a.nodes.len() so the indices match the new node array.
    let offset = a.nodes.len() as u32;
    for blk in &b.element_blocks {
        let mut shifted = ElementBlock::new(blk.element_type);
        shifted.connectivity = blk.connectivity.iter().map(|i| i + offset).collect();
        out.element_blocks.push(shifted);
    }

    out.regions = a.regions.clone();
    out.boundaries = a.boundaries.clone();
    out.recompute_stats();
    out
}

/// Snap nodes within `tolerance` to a single index, then rewrite
/// every element's connectivity to use the deduplicated index.
/// Tolerance is in the same units as the node coordinates.
///
/// Algorithm: bucket nodes into a coarse spatial hash keyed on
/// `(round(x/h), round(y/h), round(z/h))` where `h = tolerance`.
/// For each candidate node, only the bucket + 26 neighbour buckets
/// need to be searched — keeps the per-node work O(k) where k is
/// the typical bucket occupancy (≈1 for well-distributed meshes).
///
/// Returns a fresh `Mesh` with the dedup applied. The original
/// stays untouched so callers can keep both views.
pub fn merge_coincident_nodes(mesh: &Mesh, tolerance: f64) -> Mesh {
    if tolerance <= 0.0 || mesh.nodes.is_empty() {
        return mesh.clone();
    }
    let h = tolerance.max(1e-30);
    let key = |p: Vector3<f64>| -> (i64, i64, i64) {
        (
            (p.x / h).round() as i64,
            (p.y / h).round() as i64,
            (p.z / h).round() as i64,
        )
    };

    // Bucket nodes by spatial key. Each bucket holds the indices
    // of nodes that fall within it; we only search the bucket and
    // its 26 neighbours to find a coincident match.
    let mut buckets: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    let mut remap: Vec<u32> = Vec::with_capacity(mesh.nodes.len());
    let mut deduped_nodes: Vec<Vector3<f64>> = Vec::new();
    let tol_sq = tolerance * tolerance;

    for (idx, &p) in mesh.nodes.iter().enumerate() {
        let center = key(p);
        let mut found: Option<u32> = None;
        // Search 27 neighbour buckets.
        'outer: for dx in -1..=1i64 {
            for dy in -1..=1i64 {
                for dz in -1..=1i64 {
                    let nb = (center.0 + dx, center.1 + dy, center.2 + dz);
                    if let Some(candidates) = buckets.get(&nb) {
                        for &cand_idx in candidates {
                            let q = deduped_nodes[cand_idx];
                            let dx = p.x - q.x;
                            let dy = p.y - q.y;
                            let dz = p.z - q.z;
                            if dx * dx + dy * dy + dz * dz <= tol_sq {
                                found = Some(cand_idx as u32);
                                break 'outer;
                            }
                        }
                    }
                }
            }
        }
        let new_idx = match found {
            Some(i) => i,
            None => {
                let i = deduped_nodes.len();
                deduped_nodes.push(p);
                buckets.entry(center).or_default().push(i);
                i as u32
            }
        };
        remap.push(new_idx);
        let _ = idx;
    }

    // Rebuild every element block with remapped connectivity.
    //
    // R34 S1 (defense-in-depth): `remap` has exactly one entry per
    // input node, so `remap[*c]` panics if a connectivity index points
    // past the node array. The per-loader parse guards (OBJ/gmsh/
    // netgen/PLY) are the first line of defence; this seal backs them
    // so a future un-hardened loader degrades gracefully instead of
    // panicking here. We rebuild each block element-by-element and drop
    // any element that cites an out-of-range index (`.get()` + skip)
    // rather than remapping a partial, corrupt element.
    let mut out = mesh.clone();
    let node_count = remap.len();
    out.nodes = deduped_nodes;
    for blk in &mut out.element_blocks {
        let npe = blk.element_type.nodes_per_element();
        if npe == 0 {
            blk.connectivity.clear();
            continue;
        }
        let mut rebuilt: Vec<u32> = Vec::with_capacity(blk.connectivity.len());
        'elem: for chunk in blk.connectivity.chunks_exact(npe) {
            let mut out_chunk = Vec::with_capacity(npe);
            for &c in chunk {
                if (c as usize) >= node_count {
                    // Out-of-range index: drop the whole element.
                    continue 'elem;
                }
                out_chunk.push(remap[c as usize]);
            }
            rebuilt.extend_from_slice(&out_chunk);
        }
        blk.connectivity = rebuilt;
    }
    out.recompute_stats();
    out
}

/// Convenience: concatenate `a` and `b`, then merge coincident
/// nodes within `tolerance`. The most common boolean-union approx
/// pipeline for the "load two STLs that share a contact face" case.
pub fn union_concatenate(a: &Mesh, b: &Mesh, tolerance: f64) -> Mesh {
    let merged = concatenate(a, b);
    if tolerance > 0.0 {
        merge_coincident_nodes(&merged, tolerance)
    } else {
        merged
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ElementType;

    fn pt(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    /// Single-tet mesh anchored at the origin.
    fn tet_at_origin() -> Mesh {
        let mut m = Mesh::new("a");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0),
            pt(0.0, 0.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tet4);
        blk.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    /// Single-tet mesh translated by (10, 0, 0) — no spatial overlap
    /// with the origin tet so dedup is a no-op.
    fn tet_far_away() -> Mesh {
        let mut m = Mesh::new("b");
        m.nodes = vec![
            pt(10.0, 0.0, 0.0),
            pt(11.0, 0.0, 0.0),
            pt(10.0, 1.0, 0.0),
            pt(10.0, 0.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tet4);
        blk.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    /// Single-tet mesh that shares the (0,0,0) vertex with
    /// tet_at_origin. After dedup we expect 7 unique nodes
    /// (4 + 4 - 1).
    fn tet_sharing_origin_vertex() -> Mesh {
        let mut m = Mesh::new("c");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0), // shared with tet_at_origin's node 0
            pt(-1.0, 0.0, 0.0),
            pt(0.0, -1.0, 0.0),
            pt(0.0, 0.0, -1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tet4);
        blk.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    #[test]
    fn concatenate_offsets_second_meshs_connectivity() {
        let merged = concatenate(&tet_at_origin(), &tet_far_away());
        assert_eq!(merged.nodes.len(), 8); // 4 + 4
        assert_eq!(merged.element_blocks.len(), 2);
        // First block's connectivity stays at [0, 1, 2, 3].
        assert_eq!(merged.element_blocks[0].connectivity, vec![0, 1, 2, 3]);
        // Second block's connectivity must shift by 4 (a.nodes.len()).
        assert_eq!(merged.element_blocks[1].connectivity, vec![4, 5, 6, 7]);
    }

    #[test]
    fn concatenate_id_carries_both_inputs() {
        let merged = concatenate(&tet_at_origin(), &tet_far_away());
        assert_eq!(merged.id, "a+b");
    }

    #[test]
    fn concatenate_recomputes_stats() {
        let merged = concatenate(&tet_at_origin(), &tet_far_away());
        assert_eq!(merged.stats.node_count, 8);
        assert_eq!(merged.stats.element_count, 2);
    }

    #[test]
    fn merge_coincident_nodes_no_op_when_tolerance_zero() {
        let m = tet_at_origin();
        let merged = merge_coincident_nodes(&m, 0.0);
        assert_eq!(merged.nodes, m.nodes);
        assert_eq!(merged.element_blocks[0].connectivity, vec![0, 1, 2, 3]);
    }

    #[test]
    fn merge_coincident_nodes_combines_duplicate_origin_vertices() {
        // After concatenate, both tets have a vertex at (0,0,0).
        // merge_coincident_nodes should collapse the duplicate.
        let merged = concatenate(&tet_at_origin(), &tet_sharing_origin_vertex());
        assert_eq!(merged.nodes.len(), 8);
        let deduped = merge_coincident_nodes(&merged, 1e-6);
        // 4 + 4 - 1 shared vertex = 7 unique nodes.
        assert_eq!(deduped.nodes.len(), 7);
        // Both original blocks survive, but the second one's
        // connectivity now references the first one's origin index
        // (index 0) instead of the duplicate at index 4.
        assert_eq!(deduped.element_blocks.len(), 2);
        let second_conn = &deduped.element_blocks[1].connectivity;
        // Whichever index ends up as the canonical origin, the
        // second block must reuse it (and the merged mesh must NOT
        // have two distinct origin nodes).
        let origin_count = deduped.nodes.iter().filter(|p| p.norm() < 1e-12).count();
        assert_eq!(origin_count, 1, "expected exactly one (0,0,0) node");
        // The second block's first vertex must point at that
        // origin, NOT a duplicate.
        let origin_idx = deduped.nodes.iter().position(|p| p.norm() < 1e-12).unwrap() as u32;
        assert_eq!(second_conn[0], origin_idx);
    }

    #[test]
    fn merge_coincident_nodes_preserves_distant_vertices() {
        // Two well-separated tets — dedup must NOT collapse anything.
        let merged = concatenate(&tet_at_origin(), &tet_far_away());
        let deduped = merge_coincident_nodes(&merged, 1e-3);
        assert_eq!(deduped.nodes.len(), 8);
    }

    #[test]
    fn merge_coincident_nodes_respects_tolerance() {
        // Two nodes 0.01 apart; tolerance 0.05 collapses them,
        // tolerance 0.001 keeps them distinct.
        let mut m = Mesh::new("close");
        m.nodes = vec![pt(0.0, 0.0, 0.0), pt(0.01, 0.0, 0.0)];
        m.recompute_stats();
        assert_eq!(merge_coincident_nodes(&m, 0.05).nodes.len(), 1);
        assert_eq!(merge_coincident_nodes(&m, 0.001).nodes.len(), 2);
    }

    #[test]
    fn union_concatenate_collapses_shared_vertex() {
        let merged = union_concatenate(&tet_at_origin(), &tet_sharing_origin_vertex(), 1e-6);
        assert_eq!(merged.nodes.len(), 7);
    }

    #[test]
    fn union_concatenate_with_zero_tolerance_skips_dedup() {
        // tolerance = 0 -> just concatenate, no dedup.
        let merged = union_concatenate(&tet_at_origin(), &tet_sharing_origin_vertex(), 0.0);
        assert_eq!(merged.nodes.len(), 8);
    }

    /// R34 S1 (RED→GREEN): defense-in-depth sink seal. A mesh whose
    /// connectivity cites a node index past `nodes.len()` must NOT
    /// panic `merge_coincident_nodes` — the per-loader validation is
    /// the first line, but a future un-hardened loader could hand us
    /// such a mesh. Pre-fix the rewrite did `remap[*c as usize]` and
    /// panicked with "index out of bounds". Post-fix the offending
    /// element index is skipped (the element is dropped) and a result
    /// is returned.
    #[test]
    fn out_of_range_connectivity_does_not_panic() {
        let mut m = Mesh::new("hostile");
        // 2 real nodes...
        m.nodes = vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0)];
        // ...but an element citing node 5 (out of range).
        let mut blk = ElementBlock::new(ElementType::Tet4);
        blk.connectivity = vec![0, 1, 5, 0];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Must return without panicking.
        let out = merge_coincident_nodes(&m, 1e-6);
        // Both nodes are distinct (no dedup), and the malformed element
        // is dropped rather than crashing.
        assert_eq!(out.nodes.len(), 2);
        let total_conn: usize = out.element_blocks.iter().map(|b| b.connectivity.len()).sum();
        assert_eq!(
            total_conn, 0,
            "the element touching an out-of-range index must be dropped"
        );
    }

    /// R34 S1: a valid element survives the dedup remap untouched while
    /// a sibling element with an out-of-range index is dropped.
    #[test]
    fn out_of_range_element_skipped_valid_kept() {
        let merged = concatenate(&tet_at_origin(), &tet_far_away());
        // merged has 8 nodes, two valid Tet4 blocks. Inject a third
        // block citing a node past the array.
        let mut m = merged;
        let mut bad = ElementBlock::new(ElementType::Tet4);
        bad.connectivity = vec![0, 1, 2, 42];
        m.element_blocks.push(bad);
        let out = merge_coincident_nodes(&m, 1e-6);
        // The two valid blocks keep their 4-index connectivity; the bad
        // block is dropped. No panic.
        let four_index_blocks = out
            .element_blocks
            .iter()
            .filter(|b| b.connectivity.len() == 4)
            .count();
        assert_eq!(four_index_blocks, 2, "both valid tets must survive");
    }

    #[test]
    fn merge_coincident_nodes_handles_empty_mesh() {
        let m = Mesh::new("empty");
        let deduped = merge_coincident_nodes(&m, 1e-6);
        assert!(deduped.nodes.is_empty());
        assert!(deduped.element_blocks.is_empty());
    }
}
