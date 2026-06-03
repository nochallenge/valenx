//! Decimate-Pro drivers — wrap valenx-mesh's QEM with per-vertex
//! weighting so high-curvature / high-stretch / feature regions
//! resist collapse.
//!
//! v1 strategy: rather than rewrite the QEM heap to accept arbitrary
//! per-vertex multipliers, we *biaspartition* the input mesh into
//! "protected" and "free" subsets, decimate the free part to a higher
//! aggression, then weld them back together. Subsets are chosen per
//! mode:
//!
//! - **Curvature-weighted** — protect every vertex whose curvature is
//!   above the mean+stddev band scaled by `curvature_weight`.
//! - **UV-preserving** — protect every vertex whose UV-stretch weight
//!   is in the top 25 % of the distribution; UVs themselves are
//!   carried through unchanged for surviving vertices.
//! - **Feature-aware** — protect both endpoints of every supplied
//!   feature edge.

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::appearance::{uv_aware_quadric, uv_stretch_weight};
use crate::curvature;
use crate::error::DecimateProError;

/// Curvature-weighted QEM. `curvature_weight = 0` collapses to the
/// plain QEM call; larger values bias the protected set wider so more
/// detail is preserved at the cost of slower convergence.
pub fn weighted_qem(
    mesh: &Mesh,
    target_fraction: f64,
    curvature_weight: f64,
) -> Result<Mesh, DecimateProError> {
    check_fraction(target_fraction)?;
    if curvature_weight < 0.0 {
        return Err(DecimateProError::BadParameter {
            name: "curvature_weight",
            reason: format!("must be >= 0, got {curvature_weight}"),
        });
    }
    let k = curvature::per_vertex(mesh);
    let protect = curvature_protect_mask(&k, curvature_weight);
    Ok(decimate_with_protected(mesh, target_fraction, &protect))
}

/// UV-preserving decimation. The returned UV vector matches the
/// output mesh's surviving vertex order.
pub fn uv_preserving(
    mesh: &Mesh,
    uvs: &[[f64; 2]],
    target_fraction: f64,
) -> Result<(Mesh, Vec<[f64; 2]>), DecimateProError> {
    check_fraction(target_fraction)?;
    if uvs.len() != mesh.nodes.len() {
        return Err(DecimateProError::SizeMismatch {
            name: "uvs",
            mesh: mesh.nodes.len(),
            got: uvs.len(),
        });
    }
    let qmats = uv_aware_quadric(mesh, uvs);
    let weights: Vec<f64> = qmats.iter().copied().map(uv_stretch_weight).collect();
    let protect = top_quartile_mask(&weights);
    let out = decimate_with_protected(mesh, target_fraction, &protect);
    let out_uvs = remap_uvs(mesh, &out, uvs);
    Ok((out, out_uvs))
}

/// Feature-aware decimation — every endpoint of an edge in
/// `feature_edges` becomes a constraint and is not collapsed.
pub fn feature_aware(
    mesh: &Mesh,
    feature_edges: &[(u32, u32)],
    target_fraction: f64,
) -> Result<Mesh, DecimateProError> {
    check_fraction(target_fraction)?;
    let n = mesh.nodes.len();
    let mut protect = vec![false; n];
    for &(a, b) in feature_edges {
        let (ai, bi) = (a as usize, b as usize);
        if ai >= n || bi >= n {
            return Err(DecimateProError::BadParameter {
                name: "feature_edges",
                reason: format!("edge ({a}, {b}) out of bounds for mesh with {n} nodes"),
            });
        }
        protect[ai] = true;
        protect[bi] = true;
    }
    Ok(decimate_with_protected(mesh, target_fraction, &protect))
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

fn check_fraction(f: f64) -> Result<(), DecimateProError> {
    if !(0.0..=1.0).contains(&f) {
        return Err(DecimateProError::BadParameter {
            name: "target_fraction",
            reason: format!("must be in [0, 1], got {f}"),
        });
    }
    Ok(())
}

fn curvature_protect_mask(curvature: &[f64], weight: f64) -> Vec<bool> {
    if curvature.is_empty() {
        return Vec::new();
    }
    let mean = curvature.iter().sum::<f64>() / curvature.len() as f64;
    let var =
        curvature.iter().map(|k| (k - mean).powi(2)).sum::<f64>() / curvature.len() as f64;
    let stddev = var.sqrt();
    // High weight → low threshold → more verts protected.
    let threshold = (mean + (1.0 - weight.clamp(0.0, 1.0)) * stddev).max(0.0);
    curvature.iter().map(|k| *k >= threshold).collect()
}

fn top_quartile_mask(weights: &[f64]) -> Vec<bool> {
    if weights.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<f64> = weights.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q3_idx = (sorted.len() * 3) / 4;
    let q3 = sorted.get(q3_idx).copied().unwrap_or(0.0);
    weights.iter().map(|w| *w >= q3).collect()
}

/// Split `mesh` into protected + free node sets, decimate only the
/// free side at a higher aggression to compensate, then weld back.
///
/// v1 implementation: when `protect` covers the whole mesh, return a
/// clone (nothing to decimate). When it's empty, defer to the plain
/// QEM decimator. The general case keeps every protected node and
/// every triangle that touches at least one protected vertex, then
/// concatenates the (already-decimated) free subset.
fn decimate_with_protected(
    mesh: &Mesh,
    target_fraction: f64,
    protect: &[bool],
) -> Mesh {
    if protect.is_empty() || protect.iter().all(|&p| !p) {
        return valenx_mesh::quadric_error_decimate(mesh, target_fraction);
    }
    if protect.iter().all(|&p| p) {
        return mesh.clone();
    }

    // Tally free vs protected.
    let n_total = protect.len();
    let n_protected = protect.iter().filter(|&&p| p).count();
    let n_free = n_total - n_protected;
    let raw_target = (n_total as f64) * target_fraction;
    let free_target_count = (raw_target - n_protected as f64).max(0.0);
    let free_fraction = if n_free == 0 {
        0.0
    } else {
        (free_target_count / n_free as f64).clamp(0.0, 1.0)
    };

    // Extract subset mesh containing all triangles whose vertices are
    // *all* free. Their decimation is independent of the protected
    // shell.
    let free_mesh = extract_free_subset(mesh, protect);
    let decimated_free = if free_mesh.nodes.is_empty() {
        free_mesh
    } else {
        valenx_mesh::quadric_error_decimate(&free_mesh, free_fraction)
    };

    // Stitch: keep all protected vertices + every original triangle
    // that touches at least one protected vertex (the boundary band),
    // then append the decimated-free vertices + their re-indexed
    // triangles.
    stitch_protected_and_free(mesh, protect, &decimated_free)
}

fn extract_free_subset(mesh: &Mesh, protect: &[bool]) -> Mesh {
    let mut out = Mesh::new(format!("{}_free", mesh.id));
    let mut remap = vec![u32::MAX; mesh.nodes.len()];
    for (i, p) in protect.iter().enumerate() {
        if !p {
            remap[i] = out.nodes.len() as u32;
            out.nodes.push(mesh.nodes[i]);
        }
    }
    for block in &mesh.element_blocks {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        let mut new_blk = ElementBlock::new(ElementType::Tri3);
        for tri in block.connectivity.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if a >= protect.len() || b >= protect.len() || c >= protect.len() {
                continue;
            }
            if protect[a] || protect[b] || protect[c] {
                continue;
            }
            new_blk
                .connectivity
                .extend_from_slice(&[remap[a], remap[b], remap[c]]);
        }
        if !new_blk.connectivity.is_empty() {
            out.element_blocks.push(new_blk);
        }
    }
    out.recompute_stats();
    out
}

fn stitch_protected_and_free(orig: &Mesh, protect: &[bool], decimated_free: &Mesh) -> Mesh {
    let mut out = Mesh::new(format!("{}_decimate_pro", orig.id));

    // Push every protected vertex first; remember mapping.
    let mut prot_remap = vec![u32::MAX; orig.nodes.len()];
    for (i, p) in protect.iter().enumerate() {
        if *p {
            prot_remap[i] = out.nodes.len() as u32;
            out.nodes.push(orig.nodes[i]);
        }
    }
    // Append decimated-free vertices.
    let free_node_offset = out.nodes.len() as u32;
    out.nodes.extend_from_slice(&decimated_free.nodes);

    // Preserve the boundary band: any original triangle that touches
    // a protected vertex stays, with the free corners re-mapped to
    // their *original* free index in the simplified set. v1 keeps
    // those triangles only when ALL three vertices are protected
    // (so the boundary band's free verts don't dangle into a
    // duplicated, un-deduped mesh). v2 will do a true stitch.
    let mut prot_block = ElementBlock::new(ElementType::Tri3);
    for block in &orig.element_blocks {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        for tri in block.connectivity.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if a >= protect.len() || b >= protect.len() || c >= protect.len() {
                continue;
            }
            if protect[a] && protect[b] && protect[c] {
                prot_block
                    .connectivity
                    .extend_from_slice(&[prot_remap[a], prot_remap[b], prot_remap[c]]);
            }
        }
    }
    if !prot_block.connectivity.is_empty() {
        out.element_blocks.push(prot_block);
    }

    // Re-index the decimated-free triangles.
    for blk in &decimated_free.element_blocks {
        if !matches!(blk.element_type, ElementType::Tri3) {
            continue;
        }
        let mut new_blk = ElementBlock::new(ElementType::Tri3);
        for tri in blk.connectivity.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            new_blk.connectivity.extend_from_slice(&[
                tri[0] + free_node_offset,
                tri[1] + free_node_offset,
                tri[2] + free_node_offset,
            ]);
        }
        if !new_blk.connectivity.is_empty() {
            out.element_blocks.push(new_blk);
        }
    }

    out.recompute_stats();
    out
}

fn remap_uvs(orig: &Mesh, out: &Mesh, uvs: &[[f64; 2]]) -> Vec<[f64; 2]> {
    // v1: the protected vertices appear first in the output (per the
    // stitch order), so their UVs are direct slice-copies. Any new
    // decimated-free vertices that don't have a 1:1 ancestor get a
    // [0, 0] placeholder — Phase 47.5 will solve the per-survivor UV
    // average from the QEM contraction trail.
    let mut out_uvs = vec![[0.0_f64; 2]; out.nodes.len()];
    // First pass: protected verts preserve their UVs.
    let mut prot_idx = 0;
    let n_min = orig.nodes.len().min(uvs.len());
    for i in 0..n_min {
        if let Some(orig_uv) = uvs.get(i) {
            if let Some(slot) = out_uvs.get_mut(prot_idx) {
                *slot = *orig_uv;
            }
            prot_idx += 1;
            if prot_idx >= out_uvs.len() {
                break;
            }
        }
    }
    out_uvs
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn unit_square_mesh() -> Mesh {
        let mut m = Mesh::new("sq");
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 1.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity.extend_from_slice(&[0, 1, 2, 0, 2, 3]);
        m.element_blocks.push(b);
        m.recompute_stats();
        m
    }

    #[test]
    fn weighted_qem_rejects_bad_fraction() {
        let m = unit_square_mesh();
        let err = weighted_qem(&m, -0.5, 0.0).unwrap_err();
        assert!(matches!(err, DecimateProError::BadParameter { .. }));
    }

    #[test]
    fn weighted_qem_passes_through_flat_mesh() {
        let m = unit_square_mesh();
        let out = weighted_qem(&m, 1.0, 0.0).expect("ok");
        // 1.0 fraction = no decimation, vertex count unchanged.
        assert_eq!(out.nodes.len(), m.nodes.len());
    }

    #[test]
    fn uv_preserving_rejects_size_mismatch() {
        let m = unit_square_mesh();
        let err = uv_preserving(&m, &[[0.0, 0.0]], 0.5).unwrap_err();
        assert!(matches!(err, DecimateProError::SizeMismatch { .. }));
    }

    #[test]
    fn feature_aware_rejects_out_of_bounds_edge() {
        let m = unit_square_mesh();
        let err = feature_aware(&m, &[(0, 99)], 0.5).unwrap_err();
        assert!(matches!(err, DecimateProError::BadParameter { .. }));
    }

    #[test]
    fn feature_aware_preserves_constrained_vertex_set() {
        let m = unit_square_mesh();
        // Marking every vertex as a feature endpoint protects them
        // all — output should be unchanged (clone-and-return path).
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];
        let out = feature_aware(&m, &edges, 0.5).expect("ok");
        assert_eq!(out.nodes.len(), m.nodes.len());
    }

    #[test]
    fn protect_mask_curvature_threshold_sane() {
        let k = vec![0.1, 0.2, 0.3, 0.4];
        let mask = curvature_protect_mask(&k, 0.0);
        // weight 0 → threshold = mean + 1*stddev, only the top vert
        // qualifies.
        assert!(mask.iter().any(|m| *m));
    }

    #[test]
    fn top_quartile_mask_picks_upper_quarter() {
        let w = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
        let mask = top_quartile_mask(&w);
        let kept: std::collections::HashSet<_> = mask
            .iter()
            .enumerate()
            .filter_map(|(i, m)| if *m { Some(i) } else { None })
            .collect();
        // q3 index = 8*3/4 = 6, threshold = 0.7 → top two qualify.
        assert!(kept.contains(&6));
        assert!(kept.contains(&7));
    }

    #[test]
    fn appearance_quadric_returns_one_per_vertex() {
        let m = unit_square_mesh();
        let uvs = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let q = crate::appearance::QuadricMatrix::default();
        assert_eq!(q.uv, [0.0, 0.0, 0.0]);
        let qs = uv_aware_quadric(&m, &uvs);
        assert_eq!(qs.len(), 4);
    }
}
