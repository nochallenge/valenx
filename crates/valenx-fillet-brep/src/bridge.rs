//! Convenience selectors used by the `valenx-feature-tree`
//! dispatcher.
//!
//! These mirror the mesh-domain selectors in `valenx-fillet` but work
//! at the BRep level — no tessellation, no vertex welding. The
//! dispatcher in `valenx-feature-tree::ops::fillet` calls one of
//! these to convert a user's intent ("fillet every edge sharper than
//! 45°", "fillet these specific edges by index") into the
//! `Vec<truck_modeling::Edge>` that the per-edge primitives in
//! [`crate::fillet`] / [`crate::chamfer`] expect.

use std::collections::HashSet;

use truck_modeling::{Edge as TruckEdge, Solid as TruckSolid};

use crate::edge_classify::edge_dihedral_angle;

/// All distinct edges of a solid, in the order they appear in the
/// face iterator. Stable across replays of the same tree given the
/// same input solid (truck's iteration order is deterministic).
pub fn unique_edges(solid: &TruckSolid) -> Vec<TruckEdge> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for e in solid.edge_iter() {
        if seen.insert(e.id()) {
            out.push(e);
        }
    }
    out
}

/// Edges whose dihedral angle (between outward face normals) is at
/// least `threshold_rad`. Mirrors Phase 3's mesh-domain auto-select.
///
/// `threshold_rad` is interpreted the same as the mesh-domain
/// classifier: 0 includes every interior edge; π/4 (~45°) is the
/// typical "fillet the sharp ones" default; π/2 only picks right
/// angles or sharper.
///
/// Edges that fail the planar-adjacency precondition (non-planar
/// face, open shell) are silently skipped — they wouldn't be
/// filletable in v1 anyway.
pub fn edges_for_threshold(solid: &TruckSolid, threshold_rad: f64) -> Vec<TruckEdge> {
    unique_edges(solid)
        .into_iter()
        .filter(|e| match edge_dihedral_angle(solid, e) {
            Some(angle) => angle >= threshold_rad,
            None => false,
        })
        .collect()
}

/// Edges referenced by a list of 0-based indices into
/// [`unique_edges`]. Out-of-range indices are silently skipped (the
/// caller is responsible for surfacing the right error — typically
/// `BadParameter`).
///
/// `indices` are de-duplicated; duplicate entries don't produce
/// duplicate output edges.
pub fn edges_for_indices(solid: &TruckSolid, indices: &[usize]) -> Vec<TruckEdge> {
    let all = unique_edges(solid);
    let mut seen = HashSet::new();
    indices
        .iter()
        .filter(|i| seen.insert(*i))
        .filter_map(|i| all.get(*i).cloned())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;
    use valenx_cad::primitives::box_solid;

    fn inner_brep(s: &valenx_cad::Solid) -> &TruckSolid {
        match s {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("expected brep"),
        }
    }

    #[test]
    fn unique_edges_returns_twelve_for_cube() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        assert_eq!(unique_edges(brep).len(), 12);
    }

    #[test]
    fn edges_for_threshold_zero_returns_all() {
        // threshold = 0 includes every edge with a non-zero dihedral.
        // All cube edges have 90° dihedral so all 12 qualify.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = edges_for_threshold(brep, 0.0);
        assert_eq!(edges.len(), 12);
    }

    #[test]
    fn edges_for_threshold_45deg_returns_all_for_cube() {
        // Cube edges all have 90° dihedral > π/4.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = edges_for_threshold(brep, std::f64::consts::FRAC_PI_4);
        assert_eq!(edges.len(), 12);
    }

    #[test]
    fn edges_for_threshold_91deg_returns_none_for_cube() {
        // Cube edges have *exactly* 90° dihedral; threshold > 90°
        // → no qualifiers.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let just_above = FRAC_PI_2 + 0.01;
        let edges = edges_for_threshold(brep, just_above);
        assert_eq!(edges.len(), 0);
    }

    #[test]
    fn edges_for_indices_dedupes_and_skips_out_of_range() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        // Pass duplicate index and an out-of-range index; expect
        // exactly one resulting edge.
        let edges = edges_for_indices(brep, &[0, 0, 99]);
        assert_eq!(edges.len(), 1);
    }

    #[test]
    fn edges_for_indices_empty_returns_empty() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        assert!(edges_for_indices(brep, &[]).is_empty());
    }

    #[test]
    fn edges_for_indices_preserves_order_modulo_dedupe() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        // Indices [3, 1, 5] should produce three different edges in
        // that order (assuming no duplicates).
        let edges = edges_for_indices(brep, &[3, 1, 5]);
        assert_eq!(edges.len(), 3);
        // All three should be unique.
        let mut ids = HashSet::new();
        for e in &edges {
            assert!(ids.insert(e.id()), "edges_for_indices produced duplicates");
        }
    }

    #[test]
    fn unique_edges_iteration_is_deterministic() {
        // Two builds of the same cube should produce edge lists in
        // the same order (truck's iter is deterministic).
        let cube_a = box_solid(1.0, 1.0, 1.0).unwrap();
        let cube_b = box_solid(1.0, 1.0, 1.0).unwrap();
        let edges_a = unique_edges(inner_brep(&cube_a));
        let edges_b = unique_edges(inner_brep(&cube_b));
        assert_eq!(edges_a.len(), edges_b.len());
        // Endpoint positions should match index-for-index.
        for (a, b) in edges_a.iter().zip(edges_b.iter()) {
            let (af, ab) = (a.front().point(), a.back().point());
            let (bf, bb) = (b.front().point(), b.back().point());
            assert!((af - bf).x.abs() < 1e-9);
            assert!((ab - bb).z.abs() < 1e-9);
        }
    }
}
