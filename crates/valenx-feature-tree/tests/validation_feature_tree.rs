//! CAD-kernel validation suite — the parametric feature tree.
//!
//! The feature tree is the parametric heart of the kernel: a sequence
//! of features must rebuild deterministically to the *expected* solid,
//! and a parameter change must *propagate* through replay. These tests
//! assert that against **analytic ground truth** — the replayed
//! solid's measured volume must equal the closed-form value the
//! feature sequence implies.
//!
//! What is checked:
//!
//! - **Deterministic rebuild** — replaying the same tree twice gives
//!   an identical solid.
//! - **Pad volume** — a Pad of a square profile produces a prism of
//!   exactly `area · depth`.
//! - **Pad → Pocket** — pocketing a hole removes exactly the pocketed
//!   volume; the result is a valid closed solid.
//! - **Literal-parameter propagation** — changing a `Value::Literal`
//!   depth and replaying yields the volume the new depth implies.
//! - **Expression-parameter propagation** — a depth bound to a
//!   spreadsheet cell tracks the cell's value across replays.
//! - **Suppression** — a suppressed feature drops out of the rebuild.

use valenx_cad::measure::{is_closed_solid_tol, solid_volume_tol};
use valenx_feature_tree::feature::{Feature, PadParams, PocketParams, Value};
use valenx_feature_tree::{replay, replay_with_spreadsheet, FeatureTree};
use valenx_spreadsheet::{Cell, CellRef, Spreadsheet};

const FINE: f64 = 5.0e-4;

/// A square sketch centred at the origin, side `2·half`.
fn square_sketch(half: f64) -> valenx_sketch::Sketch {
    let mut s = valenx_sketch::Sketch::new();
    let a = s.add_point(-half, -half);
    let b = s.add_point(half, -half);
    let c = s.add_point(half, half);
    let d = s.add_point(-half, half);
    s.add_line(a, b).unwrap();
    s.add_line(b, c).unwrap();
    s.add_line(c, d).unwrap();
    s.add_line(d, a).unwrap();
    s
}

// ===========================================================================
// Deterministic rebuild.
// ===========================================================================

#[test]
fn replaying_a_tree_twice_gives_the_same_solid() {
    // A parametric model must rebuild deterministically — the same
    // tree replayed twice must give an identical solid (same topology,
    // same measured volume).
    let mut tree = FeatureTree::new();
    let s = tree.add_sketch(square_sketch(2.0));
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: s,
            depth: 3.0.into(),
            direction_positive: true,
        }),
        "Base Pad",
    );
    let first = replay(&tree).unwrap().expect("first replay");
    let second = replay(&tree).unwrap().expect("second replay");
    assert_eq!(
        first.faces(),
        second.faces(),
        "rebuild face count must be deterministic"
    );
    let v1 = solid_volume_tol(&first, FINE).unwrap();
    let v2 = solid_volume_tol(&second, FINE).unwrap();
    assert!(
        (v1 - v2).abs() < 1e-9,
        "rebuild volume must be deterministic: {v1} vs {v2}"
    );
}

// ===========================================================================
// Pad — analytic prism volume.
// ===========================================================================

#[test]
fn pad_of_a_square_produces_a_prism_of_exact_volume() {
    // A 4×4 square (half = 2.0) padded by depth 3 → a 4×4×3 prism,
    // volume 48. The replayed solid must measure to exactly that.
    let mut tree = FeatureTree::new();
    let s = tree.add_sketch(square_sketch(2.0));
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: s,
            depth: 3.0.into(),
            direction_positive: true,
        }),
        "Base Pad",
    );
    let solid = replay(&tree).unwrap().expect("pad replay");
    let v = solid_volume_tol(&solid, FINE).unwrap();
    assert!(
        (v - 48.0).abs() < 1e-9,
        "4×4 square padded by 3 should have volume 48, got {v}"
    );
    assert!(
        is_closed_solid_tol(&solid, FINE).unwrap(),
        "a padded prism must be a valid closed solid"
    );
}

// ===========================================================================
// Pad → Pocket — the hole removes its analytic volume.
// ===========================================================================

#[test]
fn pad_then_pocket_removes_the_pocketed_volume() {
    // Pad a 6×6 square (half = 3) by 4 → volume 144.
    // Pocket a 2×2 square (half = 1) blind to depth 2.5 →
    // removes 2·2·2.5 = 10. Result volume = 144 − 10 = 134.
    let mut tree = FeatureTree::new();
    let base = tree.add_sketch(square_sketch(3.0));
    let hole = tree.add_sketch(square_sketch(1.0));
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: base,
            depth: 4.0.into(),
            direction_positive: true,
        }),
        "Base Pad",
    );
    tree.add_feature(
        Feature::Pocket(PocketParams {
            sketch: hole,
            depth: 2.5.into(),
            direction_positive: true,
        }),
        "Blind Pocket",
    );
    let solid = replay(&tree).unwrap().expect("pad+pocket replay");
    let v = solid_volume_tol(&solid, FINE).unwrap();
    assert!(
        (v - 134.0).abs() / 134.0 < 0.005,
        "padded box minus a blind pocket should be ~134, got {v}"
    );
    assert!(
        is_closed_solid_tol(&solid, FINE).unwrap(),
        "the pocketed solid must be a valid closed solid"
    );
}

#[test]
fn pocket_punched_through_removes_a_through_hole_volume() {
    // Pad a 6×6 square by 3 → volume 108. A "through all" pocket is
    // specified by a depth that runs *past* the part (the conventional
    // CAD idiom — a coincident far cap would otherwise stall the
    // boolean kernel). Pocket a 2×2 square with depth 5 > pad depth 3:
    // the cut punches through and removes exactly the block-thickness
    // worth of material, 2·2·3 = 12. Result = 108 − 12 = 96.
    let mut tree = FeatureTree::new();
    let base = tree.add_sketch(square_sketch(3.0));
    let hole = tree.add_sketch(square_sketch(1.0));
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: base,
            depth: 3.0.into(),
            direction_positive: true,
        }),
        "Base Pad",
    );
    tree.add_feature(
        Feature::Pocket(PocketParams {
            sketch: hole,
            depth: 5.0.into(), // through all — runs past the 3-thick block
            direction_positive: true,
        }),
        "Through Pocket",
    );
    let solid = replay(&tree).unwrap().expect("through-pocket replay");
    let v = solid_volume_tol(&solid, FINE).unwrap();
    assert!(
        (v - 96.0).abs() / 96.0 < 0.005,
        "padded box minus a through-hole should be ~96, got {v}"
    );
}

// ===========================================================================
// Parameter propagation — literal depth.
// ===========================================================================

#[test]
fn changing_a_literal_pad_depth_propagates_through_replay() {
    // The defining property of a parametric model: edit a parameter,
    // replay, and the result reflects the new value. A 4×4 square
    // padded by depth d has volume 16·d. Replaying with d = 2, then
    // d = 5, must give volumes 32 and 80.
    let build_and_measure = |depth: f64| -> f64 {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(2.0));
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: depth.into(),
                direction_positive: true,
            }),
            "Pad",
        );
        let solid = replay(&tree).unwrap().expect("replay");
        solid_volume_tol(&solid, FINE).unwrap()
    };
    let v2 = build_and_measure(2.0);
    let v5 = build_and_measure(5.0);
    assert!((v2 - 32.0).abs() < 1e-9, "depth 2 → volume 32, got {v2}");
    assert!((v5 - 80.0).abs() < 1e-9, "depth 5 → volume 80, got {v5}");
    // The relationship is linear: 5/2 the depth → 5/2 the volume.
    assert!(
        (v5 / v2 - 2.5).abs() < 1e-9,
        "volume must scale linearly with depth"
    );
}

#[test]
fn editing_a_feature_in_place_then_replaying_rebuilds_correctly() {
    // Build a tree, replay, then mutate the Pad's depth in the tree's
    // feature list and replay again — the second result must reflect
    // the edit (a re-parametrisation, not a stale cached solid).
    let mut tree = FeatureTree::new();
    let s = tree.add_sketch(square_sketch(2.0)); // 4×4
    let pad_id = tree.add_feature(
        Feature::Pad(PadParams {
            sketch: s,
            depth: 1.0.into(),
            direction_positive: true,
        }),
        "Pad",
    );
    let before = replay(&tree).unwrap().unwrap();
    let v_before = solid_volume_tol(&before, FINE).unwrap();
    assert!((v_before - 16.0).abs() < 1e-9, "depth 1 → volume 16");

    // Edit the Pad's depth in place to 6.0.
    if let Feature::Pad(p) = &mut tree.features[pad_id.0].feature {
        p.depth = 6.0.into();
    } else {
        panic!("feature 0 should be a Pad");
    }
    let after = replay(&tree).unwrap().unwrap();
    let v_after = solid_volume_tol(&after, FINE).unwrap();
    assert!(
        (v_after - 96.0).abs() < 1e-9,
        "after editing depth to 6, volume should be 96, got {v_after}"
    );
}

// ===========================================================================
// Parameter propagation — spreadsheet-expression depth.
// ===========================================================================

#[test]
fn pad_depth_bound_to_a_spreadsheet_cell_tracks_the_cell() {
    // A Pad whose depth is a `Value::Expression` pointing at a
    // spreadsheet cell must rebuild to the volume the *current* cell
    // value implies — change the cell, replay, get the new volume.
    let cell = CellRef::parse("S.A1").unwrap();
    let mut tree = FeatureTree::new();
    let s = tree.add_sketch(square_sketch(2.0)); // 4×4 → area 16
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: s,
            depth: Value::cell(&cell),
            direction_positive: true,
        }),
        "Parametric Pad",
    );

    // Cell = 3 → depth 3 → volume 48.
    let mut ss = Spreadsheet::new();
    ss.add_sheet("S");
    ss.set_cell(&cell, Cell::Number(3.0)).unwrap();
    let solid = replay_with_spreadsheet(&tree, &ss)
        .unwrap()
        .expect("expression replay");
    let v = solid_volume_tol(&solid, FINE).unwrap();
    assert!(
        (v - 48.0).abs() < 1e-9,
        "cell=3 → depth 3 → volume 48, got {v}"
    );

    // Change the cell to 7 → depth 7 → volume 112. Same tree, new
    // spreadsheet state — the parametric link must track it.
    ss.set_cell(&cell, Cell::Number(7.0)).unwrap();
    let solid2 = replay_with_spreadsheet(&tree, &ss)
        .unwrap()
        .expect("expression replay 2");
    let v2 = solid_volume_tol(&solid2, FINE).unwrap();
    assert!(
        (v2 - 112.0).abs() < 1e-9,
        "cell=7 → depth 7 → volume 112, got {v2}"
    );
}

#[test]
fn pad_depth_expression_with_arithmetic_resolves_correctly() {
    // The depth expression can be a formula, not just a cell ref.
    // depth = "S.A1 * 2 + 1" with A1 = 3 → depth 7 → volume 112.
    let mut tree = FeatureTree::new();
    let s = tree.add_sketch(square_sketch(2.0)); // 4×4
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: s,
            depth: Value::Expression("S.A1 * 2 + 1".into()),
            direction_positive: true,
        }),
        "Formula Pad",
    );
    let mut ss = Spreadsheet::new();
    ss.add_sheet("S");
    ss.set_cell(&CellRef::parse("S.A1").unwrap(), Cell::Number(3.0))
        .unwrap();
    let solid = replay_with_spreadsheet(&tree, &ss).unwrap().unwrap();
    let v = solid_volume_tol(&solid, FINE).unwrap();
    assert!(
        (v - 112.0).abs() < 1e-9,
        "depth = 3*2+1 = 7 → volume 112, got {v}"
    );
}

// ===========================================================================
// Suppression.
// ===========================================================================

#[test]
fn suppressing_the_pocket_restores_the_full_padded_volume() {
    // A Pad → Pocket tree. With the Pocket live the volume is reduced;
    // suppress the Pocket and the rebuild must restore the full padded
    // volume — proof that suppression actually drops the feature.
    let mut tree = FeatureTree::new();
    let base = tree.add_sketch(square_sketch(2.0)); // 4×4
    let hole = tree.add_sketch(square_sketch(0.5)); // 1×1
    tree.add_feature(
        Feature::Pad(PadParams {
            sketch: base,
            depth: 3.0.into(),
            direction_positive: true,
        }),
        "Pad",
    );
    let pocket_id = tree.add_feature(
        Feature::Pocket(PocketParams {
            sketch: hole,
            depth: 4.0.into(), // through all — past the 3-thick pad
            direction_positive: true,
        }),
        "Pocket",
    );

    // Pocket live: 4×4×3 = 48 minus the 1×1 through-hole 1×1×3 = 3 → 45.
    let with_pocket = replay(&tree).unwrap().unwrap();
    let v_pocket = solid_volume_tol(&with_pocket, FINE).unwrap();
    assert!(
        (v_pocket - 45.0).abs() / 45.0 < 0.005,
        "padded box with pocket should be ~45, got {v_pocket}"
    );

    // Suppress the Pocket → the rebuild is just the 48-unit Pad.
    tree.set_suppressed(pocket_id, true).unwrap();
    let suppressed = replay(&tree).unwrap().unwrap();
    let v_suppressed = solid_volume_tol(&suppressed, FINE).unwrap();
    assert!(
        (v_suppressed - 48.0).abs() < 1e-9,
        "with the pocket suppressed the volume must be the full 48, got {v_suppressed}"
    );
}
