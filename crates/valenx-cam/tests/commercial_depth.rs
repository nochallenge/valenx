//! Integration tests for the commercial-depth toolpath pipeline:
//!
//! - `engagement` + `op::adaptive_constant_engagement` — modern HSM
//!   constant-engagement adaptive clearing.
//! - `arcfit` — G2/G3 fitting on a rounded path.
//! - `feedrate` — three-pass feedrate optimization with lookahead.
//! - `collision` — continuous swept CCD.
//!
//! Each test exercises one feature end-to-end against the public
//! crate API and verifies the published-quality invariants.

use nalgebra::Vector3;
use valenx_cam::{
    arcfit::{fit_arcs, ArcDir, ArcFitParams},
    collision::{
        continuous_collision_check, CollisionBody, CollisionSetup, ContinuousCollisionParams,
        Holder, SetupPartKind,
    },
    feedrate::{optimize, FeedrateParams},
    op::adaptive_constant_engagement::{generate as adaptive_ce, AdaptiveConstantEngagementParams},
    stock::Stock,
    tool::{Tool, ToolKind},
    toolpath::{Move, MoveKind, Toolpath},
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

fn cube(size: f64) -> Mesh {
    let s = size * 0.5;
    let nodes = vec![
        Vector3::new(-s, -s, -s),
        Vector3::new(s, -s, -s),
        Vector3::new(s, s, -s),
        Vector3::new(-s, s, -s),
        Vector3::new(-s, -s, s),
        Vector3::new(s, -s, s),
        Vector3::new(s, s, s),
        Vector3::new(-s, s, s),
    ];
    let conn: Vec<u32> = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5,
    ];
    let mut mesh = Mesh::new("cube");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh
}

#[test]
fn full_pipeline_adaptive_then_arcfit_then_feedrate() {
    // 1. Generate constant-engagement adaptive on a 40x40 cube pocket.
    let mesh = cube(40.0);
    let stock = Stock::new(
        Vector3::new(-25.0, -25.0, -20.0),
        Vector3::new(50.0, 50.0, 40.0),
        "alu",
    )
    .unwrap();
    let tool = Tool::new(1, "EM4", ToolKind::EndMill, 4.0, 25.0, 2, "carbide").unwrap();
    let params = AdaptiveConstantEngagementParams {
        step_down: 4.0,
        depth: 4.0,
        max_engagement_rad: 0.5,
        cell_size_mm: 0.6,
        engagement_samples: 32,
        ..Default::default()
    };
    let (tp, report) = adaptive_ce(&stock, &mesh, &params, &tool).unwrap();
    assert!(report.n_cut_moves > 50);
    // The engagement bound should hold (within sampling bucket).
    let bucket = std::f64::consts::TAU / (params.engagement_samples as f64);
    assert!(report.max_engagement_rad <= params.max_engagement_rad + 2.0 * bucket);

    // 2. Run the arc-fit pass over the resulting toolpath. The
    // adaptive path's rollovers are circular trochoidal loops —
    // good arc-fit candidates.
    let arc_params = ArcFitParams {
        chord_tol_mm: 0.05,
        ..Default::default()
    };
    let (tp_arcs, arc_report) = fit_arcs(&tp, &arc_params);
    // Some arcs should have been emitted because the rollovers are
    // perfect circles.
    assert!(
        arc_report.arcs_emitted > 0,
        "expected ≥ 1 arc from rollover loops, got {}",
        arc_report.arcs_emitted
    );
    // Output should be no larger than input.
    assert!(tp_arcs.len() <= tp.len());

    // 3. Run the feedrate optimizer. Some moves should be reduced
    // because rollovers introduce corners (G1→G2 junctions).
    let feed_params = FeedrateParams::default();
    let (tp_final, feed_report) = optimize(&tp_arcs, &feed_params);
    assert_eq!(tp_final.len(), tp_arcs.len());
    assert!(feed_report.reduced_moves > 0 || feed_report.centripetal_clamps > 0);
}

#[test]
fn arc_fit_reduces_circular_pocket_path_substantially() {
    // Hand-build a 64-segment circle path (no Z change, all Cut moves).
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(10.0, 0.0, 0.0),
        0.0,
    ));
    let r = 10.0;
    for k in 1..=64 {
        let t = (k as f64) * std::f64::consts::TAU / 64.0;
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(r * t.cos(), r * t.sin(), 0.0),
            500.0,
        ));
    }
    let before = tp.len();
    let params = ArcFitParams::default();
    let (out, report) = fit_arcs(&tp, &params);
    // Should reduce moves substantially (Mastercam typically gets
    // > 80 % reduction on a clean circle).
    let reduction = (before - out.len()) as f64 / (before as f64);
    assert!(
        reduction > 0.5,
        "expected > 50 % reduction, got {:.0} % ({} → {})",
        reduction * 100.0,
        before,
        out.len()
    );
    assert!(report.arcs_emitted >= 1);
}

#[test]
fn feedrate_optimizer_handles_arc_centripetal_correctly() {
    // Build a tight-radius arc path and confirm feed is bounded by
    // sqrt(a · r) per the centripetal-acceleration formula.
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Cut,
        Vector3::new(0.0, 0.0, 0.0),
        5000.0,
    ));
    tp.push(Move {
        kind: MoveKind::Arc {
            centre_xy: nalgebra::Vector2::new(1.0, 0.0),
            dir: ArcDir::Counterclockwise,
        },
        position: Vector3::new(2.0, 0.0, 0.0),
        feed: 5000.0,
    });
    // Limit centripetal acceleration so the bound is observable.
    let params = FeedrateParams {
        a_centripetal_max_mm_per_min2: 1_000_000.0, // 1e6
        ..Default::default()
    };
    let (out, report) = optimize(&tp, &params);
    // r = 1.0, v_max = sqrt(1e6 * 1.0) = 1000 mm/min.
    let v = out.moves[1].feed;
    assert!(
        v <= 1010.0,
        "centripetal bound failed: feed {v} > 1010 mm/min (r=1, a=1e6)",
    );
    assert_eq!(report.centripetal_clamps, 1);
}

#[test]
fn continuous_collision_catches_grazing_rapid_fixture_v1_misses() {
    let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
    let mut tp = Toolpath::new();
    // Endpoints clear the fixture, midpoint slices through it.
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(-20.0, 0.0, 5.0),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(20.0, 0.0, 5.0),
        0.0,
    ));
    let mut setup = CollisionSetup::new();
    // The clamp post must span the tool's z = 5 so the flute -- which
    // extends UP from the tool tip, `[tip_z, tip_z + flute_len]`, per the
    // #579 convention fix -- grazes it as the rapid sweeps through. (A
    // post topping out at z = 3, below the tip, sits entirely under the
    // flute and is never touched.)
    setup.push_fixture(
        Vector3::new(-2.5, -10.0, 0.0),
        Vector3::new(2.5, 10.0, 10.0),
        "clamp_post",
    );
    let hits = continuous_collision_check(
        &tp,
        &tool,
        &Holder::empty(),
        &setup,
        &ContinuousCollisionParams::default(),
    );
    assert!(
        !hits.is_empty(),
        "continuous CCD failed to catch grazing rapid"
    );
    assert_eq!(hits[0].part_kind, SetupPartKind::Fixture);
    assert_eq!(hits[0].body, CollisionBody::Flute);
}

#[test]
fn continuous_collision_reports_clean_path_clean() {
    let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
    let mut tp = Toolpath::new();
    // Path well above the fixture in Z.
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(-20.0, 0.0, 100.0),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(20.0, 0.0, 100.0),
        0.0,
    ));
    let mut setup = CollisionSetup::new();
    setup.push_fixture(
        Vector3::new(-2.5, -10.0, 0.0),
        Vector3::new(2.5, 10.0, 3.0),
        "clamp_post",
    );
    let hits = continuous_collision_check(
        &tp,
        &tool,
        &Holder::empty(),
        &setup,
        &ContinuousCollisionParams::default(),
    );
    assert!(hits.is_empty(), "clear path reported {hits:?}");
}

#[test]
fn holder_collision_detected_separately_from_flute() {
    // Short tool with a wide holder that protrudes above the flute
    // and crashes a fixture above the workpiece.
    let tool = Tool::new(1, "EM3", ToolKind::EndMill, 3.0, 8.0, 2, "").unwrap();
    let holder = Holder::cylinder_shank(12.0, 30.0);
    let mut tp = Toolpath::new();
    tp.push(Move::new(MoveKind::Cut, Vector3::new(0.0, 0.0, 0.0), 500.0));
    tp.push(Move::new(
        MoveKind::Cut,
        Vector3::new(20.0, 0.0, 0.0),
        500.0,
    ));
    // Fixture sits *above* the tool tip (z=12..18 — within the holder
    // span of 8..38).
    let mut setup = CollisionSetup::new();
    setup.push_fixture(
        Vector3::new(8.0, -10.0, 12.0),
        Vector3::new(12.0, 10.0, 18.0),
        "ledge",
    );
    let hits = continuous_collision_check(
        &tp,
        &tool,
        &holder,
        &setup,
        &ContinuousCollisionParams::default(),
    );
    let any_holder = hits.iter().any(|h| h.body == CollisionBody::Holder);
    assert!(any_holder, "expected a Holder collision in {hits:?}");
}
