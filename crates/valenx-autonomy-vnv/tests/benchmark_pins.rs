//! Benchmark-pin tests for the V&V framework.
//!
//! These pin the framework's *logic* against hand-built ground truth (per the
//! crate's honesty statement — they validate that requirements/coverage/sweep
//! compute correctly, not that any real autonomous system is safe):
//!
//! 1. A trace that violates `MinClearance` FAILS with the correct margin; a safe
//!    one PASSES (with the correct positive margin).
//! 2. A scenario engineered for a guaranteed collision FAILS `NoCollision`; one
//!    engineered safe PASSES.
//! 3. Coverage is exactly 100% when the suite spans every grid cell, and an
//!    exact fraction (< 100%) with a deliberate gap.
//! 4. A deterministic (grid and seeded-MC) sweep yields a known pass-rate.
//!
//! Plus fail-loud checks on bad configuration (empty suite, NaN params,
//! trace/requirement mismatch).

use std::collections::BTreeMap;

use nalgebra::Vector3;
use valenx_autonomy_vnv::{
    evaluate, grid_suite_auto, monte_carlo_suite, parameter_coverage, requirement_coverage,
    run_scenario, run_suite, Aabb, CommandSeq, Lidar, LidarConfig, ParamGrid, Requirement,
    RequirementSet, SampleAxis, Scenario, ScenarioSuite, SensorSet, VnvError,
};
use valenx_sensors::{Command, Plane, Scene, Sphere, VehicleState};

fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
    Vector3::new(x, y, z)
}

/// A vehicle coasting along +x at a constant velocity. With `Command::coast()`
/// (zero accel) and semi-implicit Euler, position.x after `n` ticks of `dt` is
/// exactly `vx · n · dt` — so positions are exact, no integration slop.
fn coasting_scenario(name: &str, vx: f64, dt: f64, steps: usize, scene: Scene) -> Scenario {
    let state = VehicleState {
        velocity: v(vx, 0.0, 0.0),
        ..Default::default()
    };
    Scenario::new(
        name,
        state,
        scene,
        CommandSeq::Constant {
            command: Command::coast(),
            dt,
            steps,
        },
    )
}

// ---------------------------------------------------------------------------
// Benchmark 1 — MinClearance: violating trace fails with correct margin,
// safe trace passes.
// ---------------------------------------------------------------------------

#[test]
fn min_clearance_safe_trace_passes_with_correct_margin() {
    // Stationary vehicle at the origin; a sphere (centre x=5, r=1) ⇒ the nearest
    // surface point is at x=4, so the clearance is exactly 4 m at every tick.
    let mut scene = Scene::new();
    scene.push_sphere(Sphere::new(v(5.0, 0.0, 0.0), 1.0).unwrap());
    let scenario = coasting_scenario("clear", 0.0, 0.1, 5, scene); // vx=0 ⇒ stays at origin

    let trace = run_scenario(&scenario).unwrap();
    // Sanity: it really did not move.
    assert!(trace.frames.last().unwrap().state.position.norm() < 1e-12);

    let req = Requirement::MinClearance { d: 2.0 };
    let outcome = req.evaluate(&trace).unwrap();
    assert!(outcome.pass, "4 m clearance must satisfy d=2");
    // margin = clearance(4) − d(2) = 2, exactly.
    assert!(
        (outcome.margin - 2.0).abs() < 1e-12,
        "margin = {}",
        outcome.margin
    );
    assert!(outcome.exercised);
}

#[test]
fn min_clearance_violating_trace_fails_with_correct_margin() {
    let mut scene = Scene::new();
    scene.push_sphere(Sphere::new(v(5.0, 0.0, 0.0), 1.0).unwrap());
    let scenario = coasting_scenario("too-close", 0.0, 0.1, 5, scene);
    let trace = run_scenario(&scenario).unwrap();

    // Require 5 m clearance, but only 4 m is available ⇒ violated by 1 m.
    let req = Requirement::MinClearance { d: 5.0 };
    let outcome = req.evaluate(&trace).unwrap();
    assert!(!outcome.pass, "4 m clearance cannot satisfy d=5");
    assert!(
        (outcome.margin - (-1.0)).abs() < 1e-12,
        "margin = {}",
        outcome.margin
    );
}

// ---------------------------------------------------------------------------
// Benchmark 2 — NoCollision: guaranteed collision fails, safe passes.
// ---------------------------------------------------------------------------

#[test]
fn engineered_collision_fails_no_collision() {
    // Sphere centre x=10, r=1 (surface spans x∈[9,11]). Coast at vx=2, dt=0.5 ⇒
    // tick positions x = 1,2,…,12. Ticks at x=9,10,11 lie on/inside the sphere,
    // where the surface distance clamps to 0, so NoCollision(radius=0.5) is
    // violated by 0.5.
    let mut scene = Scene::new();
    scene.push_sphere(Sphere::new(v(10.0, 0.0, 0.0), 1.0).unwrap());
    let scenario = coasting_scenario("drive-into-it", 2.0, 0.5, 12, scene);
    let trace = run_scenario(&scenario).unwrap();
    // Confirm a tick actually reaches x=10 (the centre line).
    assert!(trace
        .frames
        .iter()
        .any(|f| (f.state.position.x - 10.0).abs() < 1e-9));

    let req = Requirement::NoCollision { radius: 0.5 };
    let outcome = req.evaluate(&trace).unwrap();
    assert!(
        !outcome.pass,
        "passing through the sphere must fail NoCollision"
    );
    assert!(
        (outcome.margin - (-0.5)).abs() < 1e-9,
        "margin = {}",
        outcome.margin
    );
}

#[test]
fn engineered_safe_passes_no_collision() {
    // Same motion, but the sphere is 100 m away ⇒ never close ⇒ passes.
    let mut scene = Scene::new();
    scene.push_sphere(Sphere::new(v(100.0, 0.0, 0.0), 1.0).unwrap());
    let scenario = coasting_scenario("safe-pass", 2.0, 0.5, 12, scene);
    let trace = run_scenario(&scenario).unwrap();

    let req = Requirement::NoCollision { radius: 0.5 };
    let outcome = req.evaluate(&trace).unwrap();
    assert!(outcome.pass, "never close to the far sphere");
    assert!(outcome.margin > 0.0);
}

// ---------------------------------------------------------------------------
// StayInBounds + DetectByTime spot checks (exact margins).
// ---------------------------------------------------------------------------

#[test]
fn stay_in_bounds_pass_and_fail_have_exact_margins() {
    // Coast at vx=1, dt=1, 5 ticks ⇒ x = 1,2,3,4,5.
    let scenario = coasting_scenario("box", 1.0, 1.0, 5, Scene::new());
    let trace = run_scenario(&scenario).unwrap();

    // Box x∈[0,10], y,z∈[−1,1]. At x=5 the worst-axis slack is min(5−0,10−5,
    // 1−0,1−0) = 1 (the y/z half-width), and the tightest tick is x=1 ⇒ slack
    // min(1, 9, 1, 1) = 1. So the min signed clearance over the run is 1.
    let bounds = Aabb::new(v(0.0, -1.0, -1.0), v(10.0, 1.0, 1.0)).unwrap();
    let inside = Requirement::StayInBounds { bounds }
        .evaluate(&trace)
        .unwrap();
    assert!(inside.pass);
    assert!(
        (inside.margin - 1.0).abs() < 1e-12,
        "margin = {}",
        inside.margin
    );

    // Tight box x∈[0,3]: the vehicle reaches x=5, overshooting the far face by 2.
    let tight = Aabb::new(v(0.0, -10.0, -10.0), v(3.0, 10.0, 10.0)).unwrap();
    let out = Requirement::StayInBounds { bounds: tight }
        .evaluate(&trace)
        .unwrap();
    assert!(!out.pass);
    assert!(
        (out.margin - (-2.0)).abs() < 1e-12,
        "margin = {}",
        out.margin
    );
}

#[test]
fn detect_by_time_detects_known_target() {
    // Wall at x=5 (facing −x); a single forward LiDAR beam from the origin reads
    // 5 m, so the beam endpoint is exactly (5,0,0).
    let mut scene = Scene::new();
    scene.push_plane(Plane::new(v(5.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());
    let lidar = Lidar::new(
        LidarConfig {
            azimuth_steps: 1,
            elevation_steps: 1,
            h_fov: 0.0,
            v_fov: 0.0,
            min_range: 0.0,
            max_range: 100.0,
            range_noise_std: 0.0,
        },
        0,
    )
    .unwrap();
    let scenario = Scenario::new(
        "detect",
        VehicleState::default(),
        scene,
        CommandSeq::Constant {
            command: Command::coast(),
            dt: 0.1,
            steps: 5,
        },
    )
    .with_sensors(SensorSet::none().with_lidar(lidar));

    let trace = run_scenario(&scenario).unwrap();
    let req = Requirement::DetectByTime {
        target: v(5.0, 0.0, 0.0),
        t_max: 1.0,
        tol: 0.05,
    };
    let outcome = req.evaluate(&trace).unwrap();
    assert!(outcome.pass, "the wall is detected immediately");
    // Detected at the first tick t=0.1 ⇒ margin = t_max − 0.1 = 0.9.
    assert!(
        (outcome.margin - 0.9).abs() < 1e-9,
        "margin = {}",
        outcome.margin
    );
}

#[test]
fn detect_by_time_missing_target_fails() {
    // Same beam, but the target is somewhere the beam never reaches.
    let mut scene = Scene::new();
    scene.push_plane(Plane::new(v(5.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());
    let lidar = Lidar::new(
        LidarConfig {
            azimuth_steps: 1,
            elevation_steps: 1,
            h_fov: 0.0,
            v_fov: 0.0,
            min_range: 0.0,
            max_range: 100.0,
            range_noise_std: 0.0,
        },
        0,
    )
    .unwrap();
    let scenario = Scenario::new(
        "no-detect",
        VehicleState::default(),
        scene,
        CommandSeq::Constant {
            command: Command::coast(),
            dt: 0.1,
            steps: 5,
        },
    )
    .with_sensors(SensorSet::none().with_lidar(lidar));

    let trace = run_scenario(&scenario).unwrap();
    let req = Requirement::DetectByTime {
        target: v(0.0, 50.0, 0.0), // off to the side, never seen
        t_max: 1.0,
        tol: 0.05,
    };
    let outcome = req.evaluate(&trace).unwrap();
    assert!(!outcome.pass, "an unseen target is not detected");
    assert!(outcome.margin < 0.0, "failure margin must be negative");
}

// ---------------------------------------------------------------------------
// Benchmark 3 — coverage: 100% when the suite spans every cell; an exact
// fraction with a deliberate gap.
// ---------------------------------------------------------------------------

fn coverage_grid() -> ParamGrid {
    // 3 × 2 = 6 cells.
    ParamGrid::new()
        .with_axis("obstacle_x", vec![3.0, 5.0, 100.0])
        .with_axis("speed", vec![1.0, 2.0])
}

/// Build a scenario for a grid cell: a sphere at `obstacle_x`, vehicle coasting
/// at `speed`. (Behaviour is unimportant for the coverage count — only the
/// params tags matter — but a real scenario keeps the test honest end-to-end.)
fn cell_scenario(cell: &BTreeMap<String, f64>, i: usize) -> Result<Scenario, VnvError> {
    let ox = cell["obstacle_x"];
    let speed = cell["speed"];
    let mut scene = Scene::new();
    scene.push_sphere(Sphere::new(v(ox, 0.0, 0.0), 1.0).unwrap());
    Ok(coasting_scenario(
        &format!("cell-{i}"),
        speed,
        0.5,
        4,
        scene,
    ))
}

#[test]
fn coverage_is_100_percent_when_suite_spans_every_cell() {
    let grid = coverage_grid();
    let suite = grid_suite_auto("full", &grid, cell_scenario).unwrap();
    assert_eq!(suite.len(), 6, "full grid suite has one scenario per cell");

    let cov = parameter_coverage(&suite, &grid).unwrap();
    assert_eq!(cov.total_cells, 6);
    assert_eq!(cov.covered_cells, 6);
    assert!(cov.is_complete());
    assert!((cov.fraction() - 1.0).abs() < 1e-12);
}

#[test]
fn coverage_is_exact_fraction_with_a_deliberate_gap() {
    let grid = coverage_grid();
    let full = grid_suite_auto("full", &grid, cell_scenario).unwrap();

    // Drop one scenario ⇒ 5 of 6 cells covered ⇒ exactly 5/6.
    let mut scenarios = full.scenarios;
    scenarios.pop();
    let gapped = ScenarioSuite::new("gapped", scenarios);

    let cov = parameter_coverage(&gapped, &grid).unwrap();
    assert_eq!(cov.total_cells, 6);
    assert_eq!(cov.covered_cells, 5);
    assert!(!cov.is_complete());
    assert!(
        (cov.fraction() - 5.0 / 6.0).abs() < 1e-12,
        "fraction = {}",
        cov.fraction()
    );
}

#[test]
fn off_grid_scenarios_cover_no_cell() {
    let grid = coverage_grid();
    // A scenario whose param values are not on any axis ⇒ covers nothing.
    let mut scene = Scene::new();
    scene.push_sphere(Sphere::new(v(7.0, 0.0, 0.0), 1.0).unwrap());
    let off = coasting_scenario("off", 1.0, 0.5, 4, scene)
        .with_param("obstacle_x", 7.0) // 7.0 is not in {3,5,100}
        .with_param("speed", 1.0);
    let suite = ScenarioSuite::new("off-grid", vec![off]);

    let cov = parameter_coverage(&suite, &grid).unwrap();
    assert_eq!(cov.covered_cells, 0, "off-grid value snaps to no cell");
}

#[test]
fn requirement_coverage_tracks_exercised_and_triggered() {
    // Two runs against two requirements: one requirement that always passes
    // (never triggered) and one that fails on a close approach (triggered).
    let mut close = Scene::new();
    close.push_sphere(Sphere::new(v(5.0, 0.0, 0.0), 1.0).unwrap());
    let close_scenario = coasting_scenario("close", 0.0, 0.1, 3, close); // clearance 4

    let mut far = Scene::new();
    far.push_sphere(Sphere::new(v(100.0, 0.0, 0.0), 1.0).unwrap());
    let far_scenario = coasting_scenario("far", 0.0, 0.1, 3, far);

    let reqs = RequirementSet::new(vec![
        Requirement::MinClearance { d: 2.0 }, // 4≥2 ⇒ always passes here
        Requirement::MinClearance { d: 10.0 }, // 4<10 close ⇒ fails ⇒ triggered
    ]);

    let r1 = evaluate(&reqs, &run_scenario(&close_scenario).unwrap()).unwrap();
    let r2 = evaluate(&reqs, &run_scenario(&far_scenario).unwrap()).unwrap();
    let cov = requirement_coverage(&[r1, r2]);

    assert_eq!(cov.num_requirements(), 2);
    // Both requirements were exercised in every run.
    assert!((cov.exercised_fraction() - 1.0).abs() < 1e-12);
    // Exactly one of the two was ever driven to a failure (the d=10 one, on the
    // close run). d=2 never failed (far run has clearance 99).
    assert!(cov.triggered["MinClearance(d=10)"]);
    assert!(!cov.triggered["MinClearance(d=2)"]);
    assert!((cov.triggered_fraction() - 0.5).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// Benchmark 4 — deterministic sweep yields a known pass-rate.
// ---------------------------------------------------------------------------

#[test]
fn grid_sweep_has_known_pass_rate_and_worst_margin() {
    // Grid: obstacle_x ∈ {3, 5, 100}, speed ∈ {1,2} (6 cells). The vehicle is
    // stationary (the cell builder uses `speed` only as a tag here? no — it
    // coasts, but for the clearance check we make it stationary). To get an
    // exactly-known pass-rate we use a stationary vehicle and require ≥ 3 m
    // clearance: clearance = obstacle_x − 1 (sphere r=1). So:
    //   obstacle_x=3  ⇒ clearance 2  < 3 ⇒ FAIL (2 cells: speed 1 & 2)
    //   obstacle_x=5  ⇒ clearance 4  ≥ 3 ⇒ PASS (2 cells)
    //   obstacle_x=100⇒ clearance 99 ≥ 3 ⇒ PASS (2 cells)
    // ⇒ 4 of 6 pass ⇒ pass-rate 2/3. Worst margin = 2−3 = −1 (the x=3 cells).
    let grid = coverage_grid();
    let suite = grid_suite_auto("clearance-grid", &grid, |cell, i| {
        let ox = cell["obstacle_x"];
        let mut scene = Scene::new();
        scene.push_sphere(Sphere::new(v(ox, 0.0, 0.0), 1.0).unwrap());
        // vx = 0 ⇒ stationary, so clearance is exactly obstacle_x − 1.
        Ok(coasting_scenario(&format!("cell-{i}"), 0.0, 0.5, 3, scene))
    })
    .unwrap();

    let reqs = RequirementSet::new(vec![Requirement::MinClearance { d: 3.0 }]);
    let result = run_suite(&suite, &reqs, Some(&grid)).unwrap();

    assert_eq!(result.total_runs, 6);
    assert_eq!(result.passed_runs, 4, "x=5 and x=100 pass for both speeds");
    assert!(
        (result.pass_rate() - 2.0 / 3.0).abs() < 1e-12,
        "rate = {}",
        result.pass_rate()
    );
    // Worst margin across the whole sweep = the x=3 clearance shortfall = −1.
    assert!(
        (result.worst_margin().unwrap() - (-1.0)).abs() < 1e-12,
        "worst = {:?}",
        result.worst_margin()
    );
    // Coverage attached and complete (the suite IS the full grid).
    let cov = result.parameter_coverage.unwrap();
    assert!(cov.is_complete());
    assert!((cov.fraction() - 1.0).abs() < 1e-12);
}

#[test]
fn seeded_monte_carlo_sweep_is_reproducible_and_has_known_pass_rate() {
    // Sample obstacle_x ∈ [50, 100] (always far) for a stationary vehicle; the
    // clearance is obstacle_x − 1 ≥ 49, so MinClearance(d=3) ALWAYS passes
    // regardless of the (random) sample ⇒ pass-rate is exactly 1.0 for ANY seed.
    // This pins "a seeded sweep yields a known pass-rate" without coupling the
    // assertion to the exact PRNG stream.
    let mut axes = BTreeMap::new();
    axes.insert(
        "obstacle_x".to_string(),
        SampleAxis::new(50.0, 100.0).unwrap(),
    );

    let build = |sample: &BTreeMap<String, f64>, i: usize| {
        let ox = sample["obstacle_x"];
        let mut scene = Scene::new();
        scene.push_sphere(Sphere::new(v(ox, 0.0, 0.0), 1.0).unwrap());
        Ok(coasting_scenario(&format!("mc-{i}"), 0.0, 0.5, 3, scene))
    };

    let suite_a = monte_carlo_suite("mc", &axes, 32, 0xC0FFEE, build).unwrap();
    let suite_b = monte_carlo_suite("mc", &axes, 32, 0xC0FFEE, build).unwrap();
    // Same seed ⇒ identical sampled params (reproducible suite).
    let params_a: Vec<_> = suite_a.scenarios.iter().map(|s| s.params.clone()).collect();
    let params_b: Vec<_> = suite_b.scenarios.iter().map(|s| s.params.clone()).collect();
    assert_eq!(params_a, params_b, "same seed ⇒ identical MC suite");

    let reqs = RequirementSet::new(vec![Requirement::MinClearance { d: 3.0 }]);
    let result = run_suite(&suite_a, &reqs, None).unwrap();
    assert_eq!(result.total_runs, 32);
    assert!(result.all_passed());
    assert!((result.pass_rate() - 1.0).abs() < 1e-12);

    // A different seed still gives the all-pass rate (the requirement can't fail
    // on this envelope) but a DIFFERENT sample set — confirming the seed drives
    // the draw.
    let suite_c = monte_carlo_suite("mc", &axes, 32, 0x1234_5678, build).unwrap();
    let params_c: Vec<_> = suite_c.scenarios.iter().map(|s| s.params.clone()).collect();
    assert_ne!(params_a, params_c, "different seed ⇒ different MC suite");
    let result_c = run_suite(&suite_c, &reqs, None).unwrap();
    assert!((result_c.pass_rate() - 1.0).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// Fail-loud on bad configuration.
// ---------------------------------------------------------------------------

#[test]
fn empty_suite_fails_loud() {
    let suite = ScenarioSuite::new("empty", vec![]);
    let reqs = RequirementSet::new(vec![Requirement::NoCollision { radius: 0.5 }]);
    let err = run_suite(&suite, &reqs, None).unwrap_err();
    assert!(matches!(err, VnvError::InvalidConfig(_)), "got {err:?}");
}

#[test]
fn nan_scenario_param_fails_loud() {
    let scenario = coasting_scenario("nan", 0.0, 0.1, 3, Scene::new()).with_param("bad", f64::NAN);
    let err = scenario.validate().unwrap_err();
    assert!(matches!(err, VnvError::NonFinite(_)), "got {err:?}");
    // And the runner refuses to run it.
    assert!(run_scenario(&scenario).is_err());
}

#[test]
fn nan_initial_state_fails_loud() {
    let state = VehicleState {
        position: v(f64::NAN, 0.0, 0.0),
        ..Default::default()
    };
    let scenario = Scenario::new(
        "nan-state",
        state,
        Scene::new(),
        CommandSeq::Constant {
            command: Command::coast(),
            dt: 0.1,
            steps: 3,
        },
    );
    assert!(matches!(scenario.validate(), Err(VnvError::NonFinite(_))));
}

#[test]
fn empty_command_sequence_fails_loud() {
    let scenario = Scenario::new(
        "empty-cmd",
        VehicleState::default(),
        Scene::new(),
        CommandSeq::Explicit {
            commands: vec![],
            dt: 0.1,
        },
    );
    assert!(matches!(
        scenario.validate(),
        Err(VnvError::InvalidConfig(_))
    ));
}

#[test]
fn nonpositive_dt_fails_loud() {
    let scenario = Scenario::new(
        "bad-dt",
        VehicleState::default(),
        Scene::new(),
        CommandSeq::Constant {
            command: Command::coast(),
            dt: 0.0,
            steps: 3,
        },
    );
    assert!(matches!(
        scenario.validate(),
        Err(VnvError::InvalidConfig(_))
    ));
}

#[test]
fn detect_requirement_over_lidarless_trace_is_a_mismatch() {
    // A trace with NO LiDAR (no sensors) ⇒ DetectByTime cannot be answered ⇒
    // RequirementMismatch (not a silent fail).
    let scenario = coasting_scenario("no-lidar", 0.0, 0.1, 3, Scene::new());
    let trace = run_scenario(&scenario).unwrap();
    assert!(!trace.has_lidar());
    let req = Requirement::DetectByTime {
        target: v(5.0, 0.0, 0.0),
        t_max: 1.0,
        tol: 0.1,
    };
    let err = req.evaluate(&trace).unwrap_err();
    assert!(
        matches!(err, VnvError::RequirementMismatch(_)),
        "got {err:?}"
    );
}

#[test]
fn requirement_over_empty_trace_is_a_mismatch() {
    // Construct an empty trace directly (0 frames) ⇒ any requirement is a
    // mismatch.
    let trace = valenx_autonomy_vnv::Trace {
        scenario: "empty".to_string(),
        initial_state: VehicleState::default(),
        scene: Scene::new(),
        frames: vec![],
    };
    let req = Requirement::MinClearance { d: 1.0 };
    let err = req.evaluate(&trace).unwrap_err();
    assert!(
        matches!(err, VnvError::RequirementMismatch(_)),
        "got {err:?}"
    );
}

#[test]
fn bad_requirement_params_fail_loud() {
    // Non-finite / out-of-range requirement params are rejected at validate.
    assert!(Requirement::MinClearance { d: f64::NAN }
        .validate()
        .is_err());
    assert!(Requirement::NoCollision { radius: -1.0 }
        .validate()
        .is_err());
    assert!(Requirement::DetectByTime {
        target: v(0.0, 0.0, 0.0),
        t_max: -1.0,
        tol: 0.1
    }
    .validate()
    .is_err());
    assert!(Requirement::DetectByTime {
        target: v(0.0, 0.0, 0.0),
        t_max: 1.0,
        tol: 0.0 // tol must be > 0
    }
    .validate()
    .is_err());
}

#[test]
fn zero_length_grid_axis_fails_loud() {
    let grid = ParamGrid::new().with_axis("x", vec![]);
    assert!(matches!(grid.validate(), Err(VnvError::InvalidConfig(_))));
}

#[test]
fn nan_grid_value_fails_loud() {
    let grid = ParamGrid::new().with_axis("x", vec![1.0, f64::NAN]);
    assert!(matches!(grid.validate(), Err(VnvError::NonFinite(_))));
}

#[test]
fn zero_count_monte_carlo_fails_loud() {
    let axes: BTreeMap<String, SampleAxis> = BTreeMap::new();
    let err = monte_carlo_suite("mc", &axes, 0, 1, |_, i| {
        Ok(coasting_scenario(
            &format!("x{i}"),
            0.0,
            0.1,
            1,
            Scene::new(),
        ))
    })
    .unwrap_err();
    assert!(matches!(err, VnvError::InvalidConfig(_)), "got {err:?}");
}

#[test]
fn run_to_trace_is_byte_identical_for_a_scenario() {
    // The whole pipeline is reproducible: the same scenario yields an identical
    // trace (sensor noise is seeded). Use a noisy LiDAR to make this non-trivial.
    let mut scene = Scene::new();
    scene.push_plane(Plane::new(v(5.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)).unwrap());
    let lidar = Lidar::new(
        LidarConfig {
            azimuth_steps: 1,
            elevation_steps: 1,
            h_fov: 0.0,
            v_fov: 0.0,
            min_range: 0.0,
            max_range: 100.0,
            range_noise_std: 0.05, // noisy ⇒ exercises the seeded PRNG
        },
        7,
    )
    .unwrap();
    let scenario = Scenario::new(
        "repro",
        VehicleState::default(),
        scene,
        CommandSeq::Constant {
            command: Command {
                accel_body: v(0.3, 0.0, 0.0),
                angular_rate_body: v(0.0, 0.0, 0.01),
            },
            dt: 0.05,
            steps: 20,
        },
    )
    .with_sensors(SensorSet::none().with_lidar(lidar));

    let a = run_scenario(&scenario).unwrap();
    let b = run_scenario(&scenario).unwrap();
    assert_eq!(a, b, "same scenario ⇒ byte-identical trace");
}
