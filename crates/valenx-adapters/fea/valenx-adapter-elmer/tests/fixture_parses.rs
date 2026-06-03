//! Lock-in test: the bundled `tests/fixtures/minimal.valenx/cases/
//! heat-cube/case.toml` must parse cleanly through the Elmer
//! adapter's `case_input` module.
//!
//! Pre-fix, the fixture was authored against the schema documented
//! in `case_input.rs` but no test asserted that the adapter actually
//! accepted it. A schema tightening (e.g. adding a required field)
//! would silently drift the fixture out of sync without surfacing
//! anywhere. This integration test pins the contract.

use std::path::PathBuf;

use valenx_adapter_elmer::case_input::{BoundaryCondition, ElmerInput, Equation};

fn fixture_case_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR points at .../crates/valenx-adapters/fea/valenx-adapter-elmer
    // — walk up four levels to the workspace root.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
        .join("cases")
        .join("heat-cube")
}

#[test]
fn heat_cube_fixture_parses_through_elmer_case_input() {
    let dir = fixture_case_dir();
    assert!(
        dir.is_dir(),
        "fixture missing — expected {} to exist",
        dir.display()
    );

    let (header, input) = ElmerInput::from_case_dir(&dir).expect("parse heat-cube");

    // Header — solver tag must match what the GUI dispatches on.
    assert_eq!(header.physics, "fea");
    assert_eq!(header.solver, "elmer.heat");

    // Equation — heat conduction.
    assert!(matches!(input.equation, Equation::HeatEquation));

    // Material — aluminium with the documented properties.
    assert_eq!(input.material.name, "aluminium");
    assert!((input.material.density - 2700.0).abs() < 1e-6);
    assert!((input.material.heat_capacity - 900.0).abs() < 1e-6);
    assert!((input.material.heat_conductivity - 237.0).abs() < 1e-6);

    // Two opposing-face Dirichlet boundaries — hot at 373.15 K and
    // cold at 273.15 K.
    assert_eq!(input.boundaries.len(), 2);
    let hot = input
        .boundaries
        .iter()
        .find(|b| matches!(b, BoundaryCondition::Temperature { name, .. } if name == "hot_face"))
        .expect("hot_face boundary");
    let cold = input
        .boundaries
        .iter()
        .find(|b| matches!(b, BoundaryCondition::Temperature { name, .. } if name == "cold_face"))
        .expect("cold_face boundary");
    if let BoundaryCondition::Temperature { value, .. } = hot {
        assert!((*value - 373.15).abs() < 1e-6);
    }
    if let BoundaryCondition::Temperature { value, .. } = cold {
        assert!((*value - 273.15).abs() < 1e-6);
    }

    // Simulation control knobs.
    assert_eq!(input.simulation.max_output_level, 5);
    assert_eq!(input.simulation.steady_state_max_iterations, 50);
    assert!((input.simulation.convergence_tolerance - 1.0e-6).abs() < 1e-12);
}
